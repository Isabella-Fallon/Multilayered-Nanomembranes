%% clean_artifacts.m
%
% Removes high-amplitude artifacts from NeuroExplorer continuous bandpass
% variables (e.g. raw_1_bandpass1_10hz, raw_1_bandpass10_30hz,
% raw_1_bandpass300_3000hz) using a robust median/MAD threshold, then
% replaces flagged samples (plus padding) with linear interpolation so
% downstream spike detection isn't skewed.
%
% Output: <name>_clean, an N-by-2 matrix [voltage, timestamp(s)], one per channel.
%
% HOW TO USE
%   1. Bring the continuous variables into the MATLAB workspace as plain
%      numeric vectors named like raw_1_bandpass1_10hz.
%   2. Check the SETTINGS block (sampling rate per band, threshold, padding).
%   3. Run. A QC plot of the first channel of each band is drawn at the end.

%% ---------------- SETTINGS ----------------
% One row per band: tag = text in the variable name, fs = sampling rate of
% that variable (Hz), n_mad = threshold in scaled MADs from the median.
bands(1) = struct('tag', 'bandpass1_10hz',     'fs', 1500,  'n_mad', 20);
bands(2) = struct('tag', 'bandpass10_30hz',    'fs', 1500,  'n_mad', 20);
bands(3) = struct('tag', 'bandpass300_3000hz', 'fs', 30000, 'n_mad', 20);

pad_ms      = 2;       % ms blanked on each side of every flagged sample
name_prefix = 'raw';   % only variables starting with this are considered
t0          = 0;       % start time (s) of the recording; set if it isn't 0
clear_raw   = false;   % true = delete each raw variable from the workspace once its
                       % _clean version is saved (frees ~1/3 of the memory; re-import
                       % from NeuroExplorer if you need the raw trace again)
% --------------------------------------------

all_vars = evalin('base', 'who');
results  = struct();
qc_names = cell(numel(bands), 1);

for b = 1:numel(bands)
    tag = bands(b).tag;  fs = bands(b).fs;  n_mad = bands(b).n_mad;
    pad_samples = round(pad_ms/1000 * fs);

    % match band variables, skipping anything already cleaned or spike-detected
    is_match = startsWith(all_vars, name_prefix) & contains(all_vars, tag) ...
             & ~contains(all_vars, '_clean');
    matched = all_vars(is_match);

    if isempty(matched)
        warning('No workspace variables found for band "%s". Check names with "who".', tag);
        continue
    end
    fprintf('\n=== %s (fs = %g Hz): %d channel(s) ===\n', tag, fs, numel(matched));

    for k = 1:numel(matched)
        vname = matched{k};
        data  = double(evalin('base', vname));
        data  = data(:);

        [cleaned, bad_idx, ceiling, med] = clean_high_amplitude_mad(data, n_mad, pad_samples);
        N = numel(data);
        clear data

        % Build [voltage, time] in place to avoid holding several full-length copies
        clean_matrix = zeros(N, 2);
        clean_matrix(:, 1) = cleaned;
        clear cleaned
        chunk = 1e7;                                       % fill time column in chunks
        for i0 = 1:chunk:N
            i1 = min(i0 + chunk - 1, N);
            clean_matrix(i0:i1, 2) = t0 + (i0-1:i1-1)' / fs;   % seconds
        end
        assignin('base', [vname '_clean'], clean_matrix);
        clear clean_matrix
        if clear_raw
            evalin('base', ['clear ' vname]);
        end

        pct_flagged = 100 * numel(bad_idx) / N;
        fprintf('%-40s flagged %7d samples (%.3f%%)  ceiling = %.4g\n', ...
            vname, numel(bad_idx), pct_flagged, ceiling);

        results.(vname) = struct('band', tag, 'n_flagged', numel(bad_idx), ...
            'pct_flagged', pct_flagged, 'ceiling', ceiling, 'median', med);

        clear bad_idx
    end
    qc_names{b} = matched{1};
end

%% QC plots: first channel of each band
for b = 1:numel(bands)
    if isempty(qc_names{b}), continue, end
    vname   = qc_names{b};
    cm      = evalin('base', [vname '_clean']);
    has_raw = evalin('base', sprintf('exist(''%s'', ''var'')', vname)) == 1;
    if has_raw, raw_d = evalin('base', vname); end

    % plot at most ~2 million points so long 30 kHz recordings don't exhaust memory
    step = max(1, ceil(size(cm, 1) / 2e6));
    idx  = 1:step:size(cm, 1);

    figure('Name', ['QC: ' vname], 'Color', 'w');
    if has_raw
        plot(cm(idx,2), double(raw_d(idx)), 'Color', [0.7 0.7 0.7]); hold on;
        plot(cm(idx,2), cm(idx,1), 'b');
        legend('raw', 'cleaned');
    else
        plot(cm(idx,2), cm(idx,1), 'b');
        legend('cleaned');
    end
    xlabel('Time (s)'); ylabel('Voltage');
    title(strrep(vname, '_', '\_'), 'Interpreter', 'tex');
    clear raw_d cm idx has_raw
end

%% ---------------- core function ----------------
function [cleaned, bad_idx, ceiling, med] = clean_high_amplitude_mad(data, n_mad, pad_samples)
    % Median + scaled MAD (robust to the artifacts themselves, unlike mean/SD)
    med        = median(data, 'omitnan');
    robust_std = 1.4826 * median(abs(data - med), 'omitnan');   % ~SD for normal data
    ceiling    = n_mad * robust_std;

    bad     = abs(data - med) > ceiling | isnan(data);   % NaNs get filled too, as before
    bad_idx = find(bad);
    cleaned = data;
    if isempty(bad_idx), return, end

    % Runs of flagged samples, each widened by pad_samples on both sides.
    % Works on the (short) list of runs instead of full-length masks/index
    % arrays, so memory stays low even at 30 kHz.
    N  = numel(data);
    rs = find(bad & ~[false; bad(1:end-1)]);   % first sample of each flagged run
    re = find(bad & ~[bad(2:end); false]);     % last sample of each flagged run
    clear bad
    rs = max(rs - pad_samples, 1);
    re = min(re + pad_samples, N);

    % Merge runs that touch or overlap after padding
    new_grp = [true; rs(2:end) > re(1:end-1) + 1];
    last_in = [find(new_grp(2:end)); numel(re)];
    rs = rs(new_grp);
    re = re(last_in);

    % Linear interpolation across each gap, between the good samples either side
    for k = 1:numel(rs)
        s = rs(k);  e = re(k);
        if s == 1 && e == N
            cleaned(:) = NaN;                      % whole channel flagged
        elseif s == 1
            cleaned(s:e) = cleaned(e+1);           % gap at start: hold first good value
        elseif e == N
            cleaned(s:e) = cleaned(s-1);           % gap at end: hold last good value
        else
            a = cleaned(s-1);  b = cleaned(e+1);
            cleaned(s:e) = a + (b - a) * (1:(e-s+1))' / (e-s+2);
        end
    end
end
   save('clean_300_3000hz.mat', '-regexp', '_clean$', '-v7.3')