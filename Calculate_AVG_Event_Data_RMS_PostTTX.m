%% Event and RMS Analysis (PostTTX)
% Step 1: Compute statistics (spike rate, ISI, average waveforms, peak amplitude) per channel
% Step 2: Average sliding RMS over the epoch interval per channel
% Step 3: Visualize results

%% ========== USER CONFIGURATION ==========
% Define which channels to process (1-11)
channels_to_process = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];  % Edit this to select channels

% Variable naming. Names are built as:
%   spikes:    <prefix>_<epoch>_<ch>_Spikes     e.g. N_Poststim_1_Spikes
%   waveforms: <prefix>_<epoch>_<ch>_SpikeWf    e.g. N_Poststim_1_SpikeWf
%   RMS:       <prefix>_<epoch>_<ch>_RMS        e.g. E_Poststim_2_RMS
%              (sliding RMS already interval-filtered in NeuroExplorer)
%              falls back to RMS_raw_<ch>_<band>_clean if not found
%   interval:  <epoch>                          e.g. Poststim
epoch  = 'PostTTX';             % 'Baseline', 'Poststim' or 'PostTTX' (must match your variable names exactly)
prefix = 'E';                   % 'E' (10-30 Hz), 'E2' (1-10 Hz) or 'N' (neural 300-3000 Hz)

% Band that goes with each prefix (used to find the RMS variables)
switch prefix
    case 'E',  band = 'bandpass10_30hz';    sampling_rate = 1500;    % Hz
    case 'E2', band = 'bandpass1_10hz';     sampling_rate = 1500;
    case 'N',  band = 'bandpass300_3000hz'; sampling_rate = 30000;
    otherwise, error('Unknown prefix "%s". Use E, E2 or N.', prefix);
end
fprintf('Prefix %s -> band %s, epoch %s, fs = %d Hz\n', prefix, band, epoch, sampling_rate);

%% ========== INTERVAL DURATION ==========
% Same interval for every channel, so firing rate = spikes / interval length
if ~exist(epoch, 'var')
    error('%s not found in workspace.', epoch);
end
int_var = eval(epoch);
if size(int_var, 2) == 2          % [start end] rows (one row per interval)
    interval_bounds = int_var;
elseif numel(int_var) == 2        % [start; end] column vector
    interval_bounds = int_var(:)';
else
    error('%s has unexpected size [%s]. Expected [start end].', epoch, num2str(size(int_var)));
end
interval_duration = sum(interval_bounds(:,2) - interval_bounds(:,1));
fprintf('Interval (%s) duration: %.2f s\n', epoch, interval_duration);

%% ========== INITIALIZE OUTPUT VARIABLES ==========
% Get waveform length from the first channel that has waveform data
waveform_length = [];
for k = 1:numel(channels_to_process)
    v = sprintf('%s_%s_%d_SpikeWf', prefix, epoch, channels_to_process(k));
    if exist(v, 'var') && ~isempty(eval(v))
        waveform_length = size(eval(v), 1);
        break;
    end
end
if isempty(waveform_length)
    error('Could not find waveform data for any channel. Check variable naming (e.g. %s_%s_1_SpikeWf).', prefix, epoch);
end

% Preallocate output cells
num_channels = length(channels_to_process);
sorted_spike_times = cell(num_channels, 1);
sorted_waveforms = cell(num_channels, 1);

% Statistics output
spike_rates = cell(num_channels, 1);
isi_distributions_per_channel = cell(num_channels, 1);
average_waveforms_per_channel = zeros(waveform_length, num_channels);
std_waveforms_per_channel = zeros(waveform_length, num_channels);
avg_peaks = cell(num_channels, 1);
channel_found = false(num_channels, 1);  % true if channel data exists and passes checks

fprintf('\n========== SPIKE PROCESSING ==========\n');

% ========== PROCESS EACH CHANNEL ==========
for ch_idx = 1:num_channels
    channel_num = channels_to_process(ch_idx);

    % Build variable names
    spike_var = sprintf('%s_%s_%d_Spikes', prefix, epoch, channel_num);
    waveform_var = sprintf('%s_%s_%d_SpikeWf', prefix, epoch, channel_num);

    % Check if variables exist
    if ~exist(spike_var, 'var') || ~exist(waveform_var, 'var')
        warning('Channel %d: variable(s) not found (%s / %s). Skipping...', ...
            channel_num, spike_var, waveform_var);
        continue;
    end

    % Get spike timestamps and waveforms
    spike_times = eval(spike_var);
    waveforms = eval(waveform_var);

    % Skip if either is empty
    if isempty(spike_times) || isempty(waveforms)
        warning('Channel %d: empty spike or waveform data. Skipping...', channel_num);
        continue;
    end

    % Ensure spike_times is a column vector
    spike_times = spike_times(:);

    % Make sure each column = one spike
    if size(waveforms,2) ~= numel(spike_times) && size(waveforms,1) == numel(spike_times)
        waveforms = waveforms';
    end

    % Skip if counts still don't match
    if size(waveforms,2) ~= numel(spike_times)
        warning('Channel %d: %d waveforms but %d spike times. Skipping...', ...
            channel_num, size(waveforms,2), numel(spike_times));
        continue;
    end

    % Skip if waveform length differs from the others
    if size(waveforms,1) ~= waveform_length
        warning('Channel %d: waveform length %d (expected %d). Skipping...', ...
            channel_num, size(waveforms,1), waveform_length);
        continue;
    end

    channel_found(ch_idx) = true;   % only marked found once all checks pass

    % Convert waveforms to microvolts (each column = one spike)
    waveforms_uv = waveforms * 1000;  % Convert mV to µV

    % All spikes are kept (no peak-amplitude exclusion)
    sorted_spike_times{ch_idx} = spike_times;
    sorted_waveforms{ch_idx} = waveforms_uv;

    num_spikes = length(spike_times);

    % SPIKE RATE
    firing_rate = num_spikes / interval_duration;
    spike_rates{ch_idx} = firing_rate;

    % AVERAGE WAVEFORM
    average_waveforms_per_channel(:, ch_idx) = mean(waveforms_uv, 2);

    % STD OF WAVEFORM
    if size(waveforms_uv, 2) > 1
        std_waveforms_per_channel(:, ch_idx) = std(waveforms_uv, 0, 2);
    else
        std_waveforms_per_channel(:, ch_idx) = zeros(waveform_length, 1);
    end

    % AVERAGE AMPLITUDE PEAK
    peak_amplitudes = max(abs(waveforms_uv), [], 1);
    avg_peaks{ch_idx} = mean(peak_amplitudes);

    % ISI DISTRIBUTION
    if num_spikes > 1
        isi_distributions_per_channel{ch_idx} = diff(spike_times);
    else
        isi_distributions_per_channel{ch_idx} = [];
    end

    fprintf('Channel %d: %d spikes, %.2f Hz, Avg Peak: %.2f µV\n', ...
        channel_num, num_spikes, firing_rate, avg_peaks{ch_idx});
end

if ~any(channel_found)
    error('No channels had usable data. Check the warnings above.');
end

fprintf('\n========== STATISTICS ==========\n');
% Display summary statistics per channel
for ch_idx = 1:num_channels
    if ~channel_found(ch_idx), continue; end
    channel_num = channels_to_process(ch_idx);
    firing_rate = spike_rates{ch_idx};
    isis = isi_distributions_per_channel{ch_idx};
    avg_peak = avg_peaks{ch_idx};

    fprintf('\nChannel %d:\n', channel_num);
    fprintf('  Firing rate: %.2f Hz\n', firing_rate);
    fprintf('  Avg Peak Amplitude: %.2f µV\n', avg_peak);

    if ~isempty(isis)
        isis_ms = isis * 1000;
        fprintf('  Mean ISI: %.2f ms\n', mean(isis_ms));
        fprintf('  Median ISI: %.2f ms\n', median(isis_ms));
        fprintf('  Std ISI: %.2f ms\n', std(isis_ms));
        fprintf('  Min ISI: %.2f ms\n', min(isis_ms));
        fprintf('  Max ISI: %.2f ms\n\n', max(isis_ms));
    end
end

% ========== VISUALIZATIONS ==========
figure('Position', [100, 100, 1400, 600]);

% Left subplot: Average waveforms for all channels overlaid
subplot(1, 2, 1);
hold on;
time_axis = (0:waveform_length-1) / sampling_rate * 1000;  % Convert to ms

for ch_idx = 1:num_channels
    if ~channel_found(ch_idx), continue; end
    plot(time_axis, average_waveforms_per_channel(:, ch_idx), 'LineWidth', 2, ...
        'DisplayName', sprintf('Ch %d', channels_to_process(ch_idx)));
end

xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(sprintf('Average Waveforms by Channel (%s, %s)', prefix, epoch));
legend('Location', 'best', 'FontSize', 8);
grid on;
hold off;

% Right subplot: ISI distributions for all channels
subplot(1, 2, 2);
hold on;

for ch_idx = 1:num_channels
    if ~channel_found(ch_idx), continue; end
    isis = isi_distributions_per_channel{ch_idx};

    if ~isempty(isis)
        isis_ms = isis * 1000;
        isis_ms = isis_ms(isis_ms <= 1000);

        bin_width = 1;  % 1 ms bins
        bin_edges = 0:bin_width:1000;
        [counts, ~] = histcounts(isis_ms, bin_edges);

        percentages = (counts / sum(counts)) * 100;
        bin_centers = bin_edges(1:end-1) + bin_width/2;

        plot(bin_centers, percentages, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Ch %d', channels_to_process(ch_idx)));
    end
end

xlabel('ISI (ms)');
ylabel('Percentage (%)');
title(sprintf('ISI Distributions (%s, %s, All Channels)', prefix, epoch));
legend('Location', 'best', 'FontSize', 8);
grid on;
xlim([0, 1000]);
hold off;

%% ========== CREATE SUMMARY STATISTICS TABLE ==========
fprintf('\n========== SUMMARY STATISTICS TABLE ==========\n');
fprintf('Channel | NumSpikes | FiringRate_Hz | AvgPeak_uV | MinFR_Hz | MaxFR_Hz | MeanISI_ms | MedianISI_ms | StdISI_ms | MinISI_ms | MaxISI_ms | AvgRMS\n');
fprintf('--------|-----------|---------------|------------|----------|----------|------------|--------------|-----------|-----------|-----------|-----------\n');

summary_data = [];

for ch_idx = 1:num_channels
    if ~channel_found(ch_idx), continue; end
    channel_num = channels_to_process(ch_idx);

    num_spikes = length(sorted_spike_times{ch_idx});
    firing_rate = spike_rates{ch_idx};
    avg_peak = avg_peaks{ch_idx};
    isis = isi_distributions_per_channel{ch_idx};

    if ~isempty(isis)
        isis_ms = isis * 1000;
        mean_isi   = mean(isis_ms);
        median_isi = median(isis_ms);
        std_isi    = std(isis_ms);
        min_isi    = min(isis_ms);
        max_isi    = max(isis_ms);
        min_fr = 1000 / max_isi;
        max_fr = 1000 / min_isi;
    else
        mean_isi = NaN;
        median_isi = NaN;
        std_isi = NaN;
        min_isi = NaN;
        max_isi = NaN;
        min_fr = NaN;
        max_fr = NaN;
    end

    % AVERAGE RMS over the epoch interval
    % First choice: <prefix>_<epoch>_<ch>_RMS (already filtered to the epoch in NeuroExplorer)
    % Fallback:     RMS_raw_<ch>_<band>_clean (full recording, filtered here by time)
    rms_var     = sprintf('%s_%s_%d_RMS', prefix, epoch, channel_num);
    rms_var_raw = sprintf('RMS_raw_%d_%s_clean', channel_num, band);
    if exist(rms_var, 'var')
        avg_rms = mean_rms(eval(rms_var), [], rms_var, sampling_rate);
    elseif exist(rms_var_raw, 'var')
        avg_rms = mean_rms(eval(rms_var_raw), interval_bounds, rms_var_raw, sampling_rate);
    else
        warning('Neither %s nor %s found. AvgRMS set to NaN for channel %d.', ...
            rms_var, rms_var_raw, channel_num);
        avg_rms = NaN;
    end

    % Print channel results
    fprintf('%7d | %9d | %13.2f | %10.2f | %8.2f | %8.2f | %10.2f | %12.2f | %9.2f | %9.2f | %9.2f | %10.4f\n', ...
        channel_num, num_spikes, firing_rate, avg_peak, ...
        min_fr, max_fr, mean_isi, median_isi, std_isi, min_isi, max_isi, avg_rms);

    summary_data = [summary_data; ...
        channel_num, num_spikes, firing_rate, avg_peak, ...
        min_fr, max_fr, mean_isi, median_isi, std_isi, min_isi, max_isi, avg_rms];
end

% =============================================================
% AVERAGE ACROSS CHANNELS (each channel weighted equally)
% =============================================================
if isempty(summary_data)
    error('No channels had usable data. Check the warnings above.');
end

average_row = [NaN, mean(summary_data(:, 2:end), 1, 'omitnan')];

summary_data = [summary_data; average_row];

fprintf('\n');
fprintf('========== CHANNEL AVERAGE ==========\n');
fprintf('--------|-----------|---------------|------------|----------|----------|------------|--------------|-----------|-----------|-----------|-----------\n');
fprintf('%7s | %9.2f | %13.2f | %10.2f | %8.2f | %8.2f | %10.2f | %12.2f | %9.2f | %9.2f | %9.2f | %10.4f\n', ...
    'AVG', average_row(2:end));

% =============================================================
% CREATE SUMMARY TABLE
% =============================================================
summary_stats = array2table(summary_data, ...
    'VariableNames', {'Channel','NumSpikes', ...
    'FiringRate_Hz','AvgPeak_uV','MinFR_Hz','MaxFR_Hz', ...
    'MeanISI_ms','MedianISI_ms','StdISI_ms', ...
    'MinISI_ms','MaxISI_ms','AvgRMS'});

num_rows = size(summary_data, 1) - 1;  % channels with data (excluding average row)
row_names = cell(num_rows + 1, 1);
for i = 1:num_rows
    row_names{i} = sprintf('Ch%d', summary_data(i,1));
end
row_names{end} = 'Channel_Average';

summary_stats.Properties.RowNames = row_names;

%% ========== SAVE RESULTS TO WORKSPACE ==========
fprintf('\n========== SAVING RESULTS ==========\n');
average_waveforms = average_waveforms_per_channel;
isi_distributions = isi_distributions_per_channel;

% Also save copies named by prefix and epoch, so runs don't overwrite each other
tag = sprintf('%s_%s', prefix, epoch);   % e.g. E2_Baseline
assignin('base', 'average_waveforms', average_waveforms);
assignin('base', 'isi_distributions', isi_distributions);
assignin('base', 'summary_stats', summary_stats);
assignin('base', ['summary_stats_' tag], summary_stats);
assignin('base', ['average_waveforms_' tag], average_waveforms);
assignin('base', ['isi_distributions_' tag], isi_distributions);
assignin('base', 'channels_to_process', channels_to_process);
assignin('base', 'sampling_rate', sampling_rate);
assignin('base', 'waveform_length', waveform_length);

fprintf('Variables saved to workspace:\n');
fprintf('  - average_waveforms (and average_waveforms_%s): [samples x channels]\n', tag);
fprintf('  - isi_distributions (and isi_distributions_%s): [channels x 1] cell array\n', tag);
fprintf('  - summary_stats (and summary_stats_%s): firing rates, peak amplitudes, ISI stats, average RMS\n', tag);
fprintf('  - channels_to_process, sampling_rate, waveform_length\n');

%% ========== HELPER ==========
function avg = mean_rms(x, bounds, var_name, fs)
    % Mean of a sliding RMS variable.
    % x can be:
    %   [RMS, time (s)] or [time (s), RMS]  - time column = the one that increases steadily
    %   a plain RMS vector                  - time = (0:N-1)/fs, only used if bounds given
    % bounds = []           -> x is already restricted to the epoch; average everything
    % bounds = [start end]  -> average only samples whose time is inside the interval(s)
    if ~isnumeric(x) || isempty(x)
        warning('%s is empty or not numeric. AvgRMS = NaN.', var_name);
        avg = NaN;
        return;
    end

    if isvector(x)
        r = x(:);
        t = (0:numel(r)-1)' / fs;
    elseif size(x, 2) == 2
        col1_is_time = all(diff(x(:,1)) > 0);
        col2_is_time = all(diff(x(:,2)) > 0);
        if col2_is_time && ~col1_is_time
            r = x(:,1);  t = x(:,2);      % [RMS, time]
        elseif col1_is_time && ~col2_is_time
            t = x(:,1);  r = x(:,2);      % [time, RMS]
        else
            error(['%s: could not tell which column is time ' ...
                   '(neither or both columns increase steadily).'], var_name);
        end
    else
        error('%s should be N x 2 [RMS, time] or a vector; got size [%s].', ...
            var_name, num2str(size(x)));
    end

    if isempty(bounds)
        avg = mean(r, 'omitnan');
        return;
    end

    in_interval = false(size(t));
    for k = 1:size(bounds, 1)
        in_interval = in_interval | (t >= bounds(k,1) & t <= bounds(k,2));
    end
    if ~any(in_interval)
        warning(['%s has no samples inside the interval (time runs %.2f-%.2f s; ' ...
                 'interval starts at %.2f s). AvgRMS = NaN.'], ...
                 var_name, t(1), t(end), bounds(1,1));
        avg = NaN;
    else
        avg = mean(r(in_interval), 'omitnan');
    end
end