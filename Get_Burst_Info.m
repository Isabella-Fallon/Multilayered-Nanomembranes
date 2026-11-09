%% Burst detection from spike timestamps, across conditions
% Needs <condition>_<n>_Spikes (timestamps in seconds) in the workspace.
maxISI    = 1;      % s; a gap >= maxISI ends a burst (E_ and E2_ low-frequency events)
maxISI_N  = 1;      % s; same for N_ (neural) spikes. PLACEHOLDER: set from the neural ISI histogram
minSpikes = 2;      % spikes needed to count as a burst (1 = keep isolated spikes)
nChan     = 11;
conditions = {'E_Baseline','E_Poststim','E_PostTTX', ...
   'E2_Poststim','E2_Baseline','E2_PostTTX', ...
    'N_Baseline','N_Poststim','N_PostTTX'};
pattern   = '%s_%d_Spikes';   % e.g. E_Baseline_5_Spikes

Results   = struct();   % Results.<condition>.BurstTable / .BurstEvents
Summary   = table();    % long format: Condition | Channel | NumBursts | MeanDuration
AllEvents = table();    % one row per burst, all conditions

for c = 1:numel(conditions)
    cond = conditions{c};
    thisISI = maxISI;  if startsWith(cond, 'N_'), thisISI = maxISI_N; end
    present = false(1, nChan);
    for n = 1:nChan
        present(n) = exist(sprintf(pattern, cond, n), 'var') == 1;
    end

    startCell = cell(nChan, 1);  endCell = cell(nChan, 1);
    nspCell   = cell(nChan, 1);  durCell = cell(nChan, 1);
    counts    = zeros(nChan, 1);

    for n = 1:nChan
        if ~present(n)
            durCell{n} = NaN;  counts(n) = NaN;   % NaN = variable missing
            continue
        end
        ts = sort(eval(sprintf(pattern, cond, n)));  ts = ts(:);
        if isempty(ts)
            durCell{n} = NaN;  counts(n) = 0;     % variable exists, no spikes
            continue
        end
        gaps   = find(diff(ts) >= thisISI);        % last spike before each gap
        iFirst = [1; gaps + 1];                   % first spike of each candidate burst
        iLast  = [gaps; numel(ts)];               % last spike of each candidate burst
        nsp    = iLast - iFirst + 1;
        keep   = nsp >= minSpikes;
        iFirst = iFirst(keep);  iLast = iLast(keep);  nsp = nsp(keep);
        if isempty(iFirst)
            durCell{n} = NaN;  counts(n) = 0;     % spikes exist, no bursts
            continue
        end
        startCell{n} = ts(iFirst);
        endCell{n}   = ts(iLast);
        nspCell{n}   = nsp;
        durCell{n}   = ts(iLast) - ts(iFirst);    % first spike to last spike
        counts(n)    = numel(iFirst);
    end

    % NaN-padded duration matrix (one row per channel)
    maxN = max([counts; 1], [], 'omitnan');
    Burst_durations = NaN(nChan, maxN);
    for n = 1:nChan
        d = durCell{n}(:).';
        Burst_durations(n, 1:numel(d)) = d;
    end
    Burst_meanDur = mean(Burst_durations, 2, 'omitnan');

    % Wide table: Channel | NumBursts | Dur_1 ... Dur_maxN | MeanDuration
    Channel  = (1:nChan).';
    durNames = compose('Dur_%d', 1:maxN);
    BurstTable = [table(Channel, counts, 'VariableNames', {'Channel','NumBursts'}), ...
        array2table(Burst_durations, 'VariableNames', durNames), ...
        table(Burst_meanDur, 'VariableNames', {'MeanDuration'})];

    % One row per burst
    BurstEvents = table();
    for n = 1:nChan
        if isempty(startCell{n}), continue, end
        k = numel(startCell{n});
        BurstEvents = [BurstEvents; table(repmat(n, k, 1), startCell{n}, endCell{n}, ...
            durCell{n}, nspCell{n}, 'VariableNames', ...
            {'Channel','Start','End','Duration','NSpikes'})]; %#ok<AGROW>
    end

    Results.(cond).BurstTable  = BurstTable;
    Results.(cond).BurstEvents = BurstEvents;
    Summary = [Summary; table(repmat({cond}, nChan, 1), Channel, counts, Burst_meanDur, ...
        'VariableNames', {'Condition','Channel','NumBursts','MeanDuration'})]; %#ok<AGROW>
    if ~isempty(BurstEvents)
        AllEvents = [AllEvents; [table(repmat({cond}, height(BurstEvents), 1), ...
            'VariableNames', {'Condition'}), BurstEvents]]; %#ok<AGROW>
    end
    fprintf('\n=== %s ===\n', cond);
    disp(BurstTable);
end

disp(Summary);
% Optional export
% writetable(Summary,   'BurstSummary.csv');
% writetable(AllEvents, 'BurstEvents_all.csv');