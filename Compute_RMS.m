% Sliding RMS time series: same rows and duration as the original signal
% Output for each variable: RMS_<name> = [sliding RMS, time (s)]
dt = 1/1500;   % Sampling interval (s)

% Variable name pattern (%d = channel number) and RMS window (s)
signals = { ...
    'raw_%d_bandpass10_30hz_clean', 0.3;   % 10-30 Hz
    'raw_%d_bandpass1_10hz_clean',  1.0};  % 1-10 Hz

for s = 1:size(signals, 1)
    namePattern   = signals{s, 1};
    windowSeconds = signals{s, 2};

    % Odd window length so the window is centered on each sample
    windowSamples = round(windowSeconds / dt);
    if mod(windowSamples, 2) == 0
        windowSamples = windowSamples + 1;
    end

    fprintf('\n--- %s (window = %.3f s, %d samples) ---\n', ...
        strrep(namePattern, '%d', 'N'), windowSeconds, windowSamples);

    for n = 1:11
        varName = sprintf(namePattern, n);
        if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
            voltage = evalin('base', varName);
            voltage = voltage(:);

            % Sliding RMS, one value per input sample
            rmsValues = sqrt(movmean(voltage.^2, windowSamples, ...
                'Endpoints', 'shrink'));

            % Timestamp for every original sample
            timestamps = (0:length(voltage)-1)' * dt;

            % Column 1: sliding RMS, Column 2: time (s)
            rmsOutput = [rmsValues, timestamps];

            newName = ['RMS_' varName];
            assignin('base', newName, rmsOutput);
            fprintf('%s: %d rows\n', newName, length(rmsValues));
        else
            fprintf('Variable not found: %s\n', varName);
        end
    end
end
disp('Sliding RMS calculations complete.');