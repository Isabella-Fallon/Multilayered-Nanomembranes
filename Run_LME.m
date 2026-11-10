%% Baseline vs Post (stim / TTX): mixed-effects stats + paired plots
% Input: ephys_long.csv (one row per electrode x epoch)
% Columns: Experiment, Measure, Band, Animal, Electrode, Epoch, Value, Count, Duration_s
%
% Model (per comparison):  y ~ Post + (Post|Animal) + (1|Elec)
%   y      = log event rate  (log((Count+0.5)/Duration_s)), or log burst duration
%   Post   = 0 for Baseline, 1 for Post
%   (Post|Animal) lets each animal have its own baseline AND its own effect size
%   (1|Elec)      gives each electrode its own baseline (electrodes nested in animals)
% p-values use Satterthwaite degrees of freedom (appropriate with few animals).
% exp(estimate) = fold change Post/Baseline.

clear; close all;
T = readtable('ephys_long.csv', 'TextType', 'string');
T.Elec   = categorical(string(T.Animal) + "_" + string(T.Electrode));
T.Animal = categorical(T.Animal);
T.Post   = double(T.Epoch == "Post");

sets = unique(T(:, {'Experiment','Measure','Band'}), 'rows');
results = table();

for i = 1:height(sets)
    d = T(T.Experiment == sets.Experiment(i) & T.Measure == sets.Measure(i) & ...
          T.Band == sets.Band(i) & ~isnan(T.Value), :);

    if sets.Measure(i) == "rate"
        d.y = log((d.Count + 0.5) ./ d.Duration_s);   % +0.5 handles zero counts
    else
        d = d(d.Value > 0, :);                         % burst duration must be > 0
        d.y = log(d.Value);
    end

    lme = fitlme(d, 'y ~ Post + (Post|Animal) + (1|Elec)');
    [~, ~, st] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    st = dataset2table(st);
    r  = st(strcmp(st.Name, 'Post'), :);

    results = [results; table(sets.Experiment(i), sets.Measure(i), sets.Band(i), ...
        numel(unique(d.Animal)), numel(unique(d.Elec)), ...
        exp(r.Estimate), exp(r.Lower), exp(r.Upper), r.DF, r.pValue, ...
        'VariableNames', {'Experiment','Measure','Band','nAnimals','nElectrodes', ...
                          'FoldChange','CI95_lo','CI95_hi','DF','p'})]; %#ok<AGROW>
end

disp(results)
writetable(results, 'lme_results.csv');

%% Paired plots: faint lines = electrodes, bold = animal means
panels = {
    'Stim','rate','E2_1-10',     'Stim: 1–10 Hz event rate (Hz)'
    'Stim','rate','E_10-30',     'Stim: 10–30 Hz event rate (Hz)'
    'Stim','rate','N_300-3000',  'Stim: 300–3000 Hz event rate (Hz)'
    'Stim','burstdur','E2_1-10', 'Stim: 1–10 Hz burst duration (s)'
    'Stim','burstdur','E_10-30', 'Stim: 10–30 Hz burst duration (s)'
    'TTX', 'rate','E_10-30',     'TTX: 10–30 Hz event rate (Hz)'
    'TTX', 'rate','N_300-3000',  'TTX: 300–3000 Hz event rate (Hz)'
    'TTX', 'burstdur','E_10-30', 'TTX: 10–30 Hz burst duration (s)'
};
cols = [0.23 0.44 0.71; 0.85 0.51 0.17; 0.23 0.62 0.43];   % m1, m2, m3

figure('Color','w','Position',[100 100 1300 700]);
for k = 1:size(panels,1)
    ax(k) = subplot(2, 4, k); hold on %#ok<SAGROW>
    d = T(T.Experiment == panels{k,1} & T.Measure == panels{k,2} & T.Band == panels{k,3}, :);
    animals = categories(removecats(d.Animal));
    mm = nan(numel(animals), 2);                      % mouse means [Baseline Post]

    for a = 1:numel(animals)
        da = d(d.Animal == animals{a}, :);
        elecs = unique(da.Electrode);
        B = nan(numel(elecs),1); P = B;
        for e = 1:numel(elecs)
            B(e) = da.Value(da.Electrode == elecs(e) & da.Post == 0);
            P(e) = da.Value(da.Electrode == elecs(e) & da.Post == 1);
            plot([1 2], [B(e) P(e)], '-', 'Color', [cols(a,:) 0.25], 'LineWidth', 0.8);
        end
        mm(a,:) = [mean(B,'omitnan') mean(P,'omitnan')];
        h(a) = plot([1 2], mm(a,:), '-o', 'Color', cols(a,:), ...
            'MarkerFaceColor', cols(a,:), 'LineWidth', 2.5, 'MarkerSize', 7); %#ok<SAGROW>
    end
    % bars = mean of mouse means (n = mice), drawn behind the lines
    bh = bar([1 2], mean(mm, 1, 'omitnan'), 0.5, 'FaceColor', [0.88 0.88 0.88], 'EdgeColor', 'none');
    uistack(bh, 'bottom');

    % annotate with LME result
    r = results(results.Experiment == panels{k,1} & results.Measure == panels{k,2} & ...
                results.Band == panels{k,3}, :);
    title({panels{k,4}, sprintf('%.1f× , p = %.3g', r.FoldChange, r.p)}, ...
          'FontWeight','normal','FontSize',9);

    post = "Post-stim"; if panels{k,1} == "TTX", post = "Post-TTX"; end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Baseline', post}, 'XLim', [0.7 2.3], ...
        'Box', 'off', 'TickDir', 'out');
    if panels{k,2} == "burstdur", ylim([0 Inf]); end   % burst-duration axes start at 0
    if k == 1, legend(h, animals, 'Box', 'off', 'Location', 'northwest'); end
end

% Same y axis for the three stim event-rate panels (1-10, 10-30, 300-3000 Hz)
stimRate = find(strcmp(panels(:,1), 'Stim') & strcmp(panels(:,2), 'rate'));
yTop = max(arrayfun(@(a) a.YLim(2), ax(stimRate)));
set(ax(stimRate), 'YLim', [0 yTop]);

%% Panel E: baseline burst duration, 1-10 Hz vs 10-30 Hz (no stim)
% Uses a 300-s baseline (longer than the 180-s epochs above, for more bursts per channel)
% Input: Burst_duration_baseline300s.csv, rows 1-11 m1, 12-22 m2, 23-33 m3
B = readtable('Burst_duration_baseline300s.csv');
B.Properties.VariableNames = {'Band1_10','Band10_30'};
B.Animal    = repelem(["m1";"m2";"m3"], 11);
B.Electrode = repmat((1:11)', 3, 1);
B.Band1_10(B.Band1_10 == 0)   = NaN;     % 0 = no bursts detected
B.Band10_30(B.Band10_30 == 0) = NaN;
B.Animal = string(B.Animal);
B = B(~isnan(B.Band1_10) & ~isnan(B.Band10_30), :);   % channels with bursts in both bands

L = [table(B.Animal, B.Electrode, B.Band1_10,  zeros(height(B),1), 'VariableNames', {'Animal','Electrode','Dur','Hi'});
     table(B.Animal, B.Electrode, B.Band10_30, ones(height(B),1),  'VariableNames', {'Animal','Electrode','Dur','Hi'})];
L.Elec   = categorical(L.Animal + "_" + string(L.Electrode));
L.Animal = categorical(L.Animal);
L.y      = log(L.Dur);

% Same structure as above: each mouse gets its own band difference
lmeE = fitlme(L, 'y ~ Hi + (Hi|Animal) + (1|Elec)');
[~, ~, stE] = fixedEffects(lmeE, 'DFMethod', 'satterthwaite');
stE = dataset2table(stE);
rE  = stE(strcmp(stE.Name, 'Hi'), :);
fprintf('\nPanel E: 10-30 vs 1-10 Hz burst duration: %.2f-fold (95%% CI %.2f-%.2f), p = %.3g, n = %d channels from %d mice\n', ...
    exp(rE.Estimate), exp(rE.Lower), exp(rE.Upper), rE.pValue, height(B), numel(unique(B.Animal)));

figure('Color','w','Position',[100 100 320 360]); hold on
mice = unique(B.Animal);
mmE = nan(numel(mice), 2);
for a = 1:numel(mice)
    ci = str2double(extractAfter(mice(a), "m"));       % m1/m2/m3 -> same colors as other figure
    Ba = B(B.Animal == mice(a), :);
    for e = 1:height(Ba)
        plot([1 2], [Ba.Band1_10(e) Ba.Band10_30(e)], '-', 'Color', [cols(ci,:) 0.25], 'LineWidth', 0.8);
    end
    mmE(a,:) = [mean(Ba.Band1_10) mean(Ba.Band10_30)];
    hE(a) = plot([1 2], mmE(a,:), '-o', 'Color', cols(ci,:), ...
        'MarkerFaceColor', cols(ci,:), 'LineWidth', 2.5, 'MarkerSize', 7); %#ok<SAGROW>
end
bhE = bar([1 2], mean(mmE, 1), 0.5, 'FaceColor', [0.88 0.88 0.88], 'EdgeColor', 'none');
uistack(bhE, 'bottom');
set(gca, 'XTick', [1 2], 'XTickLabel', {'1–10 Hz','10–30 Hz'}, 'XLim', [0.7 2.3], 'Box','off', 'TickDir','out');
ylabel('Burst duration (s)'); ylim([0 Inf]);
title({'Baseline burst duration', sprintf('%.2f× , p = %.3g', exp(rE.Estimate), rE.pValue)}, ...
      'FontWeight','normal','FontSize',9);
legend(hE, mice, 'Box','off', 'Location','northeast');

%% Two-way analysis: band x stimulation (stim event rates)
% Does the stim effect differ between frequency bands?
S = T(T.Experiment == "Stim" & T.Measure == "rate", :);
S.y = log((S.Count + 0.5) ./ S.Duration_s);
S.ChanBand = categorical(string(S.Elec) + "_" + S.Band);   % each channel x band has its own baseline

% ANOVA-style F tests: band, stimulation, band x stimulation
S.Band = categorical(S.Band);
lmeB = fitlme(S, 'y ~ Post*Band + (Post|Animal) + (1|ChanBand)');
fprintf('\n===== Two-way LME (band x stimulation), stim event rates =====\n');
aB = anova(lmeB, 'DFMethod', 'satterthwaite'); disp(aB)
aB = dataset2table(aB);
pInt = aB.pValue(contains(aB.Term, ':'));   % band x stim interaction row

% Pairwise: is the stim effect larger in one band than another?
pairs = {'N_300-3000','E2_1-10'; 'N_300-3000','E_10-30'; 'E_10-30','E2_1-10'};
names = {'300-3000 vs 1-10 Hz'; '300-3000 vs 10-30 Hz'; '10-30 vs 1-10 Hz'};
pw = table();
for j = 1:size(pairs,1)
    % Use only the two bands being compared; Other = 1 for the second band
    Sj = S(ismember(string(S.Band), pairs(j,:)), :);
    Sj.Other = double(string(Sj.Band) == pairs{j,2});
    m  = fitlme(Sj, 'y ~ Post*Other + (Post|Animal) + (1|ChanBand)');
    [~, ~, st] = fixedEffects(m, 'DFMethod', 'satterthwaite');
    st = dataset2table(st);
    r  = st(contains(st.Name, ':'), :);     % the Post x Other interaction term
    % estimate = log FC(second band) - log FC(first band); flip so ratio = first/second
    pw = [pw; table(string(names{j}), exp(-r.Estimate), exp(-r.Upper), exp(-r.Lower), r.DF, r.pValue, ...
        'VariableNames', {'Comparison','RatioOfFoldChanges','CI95_lo','CI95_hi','DF','p'})]; %#ok<AGROW>
end
fprintf('\nPairwise band comparisons of the stim effect (ratio > 1 = first band has the larger effect):\n');
disp(pw)
writetable(pw, 'band_comparisons.csv');

%% Figure: before/after stim for all three bands on one graph
% Bars = mean of mouse means; faint lines = channels; bold lines = mouse means.
bandOrder = {'E2_1-10','E_10-30','N_300-3000'};
bandLab   = {'1–10 Hz','10–30 Hz','300–3000 Hz'};
mice      = categories(removecats(S.Animal));
gap = 3;                                   % x spacing between band groups
figure('Color','w','Position',[100 100 560 380]); hold on
for b = 1:numel(bandOrder)
    x = (b-1)*gap + [1 2];                 % Baseline, Post positions for this band
    mm = nan(numel(mice), 2);
    for a = 1:numel(mice)
        sb = S(string(S.Band) == bandOrder{b} & S.Animal == mice{a}, :);
        ch = unique(sb.Electrode);
        V = nan(numel(ch), 2);
        for e = 1:numel(ch)
            V(e,1) = sb.Value(sb.Electrode == ch(e) & sb.Post == 0);
            V(e,2) = sb.Value(sb.Electrode == ch(e) & sb.Post == 1);
            plot(x, V(e,:), '-', 'Color', [cols(a,:) 0.25], 'LineWidth', 0.8);
        end
        mm(a,:) = mean(V, 1, 'omitnan');
        hG(a) = plot(x, mm(a,:), '-o', 'Color', cols(a,:), 'MarkerFaceColor', cols(a,:), ...
            'LineWidth', 2.5, 'MarkerSize', 6); %#ok<SAGROW>
    end
    bg = bar(x, mean(mm, 1, 'omitnan'), 0.6, 'FaceColor', [0.88 0.88 0.88], 'EdgeColor', 'none');
    uistack(bg, 'bottom');
end
xt = reshape((0:numel(bandOrder)-1)*gap + [1; 2], 1, []);
set(gca, 'XTick', xt, 'XTickLabel', repmat({'Base','Stim'}, 1, numel(bandOrder)), ...
    'XLim', [0.4 (numel(bandOrder)-1)*gap + 2.6], 'Box', 'off', 'TickDir', 'out');
ylim([0 Inf]); ylabel('Event rate (Hz)');
yl = ylim;
for b = 1:numel(bandOrder)                 % band labels above each pair
    text((b-1)*gap + 1.5, yl(2), bandLab{b}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 9);
end
title(sprintf('band × stim interaction p = %.3g', pInt), 'FontWeight', 'normal', 'FontSize', 9);
legend(hG, mice, 'Box', 'off', 'Location', 'northwest');
% For a log axis instead (helps show the smaller 1-10 / 10-30 Hz values):
% set(gca, 'YScale', 'log'); ylim([0.003 Inf]);   % channels at 0 Hz will not be drawn