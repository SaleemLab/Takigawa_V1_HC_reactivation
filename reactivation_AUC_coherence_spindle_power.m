%% V1-HC reactivation AUC coherence analysis (Spindle power)

analysis_folder = pwd;

%%%%%%%%%%%%%%% ripple_info
load(fullfile(analysis_folder,'processed_data','ripple_info.mat'));
%%%%%%%%%%%%%%% KDE reactivation bias 
load(fullfile(analysis_folder,'processed_data','track_bias_V1.mat'))
load(fullfile(analysis_folder,'processed_data','track_bias_HC.mat'))

timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Temporal log odds AUC with different spindle power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Spindle Power

fig = figure('Name','Spindle power distribution');
fig = histogram(ripple_info.spindle_amplitude,-2:0.1:5.5,'Normalization','probability');
xlabel('Spindle Power at ripple peak')
ylabel('Proportion of events')
set(gca,'TickDir','out','Box','off','FontSize',12);


spindle_thresholds = prctile(ripple_info.spindle_percentile, 0:99.9/4:99.9);
nBins = length(spindle_thresholds) - 1;

% Time windows
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);
if ~isfile(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_spindle_power.mat'))

    for t = 1:nTime
        t0 = time_bins(t)-win_size/2;
        t1 = time_bins(t) + win_size/2;

        % Sliding V1 window (used for event selection)
        bins_to_select = bin_centers >= t0 & bin_centers < t1;

        fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

        for npower = 1:nBins
            % All events included, spindle bin applied later conditional on track side
            event_index = true(1, length(track_bias_HC));

            % Compute mean log-odds
            mean_bias_V1 = mean(track_bias_V1(bins_to_select, event_index), 'omitnan'); % selector
            mean_bias_HPC = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');       % measure
            total_events = length(mean_bias_V1);
            if total_events < 10, continue; end

            % Quantile thresholds on |V1 bias|
            thresholds = prctile(abs(mean_bias_V1), 0:10:100);
            thresholds = thresholds(1:end-1);
            nThresh = numel(thresholds);

            bias_diff_boot = NaN(nBoot, nThresh);
            bias_diff_shift_boot = NaN(nBoot, nThresh);

            parfor iBoot = 1:nBoot
                s = RandStream('philox4x32_10', 'Seed', iBoot);
                idx = randi(s, total_events, total_events, 1);
                true_idx = find(event_index);

                boot_V1  = mean_bias_V1(idx);
                boot_HPC = mean_bias_HPC(idx);

                % “Shifted”: randomise pairing between V1 & HPC
                boot_V1_shift = mean_bias_V1;
                diff_tmp = NaN(1, nThresh);
                diff_tmp_shift = NaN(1, nThresh);

                for i = 1:nThresh
                    th = thresholds(i);

                    % Identify Track 1 (positive V1 bias) and Track 2 (negative V1 bias)
                    t1 = boot_V1 >= th;     % Track 1
                    t2 = boot_V1 <= -th;    % Track 2

                    % --- Spindle power condition ---
                    % Track 1 → use right probe (2), Track 2 → left probe (1)
                    t1 = t1' + (ripple_info.spindle_percentile(true_idx(idx),2) > spindle_thresholds(npower,2) & ...
                        ripple_info.spindle_percentile(true_idx(idx),2) <= spindle_thresholds(npower+1,2)) > 1';
                    t2 = t2' + (ripple_info.spindle_percentile(true_idx(idx),1) > spindle_thresholds(npower,1) & ...
                        ripple_info.spindle_percentile(true_idx(idx),1) <= spindle_thresholds(npower+1,1)) > 1';

                    % HPC bias difference between Track 1 and Track 2
                    t1_HPC = boot_HPC(t1);
                    t2_HPC = boot_HPC(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_HPC, 'omitnan') - mean(t2_HPC, 'omitnan');
                    end

                    % --- Shifted pairing (null) ---
                    t1s = boot_V1_shift >= th;
                    t2s = boot_V1_shift <= -th;
                    t1s = t1s' + (ripple_info.spindle_percentile(true_idx,2) > spindle_thresholds(npower,2) & ...
                        ripple_info.spindle_percentile(true_idx,2) <= spindle_thresholds(npower+1,2)) > 1';
                    t2s = t2s' + (ripple_info.spindle_percentile(true_idx,1) > spindle_thresholds(npower,1) & ...
                        ripple_info.spindle_percentile(true_idx,1) <= spindle_thresholds(npower+1,1)) > 1';
                    t1_HPCs = boot_HPC(t1s);
                    t2_HPCs = boot_HPC(t2s);
                    if any(t1s) && any(t2s)
                        diff_tmp_shift(i) = mean(t1_HPCs, 'omitnan') - mean(t2_HPCs, 'omitnan');
                    end
                end

                bias_diff_boot(iBoot, :) = diff_tmp;
                bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
            end

            % Quantile-based AUC
            auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
            auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

            % Store summaries
            AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
            AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
            AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
            AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
        end
    end
    save(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_spindle_power.mat'),'AUC')
else
    load(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_spindle_power.mat'))
end

spindle_thresholds = 0:99.9/4:99.9;
nBins = length(spindle_thresholds) - 1;
% Colour scheme
colour_lines = [ ...
    241, 182, 218;
    226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC by spindle power percentile(0.1s win 0.02s step)','Position',[640 100 400 900]);
tiledlayout(nBins,1,'TileSpacing','compact');


for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(sprintf('Spindle power bin %d (%.2f–%.2f)', npower, ...
        spindle_thresholds(npower), spindle_thresholds(npower+1)));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end



% Plot temporal AUC traces
fig = figure('Name','Temporal V1 log-odds AUC low vs high spindle powers','Position',[640 100 1100/3 900/4]);
tiledlayout(nBins, 1, 'TileSpacing','compact');

for npower = [1 4]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    % Real (coloured)
    fill([tvec fliplr(tvec)], ...
        [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    yline(0, '--r');
    xlabel('Time (s relative to ripple onset)');
    ylabel('V1 bias AUC');
    title('Spindle power low vs high');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])
end

% Shifted (black)
fill([tvec fliplr(tvec)], ...
    [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

save_all_figures(fullfile(analysis_folder,'processed_data'),[])


%% Spindle power
% Spindle power binning across both probes
spindle_thresholds = [0:99.9/4:99.9];
nBins = length(spindle_thresholds) - 1;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

spindle_power_KDE_bias_difference = struct;

if ~isfile(fullfile(analysis_folder,'processed_data','spindle_power_KDE_bias_difference_based_on_V1_bias.mat'));

    fig = figure;
    fig.Position = [640 100 1100 650*2];
    fig.Name = 'KDE bias difference in HPC with different spindle powers percentile';
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        tic
        % Select events in current spindle bin
        % power_index = ripple_info.spindle_percentile(:,1) > spindle_thresholds(npower,1) & ...
        %               ripple_info.spindle_percentile(:,1) <= spindle_thresholds(npower+1,1) |...
        %               ripple_info.spindle_percentile(:,2) > spindle_thresholds(npower,2) & ...
        %               ripple_info.spindle_percentile(:,2) <= spindle_thresholds(npower+1,2);

        event_index =1:length(track_bias_HC);
        % event_index = find(power_index);

        mean_bias = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
        % mean_bias_shifted = mean(track_bias_HC(bins_to_use_shifted, event_index), 'omitnan');
        mean_bias_V1 = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
        selected_events = length(mean_bias);

        thresholds = prctile(abs(mean_bias), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = length(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        prop_events_boot = NaN(nBoot, nThresh);
        bias_diff_shifted_boot = NaN(nBoot, nThresh);
        prop_events_shifted_boot = NaN(nBoot, nThresh);


        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, selected_events, selected_events, 1);

            true_idx = find(event_index);

            bb = mean_bias(idx);
            bb_shift = mean_bias;
            boot_V1 = mean_bias_V1(idx);

            diff_tmp = NaN(1, nThresh);
            prop_tmp = NaN(1, nThresh);
            diff_tmp_shifted = NaN(1, nThresh);
            prop_tmp_shifted = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                t1 = bb >= th;
                t2 = bb <= -th;



                t1 = t1' + (ripple_info.spindle_percentile( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                    ripple_info.spindle_percentile( true_idx(idx),2) <= spindle_thresholds(npower+1,2))>1';
                t2 = t2' + (ripple_info.spindle_percentile( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                    ripple_info.spindle_percentile( true_idx(idx),1) <= spindle_thresholds(npower+1,1))>1';

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((ripple_info.spindle_percentile( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                    ripple_info.spindle_percentile( true_idx(idx),2) <= spindle_thresholds(npower+1,2))) ...
                    sum((ripple_info.spindle_percentile( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                    ripple_info.spindle_percentile( true_idx(idx),1) <= spindle_thresholds(npower+1,1)))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s' + (ripple_info.spindle_percentile( true_idx,2) > spindle_thresholds(npower,2) &...
                    ripple_info.spindle_percentile( true_idx,2) <= spindle_thresholds(npower+1,2))>1';
                t2s = t2s' + (ripple_info.spindle_percentile( true_idx,1) > spindle_thresholds(npower,1) &...
                    ripple_info.spindle_percentile( true_idx,1) <= spindle_thresholds(npower+1,1))>1';
                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            prop_events_boot(iBoot, :) = prop_tmp;
            bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
            prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
        end

        % Stats
        bias_mean = mean(bias_diff_boot, 1, 'omitnan');
        bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
        bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
        prop_mean = mean(prop_events_boot, 1, 'omitnan');
        prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
        prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
        bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
        bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
        bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
        prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
        prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
        prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

        % Store
        spindle_power_KDE_bias_difference(npower).power_range = ...
            [spindle_thresholds(npower), spindle_thresholds(npower+1)];
        spindle_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        spindle_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        spindle_power_KDE_bias_difference(npower).prop_mean = prop_mean;
        spindle_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        spindle_power_KDE_bias_difference(npower).thresholds = thresholds;
        spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        spindle_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        spindle_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

        % store AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

        spindle_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        spindle_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        spindle_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        spindle_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

        % ----------- PLOTS (A, B, C) ----------------

        % A: Bias vs. Threshold
        nexttile((npower-1)*3 + 1);
        hold on;
        fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
            colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([thresholds, fliplr(thresholds)], ...
            [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
            [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
        ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
        xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
        title(sprintf('Spindle bin %d: %.2f–%.2f', npower, ...
            spindle_thresholds(npower), spindle_thresholds(npower+1)));
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

        % B: Proportion vs Bias Diff
        nexttile((npower-1)*3 + 2);
        hold on;
        valid = isfinite(bias_mean) & isfinite(prop_mean);
        fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
            [prop_mean(valid), fliplr(prop_mean(valid))], ...
            colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
            [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
            [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Event Proportion vs. Bias Difference');
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

        % C: Proportion CI vs Bias Diff
        nexttile((npower-1)*3 + 3);
        hold on;
        fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
            [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
            colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
            [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
            [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Proportion vs. HPC Bias Difference');
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
        toc
    end
    save(fullfile(analysis_folder,'processed_data','spindle_power_KDE_bias_difference_based_on_V1_bias.mat'),'spindle_power_KDE_bias_difference');
else
    load(fullfile(analysis_folder,'processed_data','spindle_power_KDE_bias_difference_based_on_V1_bias.mat'),'spindle_power_KDE_bias_difference');
end
clear Fill


fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC low vs high spindle power percentile';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = spindle_power_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = spindle_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = spindle_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC spindle power percentile low vs high';
data = spindle_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.3;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    % --- 1. Plot Shuffled Data (Left Bar) ---'
    if i ==1
        x_shuf = i - group_offset;
        y_shuf = data(1).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    end

    % Plot Error Bar
    E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
        'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


    % --- 2. Plot Real Data (Right Bar) ---
    x_real = i + group_offset;
    y_real = data(i).AUC_mean;

    neg_err_real = y_real - data(i).AUC_CI(1);
    pos_err_real = data(i).AUC_CI(2) - y_real;

    % Get specific color for this power bin
    % Assuming colour_lines is size [n_bins x 3]
    this_color = colour_lines(i, :);

    % Plot Bar
    BAR(i+1) = bar(x_real, y_real, bar_width, ...
        'FaceColor', this_color, ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none');

    errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
        'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    % end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','Low','High'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'processed_data'),[])


%% Spindle power PRE
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_percentile, 1);  % avg of probe 1 and 2
spindle_thresholds = prctile(ripple_info.spindle_percentile, 0:99.9/4:99.9);
nBins = length(spindle_thresholds) - 1;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;

spindle_power_KDE_bias_difference = struct;
if ~isfile(fullfile(analysis_folder,'processed_data','spindle_power_KDE_bias_difference_based_on_PRE_V1_bias.mat'));

    fig = figure;
    fig.Position = [640 100 1100 650*2];
    fig.Name = 'KDE bias difference in HPC with different spindle powers percentile (PRE ripple)';
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        tic
        % Select events in current spindle bin
        % power_index = ripple_info.spindle_percentile(:,1) > spindle_thresholds(npower,1) & ...
        %               ripple_info.spindle_percentile(:,1) <= spindle_thresholds(npower+1,1) |...
        %               ripple_info.spindle_percentile(:,2) > spindle_thresholds(npower,2) & ...
        %               ripple_info.spindle_percentile(:,2) <= spindle_thresholds(npower+1,2);

        event_index =1:length(track_bias_HC);
        % event_index = find(power_index);

        mean_bias = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
        % mean_bias_shifted = mean(track_bias_HC(bins_to_use_shifted, event_index), 'omitnan');
        mean_bias_V1 = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
        selected_events = length(mean_bias);

        thresholds = prctile(abs(mean_bias), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = length(thresholds);

        bias_diff_boot = NaN(nBoot, nThresh);
        prop_events_boot = NaN(nBoot, nThresh);
        bias_diff_shifted_boot = NaN(nBoot, nThresh);
        prop_events_shifted_boot = NaN(nBoot, nThresh);


        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, selected_events, selected_events, 1);

            true_idx = find(event_index);

            bb = mean_bias(idx);
            bb_shift = mean_bias;
            boot_V1 = mean_bias_V1(idx);

            diff_tmp = NaN(1, nThresh);
            prop_tmp = NaN(1, nThresh);
            diff_tmp_shifted = NaN(1, nThresh);
            prop_tmp_shifted = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                t1 = bb >= th;
                t2 = bb <= -th;



                t1 = t1' + (ripple_info.spindle_percentile( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                    ripple_info.spindle_percentile( true_idx(idx),2) <= spindle_thresholds(npower+1,2))>1';
                t2 = t2' + (ripple_info.spindle_percentile( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                    ripple_info.spindle_percentile( true_idx(idx),1) <= spindle_thresholds(npower+1,1))>1';

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((ripple_info.spindle_percentile( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                    ripple_info.spindle_percentile( true_idx(idx),2) <= spindle_thresholds(npower+1,2))) ...
                    sum((ripple_info.spindle_percentile( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                    ripple_info.spindle_percentile( true_idx(idx),1) <= spindle_thresholds(npower+1,1)))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s' + (ripple_info.spindle_percentile( true_idx,2) > spindle_thresholds(npower,2) &...
                    ripple_info.spindle_percentile( true_idx,2) <= spindle_thresholds(npower+1,2))>1';
                t2s = t2s' + (ripple_info.spindle_percentile( true_idx,1) > spindle_thresholds(npower,1) &...
                    ripple_info.spindle_percentile( true_idx,1) <= spindle_thresholds(npower+1,1))>1';
                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            prop_events_boot(iBoot, :) = prop_tmp;
            bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
            prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
        end

        % Stats
        bias_mean = mean(bias_diff_boot, 1, 'omitnan');
        bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
        bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
        prop_mean = mean(prop_events_boot, 1, 'omitnan');
        prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
        prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
        bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
        bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
        bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
        prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
        prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
        prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

        % Store
        spindle_power_KDE_bias_difference(npower).power_range = ...
            [spindle_thresholds(npower), spindle_thresholds(npower+1)];
        spindle_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        spindle_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        spindle_power_KDE_bias_difference(npower).prop_mean = prop_mean;
        spindle_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        spindle_power_KDE_bias_difference(npower).thresholds = thresholds;
        spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        spindle_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        spindle_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

        % store AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

        spindle_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        spindle_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        spindle_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        spindle_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

        % ----------- PLOTS (A, B, C) ----------------

        % A: Bias vs. Threshold
        nexttile((npower-1)*3 + 1);
        hold on;
        fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
            colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([thresholds, fliplr(thresholds)], ...
            [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
            [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
        ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
        xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
        title(sprintf('Spindle bin %d: %.2f–%.2f', npower, ...
            spindle_thresholds(npower), spindle_thresholds(npower+1)));
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

        % B: Proportion vs Bias Diff
        nexttile((npower-1)*3 + 2);
        hold on;
        valid = isfinite(bias_mean) & isfinite(prop_mean);
        fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
            [prop_mean(valid), fliplr(prop_mean(valid))], ...
            colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
            [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
            [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Event Proportion vs. Bias Difference');
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

        % C: Proportion CI vs Bias Diff
        nexttile((npower-1)*3 + 3);
        hold on;
        fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
            [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
            colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
            [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
            [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Proportion vs. HPC Bias Difference');
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
        toc
    end

    save(fullfile(analysis_folder,'processed_data','spindle_power_KDE_bias_difference_based_on_PRE_V1_bias.mat'),'spindle_power_KDE_bias_difference');
else
    load(fullfile(analysis_folder,'processed_data','spindle_power_KDE_bias_difference_based_on_PRE_V1_bias.mat'),'spindle_power_KDE_bias_difference');
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC low vs high spindle power percentile (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = spindle_power_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = spindle_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = spindle_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC spindle power percentile low vs high (PRE)';
data = spindle_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.3;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    % --- 1. Plot Shuffled Data (Left Bar) ---'
    if i ==1
        x_shuf = i - group_offset;
        y_shuf = data(1).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    end

    % Plot Error Bar
    E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
        'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


    % --- 2. Plot Real Data (Right Bar) ---
    x_real = i + group_offset;
    y_real = data(i).AUC_mean;

    neg_err_real = y_real - data(i).AUC_CI(1);
    pos_err_real = data(i).AUC_CI(2) - y_real;

    % Get specific color for this power bin
    % Assuming colour_lines is size [n_bins x 3]
    this_color = colour_lines(i, :);

    % Plot Bar
    BAR(i+1) = bar(x_real, y_real, bar_width, ...
        'FaceColor', this_color, ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none');

    errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
        'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    % end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','Low','High'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



save_all_figures(fullfile(analysis_folder,'processed_data'),[])
