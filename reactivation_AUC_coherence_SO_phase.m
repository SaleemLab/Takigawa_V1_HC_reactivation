%% V1-HC reactivation AUC coherence analysis (SO phase)

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

%% Temporal log odds AUC with SO peak vs SO trough
timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;

SO_thresholds = {'peak','trough'};
nBins = 2;
nBoot = 1000;

win_size = 0.1;
step_size = 0.02;

time_bins = -1:step_size:1;
nTime = numel(time_bins);
bins_to_use = bin_centers>0 & bin_centers<0.1;

AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
    ] / 256;

if ~isfile(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_SO_phase.mat'))
    for t = 1:nTime
        t0 = time_bins(t);
        t1 = t0 + win_size;
        bins_to_select = bin_centers >= t0 & bin_centers < t1;

        fprintf('SO phase: processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);
        tic
        for npower = 1:nBins
            event_index =1:length(track_bias_HC);
            % event_index = find(power_index);


            mean_bias_V1   = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
            % mean_bias_shifted = mean(track_bias_HC(bins_to_use_shifted, event_index), 'omitnan');
            mean_bias_HPC = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
            mean_bias = mean_bias_V1;

            selected_events = length(mean_bias);

            thresholds = prctile(abs(mean_bias), 0:10:100);
            thresholds = thresholds(1:end-1);
            nThresh = length(thresholds);

            bias_diff_boot = NaN(nBoot, nThresh);
            bias_diff_shifted_boot = NaN(nBoot, nThresh);
            % event_phase = ripple_info.SO_phase';

            parfor iBoot = 1:nBoot
                s = RandStream('philox4x32_10', 'Seed', iBoot);
                idx = randi(s, selected_events, selected_events, 1);

                true_idx = find(event_index);
                boot_HPC = mean_bias_HPC(idx);

                bb_shift = mean_bias_V1;
                boot_V1 = mean_bias_V1(idx);

                diff_tmp = NaN(1, nThresh);
                % prop_tmp = NaN(1, nThresh);
                diff_tmp_shifted = NaN(1, nThresh);
                % prop_tmp_shifted = NaN(1, nThresh);
                event_phase = ripple_info.SO_phase(true_idx(idx),:)';
                event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

                for i = 1:nThresh
                    th = thresholds(i);

                    t1 = boot_V1 >= th;
                    t2 = boot_V1 <= -th;

                    if npower == 1 % if phase peak


                        t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                        t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                        t1_V1 = boot_HPC(t1);
                        t2_V1 = boot_HPC(t2);
                        if any(t1) && any(t2)
                            diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                        end

                        % total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                        %     sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                        % prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;


                        t1s = bb_shift >= th;
                        t2s = bb_shift <= -th;
                        t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                        t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                        t1_V1 = boot_HPC(t1s);
                        t2_V1 = boot_HPC(t2s);

                        if any(t1s) && any(t2s)
                            diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                        end
                        % prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                    elseif npower == 2 % if phase trough

                        t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                        t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                        t1_V1 = boot_HPC(t1);
                        t2_V1 = boot_HPC(t2);
                        if any(t1) && any(t2)
                            diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                        end

                        total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                            sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                        % prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

                        t1s = bb_shift >= th;
                        t2s = bb_shift <= -th;
                        t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                        t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                        t1_V1 = boot_HPC(t1s);
                        t2_V1 = boot_HPC(t2s);

                        if any(t1s) && any(t2s)
                            diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                        end
                        % prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                    end
                end
                bias_diff_boot(iBoot, :) = diff_tmp;
                bias_diff_shift_boot(iBoot, :) = diff_tmp_shifted;
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
    save(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_SO_phase.mat'),'AUC')
else
    % save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_SO_phase.mat'),'AUC')
    load(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_SO_phase.mat'))
end

% -------- Plot --------
fig = figure('Name','Temporal HPC log-odds AUC SO peak vs SO trough (0.1s win 0.02s step)','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');

% figure
for npower = 1:nBins
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),2);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),2,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    fill([tvec fliplr(tvec)], [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);


    yline(0,'--r');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title('SO peak vs SO trough');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.25])

    xline(0,'--k');
end

fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);


% Save results
save_all_figures(fullfile(analysis_folder,'processed_data'),[])


%% SO phases
%%%%%%%%%%%%%%%%
% SO phase binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak','trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
    ] / 256;

if ~isfile(fullfile(analysis_folder,'processed_data','SO_phase_KDE_bias_difference_based_on_V1_bias.mat'))
    fig = figure;
    fig.Position = [640 100 1100 650];
    fig.Name = 'KDE bias difference in HPC with different SO phase';
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        tic

        event_index =1:length(track_bias_HC);
        % event_index = find(power_index);

        mean_bias = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
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
            event_phase = ripple_info.SO_phase(true_idx(idx),:)';
            event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

            for i = 1:nThresh
                th = thresholds(i);

                t1 = bb >= th;
                t2 = bb <= -th;

                if npower == 1 % if phase peak


                    t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                    t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                    t1_V1 = boot_V1(t1);
                    t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                    end

                    total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                        sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;


                    t1s = bb_shift >= th;
                    t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                    t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                    t1_V1 = boot_V1(t1s);
                    t2_V1 = boot_V1(t2s);

                    if any(t1s) && any(t2s)
                        diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                    end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                elseif npower == 2 % if phase trough

                    t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                    t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                    t1_V1 = boot_V1(t1);
                    t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                    end

                    total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                        sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                    t1s = bb_shift >= th;
                    t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                    t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                    t1_V1 = boot_V1(t1s);
                    t2_V1 = boot_V1(t2s);

                    if any(t1s) && any(t2s)
                        diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                    end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                end



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
        SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
            SO_thresholds{npower};
        SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
        SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

        % store AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

        SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);


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
        title(sprintf('SO phase bin %s', SO_thresholds{npower}));
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
    save(fullfile(analysis_folder,'processed_data','SO_phase_KDE_bias_difference_based_on_V1_bias.mat'),'SO_phase_KDE_bias_difference');
else
    load(fullfile(analysis_folder,'processed_data','SO_phase_KDE_bias_difference_based_on_V1_bias.mat'),'SO_phase_KDE_bias_difference');
end
clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC SO phase';

% Select bins (1 = low, 4 = high)
subplot(1,2,1);
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

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

bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');



%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO peak vs trough';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.3;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:2
    % --- 1. Plot Shuffled Data (Left Bar) ---
    if i == 1
        x_shuf = i - group_offset;
        y_shuf = data(end).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(end).AUC_CI_shuffled(2) - y_shuf;


        BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end




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

end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO peak','SO trough'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



save_all_figures(fullfile(analysis_folder,'processed_data'),[])



%%%%%%%%%%%%%%%%%%%%%%%% SO phase binning PRE
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak','trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_phase_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
    ] / 256;

if ~isfile(fullfile(analysis_folder,'processed_data','SO_phase_KDE_bias_difference_based_on_PRE_V1_bias.mat'))

    fig = figure;
    fig.Position = [640 100 1100 650];
    fig.Name = 'KDE bias difference in HPC with different SO phase (PRE ripple)';
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        tic
        % Select events in current spindle bin
        % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
        %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
        %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
        %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

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


        % event_phase = ripple_info.SO_phase';


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
            event_phase = ripple_info.SO_phase(true_idx(idx),:)';
            event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

            for i = 1:nThresh
                th = thresholds(i);

                t1 = bb >= th;
                t2 = bb <= -th;

                if npower == 1 % if phase peak


                    t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                    t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                    t1_V1 = boot_V1(t1);
                    t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                    end

                    total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                        sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                    t1s = bb_shift >= th;
                    t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                    t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                    t1_V1 = boot_V1(t1s);
                    t2_V1 = boot_V1(t2s);

                    if any(t1s) && any(t2s)
                        diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                    end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                elseif npower == 2 % if phase trough

                    t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                    t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                    t1_V1 = boot_V1(t1);
                    t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2)
                        diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                    end

                    total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                        sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                    t1s = bb_shift >= th;
                    t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                    t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                    t1_V1 = boot_V1(t1s);
                    t2_V1 = boot_V1(t2s);

                    if any(t1s) && any(t2s)
                        diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                    end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                end



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
        SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
            SO_thresholds{npower};
        SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
        SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];


        % store AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

        SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);


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
        title(sprintf('SO phase bin %s', SO_thresholds{npower}));
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
    save(fullfile(analysis_folder,'processed_data','SO_phase_KDE_bias_difference_based_on_PRE_V1_bias.mat'),'SO_phase_KDE_bias_difference');
else
    load(fullfile(analysis_folder,'processed_data','SO_phase_KDE_bias_difference_based_on_PRE_V1_bias.mat'),'SO_phase_KDE_bias_difference');
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO phase (PRE ripple)';

% Select bins (1 = low, 4 = high)
subplot(1,2,1);
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

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

bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');



%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO peak vs trough (PRE)';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.3;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:2
    if i == 1
        x_shuf = i - group_offset;
        y_shuf = data(end).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(end).AUC_CI_shuffled(2) - y_shuf;


        BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
        
        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end



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
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO peak','SO trough'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'processed_data'),[])



%% Plot V1 MUA SO phase plot

%%%%% Plot MUA SO phase relationship
load(fullfile(analysis_folder,'processed_data','SO_phase_MUA_spike_rate.mat'),'SO_phase_MUA_spike_rate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_phase_bins = 10;
phase_bin_edges = -pi:pi/no_phase_bins:pi;
phase_bins = -pi+pi/no_phase_bins/2:pi/no_phase_bins:pi-pi/no_phase_bins/2;


V1_phase = [];

for hemi = 1:2
    
    for nsession = 1:22
        V1_phase{hemi}(1,nsession,:) = SO_phase_MUA_spike_rate{nsession}(hemi,1,:);
        V1_phase{hemi}(2,nsession,:) = SO_phase_MUA_spike_rate{nsession}(hemi,2,:);
        %         V1_phase{hemi}(1,nession,:) = normalize(SO_phase_MUA_spike_rate{nsession}(hemi,1,:),'range');
        % V1_phase{hemi}(2,nession,:) = normalize(SO_phase_MUA_spike_rate{nsession}(hemi,2,:),'range');
    end
end

mean_phases_ipsi = mean([squeeze((V1_phase{1}(1,:,:))); squeeze((V1_phase{2}(2,:,:)))]);
mean_phases_contra = mean([squeeze((V1_phase{1}(2,:,:))); squeeze((V1_phase{2}(1,:,:)))]);

SE_phases_ipsi = std([squeeze((V1_phase{1}(1,:,:))); squeeze((V1_phase{2}(2,:,:)))])/sqrt(length(SO_phase_MUA_spike_rate{nsession}));
SE_phases_contra = std([squeeze((V1_phase{1}(2,:,:))); squeeze((V1_phase{2}(1,:,:)))])/sqrt(length(SO_phase_MUA_spike_rate{nsession}));
% 
%
%%% Plot SO V1 phase

fig = figure;
fig.Name = 'V1 SO phase Peak vs Trough';
% fig.Name = 'V1 SO phase Peak vs Trough (DOWN 2Hz)';
sgtitle(fig.Name)
fig.Position= [844 66 560 906];


% 1. Setup Colors (Example: Teal for Ipsi, Orange for Contra)
c_ipsi = [35,139,69]/256;
c_contra = [106,81,163]/256;

% 1. Align Data and Phase Bins
% Find the index for -pi/2 to center the trough at 0
[~, shift_idx] = min(abs(phase_bins - (-pi/2)));

% Shift the Firing Rate data
m_i_s = circshift(mean_phases_ipsi, [0, -(shift_idx-1)]);
m_c_s = circshift(mean_phases_contra, [0, -(shift_idx-1)]);
s_i_s = circshift(SE_phases_ipsi, [0, -(shift_idx-1)]);
s_c_s = circshift(SE_phases_contra, [0, -(shift_idx-1)]);

% 2. Construct the Continuous X-Axis
% Calculate actual bin width to prevent "drift" over multiple cycles
dx = mean(diff(phase_bins));

% Create a single cycle
x_single = ((0:length(m_i_s)-1) * dx);

% Concatenate for 2 full cycles (ensuring the second cycle starts exactly one bin after the first ends)
x_2cyc = [x_single, x_single + (x_single(end) + dx)];

% Concatenate the Data
mean_ipsi_2cyc = [m_i_s,m_i_s];
mean_contra_2cyc = [m_c_s,m_c_s];
SE_ipsi_2cyc = [s_i_s,s_i_s];
SE_contra_2cyc = [s_c_s,s_c_s];

% 3. Plotting with Precise Bounds
ipsi_upper = mean_ipsi_2cyc + SE_ipsi_2cyc;
ipsi_lower = mean_ipsi_2cyc - SE_ipsi_2cyc;
contra_upper = mean_contra_2cyc + SE_contra_2cyc;
contra_lower = mean_contra_2cyc - SE_contra_2cyc;

subplot(3,2,1)
hold on;
% Shading
fill([x_2cyc, fliplr(x_2cyc)], [ipsi_lower, fliplr(ipsi_upper)], ...
    c_ipsi, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
fill([x_2cyc, fliplr(x_2cyc)], [contra_lower, fliplr(contra_upper)], ...
    c_contra, 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Mean Lines
h1 = plot(x_2cyc, mean_ipsi_2cyc, 'Color', c_ipsi, 'LineWidth', 2.5);
h2 = plot(x_2cyc, mean_contra_2cyc, 'Color', c_contra, 'LineWidth', 2.5);

% Vertical Grid lines at intervals of pi
for v = 0:pi:3*pi
    line([v v], [1 6], 'Color', [0.7 0.7 0.7], 'LineStyle', '--', 'HandleVisibility', 'off');
end

% Formatting
xticks(pi/2:pi/2:5*pi/2);
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2'});
xlim([pi/2 5*pi/2]);
ylim([1 6]);
ylabel('Firing Rate (Hz)');
xlabel('Phase (rad)');
set(gca, 'TickDir', 'out', 'Box', 'off');
hold off;


% SO_phase_DOWN_MUA_spike_rate

%%%%%% paired t test
for nsession = 1:22
    for hemi = 1:2
        temp = squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,1,:));
        V1_SO_FR(nsession,hemi,1,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
        temp = squeeze(SO_phase_MUA_spike_rate{nsession}(hemi,2,:));
        V1_SO_FR(nsession,hemi,2,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];

        % temp = squeeze(SO_phase_DOWN_HPC_MUA_spike_rate{nsession}(hemi,1,:));
        % HPC_SO_FR(nsession,hemi,1,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
        % temp = squeeze(SO_phase_DOWN_HPC_MUA_spike_rate{nsession}(hemi,2,:));
        % HPC_SO_FR(nsession,hemi,2,:) = [mean(temp(phase_bins>-pi/2 & phase_bins<pi/2)) mean(temp(phase_bins<-pi/2 | phase_bins>pi/2))];
    end
end

peak_vs_trough_ipsi = [];
peak_vs_trough_contra = [];
peak_vs_trough_ipsi(:,1) =  mean([squeeze(V1_SO_FR(:,1,1,1)) squeeze(V1_SO_FR(:,2,2,1))],2);
peak_vs_trough_ipsi(:,2) =  mean([squeeze(V1_SO_FR(:,1,1,2)) squeeze(V1_SO_FR(:,2,2,2))],2);

peak_vs_trough_contra(:,1) =  mean([squeeze(V1_SO_FR(:,1,2,1)) squeeze(V1_SO_FR(:,2,1,1))],2);
peak_vs_trough_contra(:,2) =  mean([squeeze(V1_SO_FR(:,1,2,2)) squeeze(V1_SO_FR(:,2,1,2))],2);


p=[];
p(1) =signrank(peak_vs_trough_ipsi(:,1),peak_vs_trough_ipsi(:,2));
p(2) =signrank(peak_vs_trough_contra(:,1),peak_vs_trough_contra(:,2));

p(3) =signrank(peak_vs_trough_ipsi(:,1),peak_vs_trough_contra(:,1));% ipsi vs contra peak
p(4) =signrank(peak_vs_trough_ipsi(:,2),peak_vs_trough_contra(:,2));% ipsi vs contra trough

% peak_vs_trough =  squeeze(V1_SO_FR(:,2,1,:));
% peak_vs_trough =  squeeze(V1_SO_FR(:,2,2,:));


%%% Plot peak vs trough phases
% 1. Setup the figure and colors
peak_vs_trough = {peak_vs_trough_ipsi,peak_vs_trough_contra,[peak_vs_trough_ipsi(:,1) peak_vs_trough_contra(:,1)],[peak_vs_trough_ipsi(:,2) peak_vs_trough_contra(:,2)]};
conditions = {'V1 ipsi','V1 contra','V1 ipsi vs contra peak','V1 ipsi vs contra trough'};

% fig = figure;
% fig.Name = 'V1 SO phase Peak vs Trough';
% sgtitle(fig.Name)
% fig.Position= [844 66 560 906];
% counter = 1;

for iplot = 1:4
    subplot(3,2,iplot+2)
    hold on;

    % Define colors
    color_col1 = [0.1 0.5 0.9];  % Blue for condition 1
    color_col2 = [0.9 0.4 0.1];  % Orange for condition 2
    color_link = [0.5 0.5 0.5];  % Gray linking lines

    % Get the number of pairs (rows)
    num_pairs = size(peak_vs_trough{iplot}, 1);

    rng(1);
    jitter_col1 = (rand(num_pairs, 1) - 0.5) * 0.2;
    rng(2);
    jitter_col2 = (rand(num_pairs, 1) - 0.5) * 0.2;

    % 2. Draw the linking lines (One line per pair)
    for i = 1:num_pairs
        % X-values are fixed at [1 2] (the categories)
        % Y-values are the data [value_at_1 value_at_2]
        plot([1 + jitter_col1(i), 2 + jitter_col2(i)], peak_vs_trough{iplot}(i, :), ...
            'Color', color_link, ...
            'LineWidth', 0.5);
    end

    % 3. Plot the scatter points on top
    % Scatter for Condition 1 (at X=1)
    scatter(ones(num_pairs, 1)+jitter_col1, peak_vs_trough{iplot}(:, 1), ...
        70, color_col1, 'filled', ...
        'DisplayName', 'Condition 1');

    % Scatter for Condition 2 (at X=2)
    scatter(2 * ones(num_pairs, 1)+jitter_col2, peak_vs_trough{iplot}(:, 2), ...
        70, color_col2, 'filled', ...
        'DisplayName', 'Condition 2');


    % 4. Formatting and Labels
    ylabel('Mean firing rate (Hz)');
    % xlabel('Condition');
    title(conditions{iplot});

    % Set X-axis to show categories 1 and 2 clearly
    xlim([0.5 2.5]);
    xticks([1 2]);
    if iplot <3
        xticklabels({'Peak', 'Trough'}); % Use descriptive labels
    else
        xticklabels({'Ipsi', 'contra'}); % Use descriptive labels
    end

    text(1.5,4,sprintf('p = %.3e',p(iplot)))
    % Clean up the axes
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    grid on;

    legend('off'); % Turn off legend since the categories are labeled on the X-axis
    hold off;
end

save_all_figures(fullfile(analysis_folder,'processed_data'),[])

