%% V1-HC reactivation AUC coherence analysis (SO Trough Synchrony)

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


%% Distribution of SO phases (matched/dominant vs non-matched/non-dominant)
bins_to_use = bin_centers > 0 & bin_centers < 0.1;
V1_bias = mean(track_bias_V1(bin_centers > -0.2 & bin_centers < 0.0, :), 'omitnan');

% Indices for Matched (Dominant) and Non-matched (non-dominant)
matched_trough = {find(V1_bias < 0 & (ripple_info.SO_phase(:,1)' < -pi/2 | ripple_info.SO_phase(:,1)' > pi/2)), ...
                  find(V1_bias > 0 & (ripple_info.SO_phase(:,2)' < -pi/2 | ripple_info.SO_phase(:,2)' > pi/2))};

matched_peak = {find(V1_bias < 0 & (ripple_info.SO_phase(:,1)' < pi/2 & ripple_info.SO_phase(:,1)' > -pi/2)), ...
                find(V1_bias > 0 & (ripple_info.SO_phase(:,2)' < pi/2 & ripple_info.SO_phase(:,2)' > -pi/2))};

matched = {[matched_trough{1} matched_peak{1}], [matched_trough{2} matched_peak{2}]};

% Shifting by PI to center the Trough
% Map original -pi:pi to 0:2pi. Original +/-pi becomes 0/2pi, original 0 becomes pi.
% To put the TROUGH (pi) in the center, we shift the raw data so that 
% the region of interest is centered on the plot limits.
shifted_phase_matched = [mod(ripple_info.SO_phase(matched{1},1), 2*pi); ...
                         mod(ripple_info.SO_phase(matched{2},2), 2*pi)];
                     
shifted_phase_non_matched = [mod(ripple_info.SO_phase(matched{1},2), 2*pi); ...
                             mod(ripple_info.SO_phase(matched{2},1), 2*pi)];

% Binning and Histograms
bin_width = pi/6;
edges = 0:bin_width:2*pi;
[N, Xedges, Yedges] = histcounts2(shifted_phase_matched, shifted_phase_non_matched, edges, edges);

% Center the bin labels for imagesc
centers = edges(1:end-1) + bin_width/2;

% Plotting
figure('Name','Distribution of Bilateral V1 SO phases dominant vs non-dominant (Trough Centered)')
imagesc(centers, centers, N/sum(N(:))); 
colorbar; 
colormap(flipud(gray));
clim([0 0.033]);

% Set ticks to show the shift (Centering PI)
xticks(0:pi/2:2*pi);
xticklabels({'0', 'pi/2', 'pi (Trough)', '3pi/2', '2pi'});
yticks(0:pi/2:2*pi);
yticklabels({'0', 'pi/2', 'pi (Trough)', '3pi/2', '2pi'});

xlabel('Dominant V1 SO phase (ripple)');
ylabel('Non-dominant V1 SO phase (ripple)');

% Draw lines to highlight the central trough quadrant (pi/2 to 3pi/2)
xline(pi/2, 'r-', 'LineWidth', 2); xline(3*pi/2, 'r-', 'LineWidth', 2);
yline(pi/2, 'r-', 'LineWidth', 2); yline(3*pi/2, 'r-', 'LineWidth', 2);

% Formatting
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
set(gca, 'YDir', 'normal'); % Standard Cartesian orientation

%% Temporal log odds AUC with SO trough (Unilateral vs Bilateral)
SO_thresholds = {'trough unique','trough both'};
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

if ~isfile(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_SO_trough_synchrony.mat'))
    for t = 1:nTime
        t0 = time_bins(t);
        t1 = t0 + win_size;
        bins_to_select = bin_centers >= t0 & bin_centers < t1;

        fprintf('SO Trough Synchrony: processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);
        tic
        for npower = 1:nBins
            event_index = 1:length(track_bias_HC);

            mean_bias_V1 = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
            mean_bias_HPC = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
            mean_bias = mean_bias_V1;

            selected_events = length(mean_bias);

            thresholds = prctile(abs(mean_bias), 0:10:100);
            thresholds = thresholds(1:end-1);
            nThresh = length(thresholds);

            bias_diff_boot = NaN(nBoot, nThresh);
            bias_diff_shift_boot = NaN(nBoot, nThresh);

            parfor iBoot = 1:nBoot
                s = RandStream('philox4x32_10', 'Seed', iBoot);
                idx = randi(s, selected_events, selected_events, 1);

                true_idx = find(event_index);
                boot_HPC = mean_bias_HPC(idx);

                bb_shift = mean_bias;
                boot_V1 = mean_bias_V1(idx);

                diff_tmp = NaN(1, nThresh);
                diff_tmp_shifted = NaN(1, nThresh);
                event_phase = ripple_info.SO_phase(true_idx(idx),:)';
                event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

                for i = 1:nThresh
                    th = thresholds(i);

                    t1 = boot_V1 >= th;
                    t2 = boot_V1 <= -th;

                    if npower == 1 % if trough not sync (Unilateral)
                        t1 = t1 & (event_phase(1,:) > -pi/2 & event_phase(1,:) < pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);
                        t2 = t2 & (event_phase(2,:) > -pi/2 & event_phase(2,:) < pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);

                        t1_V1 = boot_HPC(t1);
                        t2_V1 = boot_HPC(t2);
                        if any(t1) && any(t2)
                            diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                        end

                        t1s = bb_shift >= th;
                        t2s = bb_shift <= -th;
                        t1s = t1s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) > pi/2 & event_phase_shifted(2,:) <= pi);
                        t2s = t2s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) < -pi/2 | event_phase_shifted(1,:) > pi/2 & event_phase_shifted(1,:) <= pi);

                        t1_V1 = boot_HPC(t1s);
                        t2_V1 = boot_HPC(t2s);

                        if any(t1s) && any(t2s)
                            diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                        end

                    elseif npower == 2 % if trough sync (Bilateral)
                        t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                        t2 = t2 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                        t1_V1 = boot_HPC(t1);
                        t2_V1 = boot_HPC(t2);
                        if any(t1) && any(t2)
                            diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                        end

                        t1s = bb_shift >= th;
                        t2s = bb_shift <= -th;
                        t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                        t2s = t2s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                        t1_V1 = boot_HPC(t1s);
                        t2_V1 = boot_HPC(t2s);

                        if any(t1s) && any(t2s)
                            diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                        end
                    end
                end
                bias_diff_boot(iBoot, :) = diff_tmp;
                bias_diff_shift_boot(iBoot, :) = diff_tmp_shifted;
            end

            % Quantile-based AUC calculation
            auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
            auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

            % Store Moving Estimates
            AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
            AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
            AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
            AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
        end
        toc
    end
    save(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_SO_trough_synchrony.mat'),'AUC')
else
    load(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_SO_trough_synchrony.mat'))
end

% -------- Moving Windows Visualization --------
fig = figure('Name','Temporal HPC log-odds AUC SO Trough Synchrony','Position',[640 100 400 900/4]);
tiledlayout(nBins,1,'TileSpacing','compact');
clear p
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
    p(npower)=plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    fill([tvec fliplr(tvec)], [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0,'--r');
    xline(0,'--k');
    xlabel('Time (s relative to ripple)');
    ylabel('HPC bias AUC');
    title(SO_thresholds{npower});
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.26]);
end
legend([p(1:2)],SO_thresholds{1},SO_thresholds{2},'box','off')
save_all_figures(fullfile(analysis_folder,'processed_data'),[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SO Trough Synchrony Binning (POST-Ripple Window Analysis)

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;

SO_phase_KDE_bias_difference = struct;

if ~isfile(fullfile(analysis_folder,'processed_data','SO_trough_KDE_bias_difference_based_on_V1_bias.mat'))
    fig = figure('Name','KDE bias difference in HPC with different SO trough synchrony','Position',[640 100 1100 650]);
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        tic
        event_index = 1:length(track_bias_HC);

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

                if npower == 1 % trough unilateral
                    t1 = t1 & (event_phase(1,:) > -pi/2 & event_phase(1,:) < pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);
                    t2 = t2 & (event_phase(2,:) > -pi/2 & event_phase(2,:) < pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1); t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2), diff_tmp(i) = mean(t1_V1) - mean(t2_V1); end

                    total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                        sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);
                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

                    t1s = bb_shift >= th; t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) > pi/2 & event_phase_shifted(2,:) <= pi);
                    t2s = t2s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) < -pi/2 | event_phase_shifted(1,:) > pi/2 & event_phase_shifted(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1s); t2_V1 = boot_V1(t2s);
                    if any(t1s) && any(t2s), diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1); end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                elseif npower == 2 % trough bilateral
                    t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                    t2 = t2 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1); t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2), diff_tmp(i) = mean(t1_V1) - mean(t2_V1); end

                    total_events = sum((event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi));
                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

                    t1s = bb_shift >= th; t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                    t2s = t2s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1s); t2_V1 = boot_V1(t2s);
                    if any(t1s) && any(t2s), diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1); end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
                end
            end
            bias_diff_boot(iBoot, :) = diff_tmp;
            prop_events_boot(iBoot, :) = prop_tmp;
            bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
            prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
        end

        % Extract Structural Parameter Statistics
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

        % Save Elements to Structure array
        SO_phase_KDE_bias_difference(npower).phase_range{npower} = SO_thresholds{npower};
        SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
        SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

        % Calculate Post-Integration Summaries
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';
        SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

        % ----------- 3-Panel Diagnostics Panels (A, B, C) ----------------
        % A: Bias Performance vs Threshold
        nexttile((npower-1)*3 + 1); hold on;
        fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([thresholds, fliplr(thresholds)], [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
        ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
        xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
        title(sprintf('SO phase bin %s', SO_thresholds{npower}));
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12);

        % B: Volumetric Yield Rates
        nexttile((npower-1)*3 + 2); hold on;
        valid = isfinite(bias_mean) & isfinite(prop_mean);
        fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], [prop_mean(valid), fliplr(prop_mean(valid))], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Event Proportion vs. Bias Difference');
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12);

        % C: Structural Yield Variation Confidence Intervals
        nexttile((npower-1)*3 + 3); hold on;
        fill([bias_mean(valid), fliplr(bias_mean(valid))], [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Proportion vs. HPC Bias Difference');
        set(gca,"TickDir","out",'box','off','Color','none','FontSize',12);
        toc
    end
    save(fullfile(analysis_folder,'processed_data','SO_trough_KDE_bias_difference_based_on_V1_bias.mat'),'SO_phase_KDE_bias_difference');
else
    load(fullfile(analysis_folder,'processed_data','SO_trough_KDE_bias_difference_based_on_V1_bias.mat'),'SO_phase_KDE_bias_difference');
end

% -------- Comparative Line Profiles Mapping --------
clear Fill;
fig = figure('Name','KDE bias difference in HPC SO trough synchrony','Position',[640 100 2*1100/3 650/2]);
tiledlayout(1,2,'TileSpacing','compact');

subplot(1,2,1)
hold on
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    Fill(npower) = fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
end
bias_mean_shuf = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo_shuf = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi_shuf = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:);

Fill(3) = fill([thresholds, fliplr(thresholds)], [bias_CI_lo_shuf, fliplr(bias_CI_hi_shuf)], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean_shuf, 'Color','k', 'LineWidth', 2);
yline(0, '--r');
xlabel('V1 Bias threshold'); ylabel('HPC bias diff (T1 - T2)');
ylim([-0.1 0.35]); set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
legend(Fill([1 2 3]), {'SO trough unilateral', 'SO trough bilateral','Shuffled'}, 'box', 'off');


% -------- Categorical Integral Summary Bars --------
fig = figure('Name','KDE bias V1 AUC SO trough synchrony','Position',[640 100 281 325]);
data = SO_phase_KDE_bias_difference;
bar_width = 0.3; group_offset = 0.3; hold on;
clear BAR

for i = 1:2
    hold on
    % --- 1. Shuffled Baseline Metrics Execution ---
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


    % --- 2. Observed Sample Performance Profiles Execution ---
    x_real = i + group_offset;
    y_real = data(i).AUC_mean;
    neg_err_real = y_real - data(i).AUC_CI(1);
    pos_err_real = data(i).AUC_CI(2) - y_real;

    BAR(i+1) = bar(x_real, y_real, bar_width, 'FaceColor', colour_lines(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    errorbar(x_real, y_real, neg_err_real, pos_err_real, 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
end
hold off;
set(gca, 'XTick', 1:nBins, 'XTickLabel', {'Unilateral', 'Bilateral'});
xlim([0.5, nBins + 0.5]); ylabel('V1 bias AUC'); xlabel('Synchrony States');
legend(BAR, {'Shuffled', 'SO trough unilateral', 'SO trough bilateral'}, 'box', 'off');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12);

save_all_figures(fullfile(analysis_folder,'processed_data'),[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SO Trough Synchrony Binning (PRE-Ripple Window Analysis)

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;

if ~isfile(fullfile(analysis_folder,'processed_data','SO_trough_KDE_bias_difference_based_on_PRE_V1_bias.mat'))
    fig = figure('Name','KDE bias difference in HPC with different SO trough synchrony (PRE ripple)','Position',[640 100 1100 650]);
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        tic
        event_index = 1:length(track_bias_HC);

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

                if npower == 1 % trough unilateral
                    t1 = t1 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                    t2 = t2 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1); t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2), diff_tmp(i) = mean(t1_V1) - mean(t2_V1); end

                    total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                        sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);
                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

                    t1s = bb_shift >= th; t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                    t2s = t2s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1s); t2_V1 = boot_V1(t2s);
                    if any(t1s) && any(t2s), diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1); end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

                elseif npower == 2 % trough bilateral
                    t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                    t2 = t2 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1); t2_V1 = boot_V1(t2);
                    if any(t1) && any(t2), diff_tmp(i) = mean(t1_V1) - mean(t2_V1); end

                    total_events = sum((event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi));
                    prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

                    t1s = bb_shift >= th; t2s = bb_shift <= -th;
                    t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                    t2s = t2s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                    
                    t1_V1 = boot_V1(t1s); t2_V1 = boot_V1(t2s);
                    if any(t1s) && any(t2s), diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1); end
                    prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
                end
            end
            bias_diff_boot(iBoot, :) = diff_tmp;
            prop_events_boot(iBoot, :) = prop_tmp;
            bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
            prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
        end

        % Extract Non-parametric statistics summaries
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

        % Array allocation
        SO_phase_KDE_bias_difference(npower).phase_range{npower} = SO_thresholds{npower};
        SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
        SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';
        SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

        % ----------- 3-Panel Diagnostics Panels PRE ----------------
        nexttile((npower-1)*3 + 1); hold on;
        fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([thresholds, fliplr(thresholds)], [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
        ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
        xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
        title(sprintf('SO phase bin %s', SO_thresholds{npower})); set(gca,"TickDir","out",'box','off','Color','none','FontSize',12);

        nexttile((npower-1)*3 + 2); hold on;
        valid = isfinite(bias_mean) & isfinite(prop_mean);
        fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], [prop_mean(valid), fliplr(prop_mean(valid))], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Event Proportion vs. Bias Difference'); set(gca,"TickDir","out",'box','off','Color','none','FontSize',12);

        nexttile((npower-1)*3 + 3); hold on;
        fill([bias_mean(valid), fliplr(bias_mean(valid))], [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
        fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
        xlim([-0.1 0.35]); xline(0, '--r');
        xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
        title('Proportion vs. HPC Bias Difference'); set(gca,"TickDir","out",'box','off','Color','none','FontSize',12);
        toc
    end
    save(fullfile(analysis_folder,'processed_data','SO_trough_KDE_bias_difference_based_on_PRE_V1_bias.mat'),'SO_phase_KDE_bias_difference');
else
    load(fullfile(analysis_folder,'processed_data','SO_trough_KDE_bias_difference_based_on_PRE_V1_bias.mat'),'SO_phase_KDE_bias_difference');
end

% -------- Comparative Line Profiles Mapping PRE --------
clear Fill;
fig = figure('Name','KDE bias difference in HPC SO trough synchrony (PRE ripple)','Position',[640 100 2*1100/3 650/2]);
tiledlayout(1,2,'TileSpacing','compact');

subplot(1,2,1); hold on;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    Fill(npower) = fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
end
bias_mean_shuf = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo_shuf = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi_shuf = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:);

Fill(3) = fill([thresholds, fliplr(thresholds)], [bias_CI_lo_shuf, fliplr(bias_CI_hi_shuf)], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean_shuf, 'Color','k', 'LineWidth', 2);
yline(0, '--r');
xlabel('V1 Bias threshold'); ylabel('HPC bias diff (T1 - T2)');
ylim([-0.1 0.35]); set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
legend(Fill([1 2 3]), {'SO trough unilateral', 'SO trough bilateral','Shuffled'}, 'box', 'off');

% -------- Categorical Integral Summary Bars PRE --------
fig = figure('Name','KDE bias V1 AUC SO trough synchrony (PRE)','Position',[640 100 281 325]);
data = SO_phase_KDE_bias_difference;
bar_width = 0.3; group_offset = 0.15; hold on;
clear BAR

for i = 1:2
    % --- 1. Shuffled Baseline Metrics Execution ---
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

    x_real = i + group_offset;
    y_real = data(i).AUC_mean;
    neg_err_real = y_real - data(i).AUC_CI(1);
    pos_err_real = data(i).AUC_CI(2) - y_real;

    BAR(i+1) = bar(x_real, y_real, bar_width, 'FaceColor', colour_lines(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    errorbar(x_real, y_real, neg_err_real, pos_err_real, 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
end
hold off;
set(gca, 'XTick', 1:nBins, 'XTickLabel', {'Unilateral', 'Bilateral'});
xlim([0.5, nBins + 0.5]); ylabel('V1 bias AUC'); xlabel('Synchrony States');
legend(BAR, {'Shuffled', 'SO trough unilateral', 'SO trough bilateral'}, 'box', 'off');
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12);

save_all_figures(fullfile(analysis_folder,'processed_data'),[])