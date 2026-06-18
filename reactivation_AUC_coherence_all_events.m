%% V1-HC reactivation AUC coherence analysis (All events)

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

%% Plot HC track bias distribution
% % event_index = [T1_events T2_events];
% % z_event = nanmean(z_bias(bins_to_use,[T1_events T2_events]));
% 
track1_color = [0 0 1]; % blue
track2_color = [1 0 0]; % red

n_half = 80;     % Number of steps from red to white and white to blue
n_white = 60;     % Number of extra white steps to make the white region broader

% Red to white
r2w = [linspace(1,1,n_half)', linspace(0,1,n_half)', linspace(0,1,n_half)'];
% White to blue
w2b = [linspace(1,0,n_half)', linspace(1,0,n_half)', linspace(1,1,n_half)'];
% Extra white section
white = ones(n_white, 3);
% Combine
red_white_blue = [r2w; white; w2b];
% 

nfig = figure;
% nfig.Name = 'KDE Bias PSTH all events (no thresholding blue red white)'; 
nfig.Name = 'HC track bias PSTH all events'; 
bins_to_use = bin_centers>0 & bin_centers<0.1;

event_index = 1:length(track_bias_HC);
z_event = nanmean(track_bias_HC(bins_to_use,event_index));
[~,sorted_index] = sort(z_event,'descend');

h = imagesc(bin_centers,[],track_bias_HC(:,event_index(sorted_index))');clim([-3 3])
yline([0.75*length(track_bias_HC)],'r','LineWidth',1)
yline([0.25*length(track_bias_HC)],'b','LineWidth',1)
set(h, 'AlphaData', ~isnan(track_bias_HC(:, event_index(sorted_index))'));  % Hide NaNs (make them transparent)
xlim([-0.5 0.5])
colorbar;
colormap((red_white_blue))
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
% imagesc(z_bias_V1);colorbar;
xlabel('Time (s)')
ylabel('Ripple event')


%% Temporal log odds AUC for all ripple events

% Time window parameters
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 20 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC selection window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme (Using the darkest pink from your original scheme for all events)
colour_line = [231,  41, 138] / 256; 

% Storage (Simplified to 1D over time)
AUC.mean = nan(nTime, 1);
AUC.ci = nan(nTime, 2);
AUC.shifted_mean = nan(nTime, 1);
AUC.shifted_ci = nan(nTime, 2);

if isfile(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_all_ripples.mat'))==0;

    for t = 1:nTime
        t0 = time_bins(t) - win_size/2;
        t1 = time_bins(t) + win_size/2;

        % Sliding V1 window
        bins_to_select = bin_centers >= t0 & bin_centers < t1;
        fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s) [All Ripples]\n', t0, t1);

        % Select ALL ripple events
        event_index = true(size(ripple_info.ripple_power));
        total_events = sum(event_index);

        if total_events < 10, continue; end

        % Mean log-odds
        mean_bias = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');        % HPC fixed
        mean_bias_V1 = mean(track_bias_V1(bins_to_select, event_index), 'omitnan'); % V1 sliding

        % Quantile thresholds based on |HPC bias|
        thresholds = prctile(abs(mean_bias), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        % Bootstrap arrays
        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1); % resample event IDs

            boot_HPC = mean_bias(idx);
            boot_V1 = mean_bias_V1(idx);
            boot_bias_shifted = mean_bias; % event id randomised

            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Real HPC–V1 pairing
                t1 = boot_HPC >= th;
                t2 = boot_HPC <= -th;
                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if ~isempty(t1_V1) && ~isempty(t2_V1)
                    diff_tmp(i) = mean(t1_V1, 'omitnan') - mean(t2_V1, 'omitnan');
                end

                % Shifted pairing (HPC–V1 decoupled)
                t1s = boot_bias_shifted >= th;
                t2s = boot_bias_shifted <= -th;
                t1_V1s = boot_V1(t1s);
                t2_V1s = boot_V1(t2s);
                if ~isempty(t1_V1s) && ~isempty(t2_V1s)
                    diff_tmp_shift(i) = mean(t1_V1s, 'omitnan') - mean(t2_V1s, 'omitnan');
                end
            end
            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC (mean across thresholds)
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t) = mean(auc_boot, 'omitnan');
        AUC.ci(t, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end

    %%% Save Data
    save_all_figures(fullfile(analysis_folder,'processed_data'),[]);
    save(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_all_ripples.mat'), 'AUC');
else
    load(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_all_ripples.mat'), 'AUC');
end

%%% Plot Temporal AUC Trace (All Events vs Shuffled)
clear p
fig = figure('Name', 'Temporal V1 log-odds AUC All Ripples', 'Position', [640 100 1100/3 900/4]);
% fig = figure('Name','Temporal V1 log-odds AUC low vs high ripple powers','Position',[640 100 1100/3 900/4]);
tiledlayout(4, 1, 'TileSpacing','compact');
hold on;

m  = AUC.mean;
ci = AUC.ci;
m_shift  = AUC.shifted_mean;
ci_shift = AUC.shifted_ci;
tvec = time_bins;

% Remove any NaN times where there weren't enough events
valid_idx = ~isnan(m);
tvec = tvec(valid_idx);
m = m(valid_idx);
ci = ci(valid_idx, :);
m_shift = m_shift(valid_idx);
ci_shift = ci_shift(valid_idx, :);

% Plot Shifted / Shuffled Baseline (Black)
fill([tvec, fliplr(tvec)], ...
     [ci_shift(:,1)', fliplr(ci_shift(:,2)')], ...
     [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
p(1)=plot(tvec, m_shift, 'k', 'LineWidth', 1.2, 'DisplayName', 'Shuffled Baseline');

% Plot Real Data (Coloured Pink)
fill([tvec, fliplr(tvec)], ...
     [ci(:,1)', fliplr(ci(:,2)')], ...
     colour_line, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
p(2)=plot(tvec, m, 'Color', colour_line, 'LineWidth', 2, 'DisplayName', 'Real Data');

% Formatting
yline(0, '--r', 'LineWidth', 1);
xlabel('Time (s relative to ripple onset)');
ylabel('V1 bias AUC');
title('Temporal V1 Log-Odds AUC (All Ripple Events)');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-0.5 0.5]);
ylim([-0.1 0.2]);
legend([p(1:2)],'Shuffle','All events','Box', 'off');

% Save Figures
save_all_figures(fullfile(analysis_folder,'processed_data'),[]);


%% All ripples (POST)
% Selection windows
bins_to_use = bin_centers > 0 & bin_centers < 0.1;
bins_to_select = bin_centers > 0 & bin_centers < 0.2;
nBoot = 1000;

% Color scheme (Dark pink from your original palette used for the combined data)
main_color = [231, 41, 138] / 256; 

%%% Process All Ripple Events Combined
% Use all events by assigning a true mask everywhere
event_index = true(size(ripple_info.ripple_power));

mean_bias = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
mean_bias_V1 = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
total_events = length(mean_bias);

% Thresholds for bias
thresholds = prctile(abs(mean_bias), 0:10:100);
thresholds = thresholds(1:end-1);
nThresh = length(thresholds);

% Bootstrap storage
bias_diff_boot = NaN(nBoot, nThresh);
prop_events_boot = NaN(nBoot, nThresh);
bias_diff_shifted_boot = NaN(nBoot, nThresh);
prop_events_shifted_boot = NaN(nBoot, nThresh);

if isfile(fullfile(analysis_folder,'processed_data','all_ripples_KDE_bias_difference.mat'))==0;

    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, total_events, total_events, 1);

        boot_bias_shifted = mean_bias; % Real baseline template
        boot_bias = mean_bias(idx);
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);

        for i = 1:nThresh
            th = thresholds(i);

            % --- Real Analysis ---
            t1 = boot_bias >= th;
            t2 = boot_bias <= -th;
            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if ~isempty(t1_V1) && ~isempty(t2_V1)
                diff_tmp(i) = mean(t1_V1, 'omitnan') - mean(t2_V1, 'omitnan');
            end
            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

            % --- Shuffled Analysis ---
            t1s = boot_bias_shifted >= th;
            t2s = boot_bias_shifted <= -th;
            t1_V1s = boot_V1(t1s);
            t2_V1s = boot_V1(t2s);
            if ~isempty(t1_V1s) && ~isempty(t2_V1s)
                diff_tmp_shifted(i) = mean(t1_V1s, 'omitnan') - mean(t2_V1s, 'omitnan');
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    %%% Compute Statistics
    all_ripples_res.bias_diff_mean = mean(bias_diff_boot, 1, 'omitnan');
    all_ripples_res.bias_diff_CI = [prctile(bias_diff_boot, 2.5, 1); prctile(bias_diff_boot, 97.5, 1)];
    all_ripples_res.prop_mean = mean(prop_events_boot, 1, 'omitnan');
    all_ripples_res.prop_CI = [prctile(prop_events_boot, 2.5, 1); prctile(prop_events_boot, 97.5, 1)];
    all_ripples_res.thresholds = thresholds;

    all_ripples_res.bias_diff_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    all_ripples_res.bias_diff_shifted_CI = [prctile(bias_diff_shifted_boot, 2.5, 1); prctile(bias_diff_shifted_boot, 97.5, 1)];
    all_ripples_res.prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    all_ripples_res.prop_shifted_CI = [prctile(prop_events_shifted_boot, 2.5, 1); prctile(prop_events_shifted_boot, 97.5, 1)];

    % Area Under Curve (AUC) calculations
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds) - min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds) - min(thresholds)))';

    all_ripples_res.AUC_mean = mean(auc_boot, 'omitnan');
    all_ripples_res.AUC_CI = prctile(auc_boot, [2.5 97.5]);
    all_ripples_res.AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    all_ripples_res.AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % Save calculation results
    save(fullfile(analysis_folder,'processed_data','all_ripples_KDE_bias_difference.mat'), 'all_ripples_res');

else
    load(fullfile(analysis_folder,'processed_data','all_ripples_KDE_bias_difference.mat'), 'all_ripples_res');
end

%%% Plotting PLS-KDE track bias difference in V1 based on HC track bias
%%% selection
fig1 = figure('Name', 'KDE bias difference in V1 - All Ripples', 'Position', [640 100 1100 350]);
tiledlayout(1, 3, 'TileSpacing', 'compact');

% Setup local variables for easy plotting syntax
b_mean  = all_ripples_res.bias_diff_mean;
b_CI_lo = all_ripples_res.bias_diff_CI(1,:);
b_CI_hi = all_ripples_res.bias_diff_CI(2,:);
p_mean  = all_ripples_res.prop_mean;
p_CI_lo = all_ripples_res.prop_CI(1,:);
p_CI_hi = all_ripples_res.prop_CI(2,:);

bs_mean  = all_ripples_res.bias_diff_shifted_mean;
bs_CI_lo = all_ripples_res.bias_diff_shifted_CI(1,:);
bs_CI_hi = all_ripples_res.bias_diff_shifted_CI(2,:);
ps_mean  = all_ripples_res.prop_shifted_mean;
ps_CI_lo = all_ripples_res.prop_shifted_CI(1,:);
ps_CI_hi = all_ripples_res.prop_shifted_CI(2,:);

% ---- Panel A: Bias difference vs. Threshold ----
clear p
nexttile; hold on;
fill([thresholds, fliplr(thresholds)], [b_CI_lo, fliplr(b_CI_hi)], main_color, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
p(1) = plot(thresholds, b_mean, 'Color', main_color, 'LineWidth', 2);
fill([thresholds, fliplr(thresholds)], [bs_CI_lo, fliplr(bs_CI_hi)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
p(2) = plot(thresholds, bs_mean, 'k-', 'LineWidth', 1.5);
ylim([-0.15 0.35]); xlim([0 1.4]); yline(0, '--r');
xlabel('HPC bias threshold'); ylabel('V1 bias diff (T1 - T2)');
title('All Ripples Combined');
set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);
legend([p(1:2)],'Shuffle','All events','Box', 'off');

% ---- Panel B: Proportion vs. Bias Difference (Shaded on X-axis) ----
nexttile; hold on;
valid_idx = isfinite(b_mean) & isfinite(p_mean);
fill([b_CI_lo(valid_idx), fliplr(b_CI_hi(valid_idx))], [p_mean(valid_idx), fliplr(p_mean(valid_idx))], main_color, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(b_mean(valid_idx), p_mean(valid_idx), '-', 'Color', main_color, 'LineWidth', 2);

fill([bs_CI_lo(valid_idx), fliplr(bs_CI_hi(valid_idx))], [ps_mean(valid_idx), fliplr(ps_mean(valid_idx))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot(bs_mean(valid_idx), ps_mean(valid_idx), 'k-', 'LineWidth', 1.5);
xlim([-0.1 0.35]); xline(0, '--r');
xlabel('V1 bias diff (T1 - T2)'); ylabel('Proportion of events detected');
title('Event Proportion vs. Bias Difference');
set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);
nexttile; 

%%% Figure 2: AUC Single Group Bar Chart
fig2 = figure('Name', 'KDE bias V1 AUC All Ripples', 'Position', [640 100 240 325]);
hold on;

bar_width = 0.4;
% Shuffled Bar
y_shuf = all_ripples_res.AUC_mean_shuffled;
neg_err_shuf = y_shuf - all_ripples_res.AUC_CI_shuffled(1);
pos_err_shuf = all_ripples_res.AUC_CI_shuffled(2) - y_shuf;

h_shuf = bar(1, y_shuf, bar_width, 'FaceColor', 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
errorbar(1, y_shuf, neg_err_shuf, pos_err_shuf, 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');

% Real Bar
y_real = all_ripples_res.AUC_mean;
neg_err_real = y_real - all_ripples_res.AUC_CI(1);
pos_err_real = all_ripples_res.AUC_CI(2) - y_real;

h_real = bar(2, y_real, bar_width, 'FaceColor', main_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
errorbar(2, y_real, neg_err_real, pos_err_real, 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');

hold off;

% Configuration
set(gca, 'XTick', [1 2], 'XTickLabel', {'Shuffled', 'All Ripples'});
xlim([0.3, 2.7]);
ylim([-0.02 0.15])
ylabel('V1 bias AUC');
set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);


%% All ripples (PRE)
% Selection windows
bins_to_use = bin_centers > 0 & bin_centers < 0.1;
bins_to_select = bin_centers > -0.2 & bin_centers < 0;
% bins_to_select = bin_centers > -0.05 & bin_centers < 0.05;

nBoot = 1000;

% Color scheme (Dark pink from your original palette used for the combined data)
main_color = [231, 41, 138] / 256; 

%%% Process All Ripple Events Combined
% Use all events by assigning a true mask everywhere
event_index = true(size(ripple_info.ripple_power));

mean_bias = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
mean_bias_V1 = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
total_events = length(mean_bias);

% Thresholds for bias
thresholds = prctile(abs(mean_bias), 0:10:100);
thresholds = thresholds(1:end-1);
nThresh = length(thresholds);

% Bootstrap storage
bias_diff_boot = NaN(nBoot, nThresh);
prop_events_boot = NaN(nBoot, nThresh);
bias_diff_shifted_boot = NaN(nBoot, nThresh);
prop_events_shifted_boot = NaN(nBoot, nThresh);

if isfile(fullfile(analysis_folder,'processed_data','all_ripples_KDE_bias_difference.mat'))==0;
    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, total_events, total_events, 1);

        boot_bias_shifted = mean_bias; % Real baseline template
        boot_bias = mean_bias(idx);
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);

        for i = 1:nThresh
            th = thresholds(i);

            % --- Real Analysis ---
            t1 = boot_bias >= th;
            t2 = boot_bias <= -th;
            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if ~isempty(t1_V1) && ~isempty(t2_V1)
                diff_tmp(i) = mean(t1_V1, 'omitnan') - mean(t2_V1, 'omitnan');
            end
            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

            % --- Shuffled Analysis ---
            t1s = boot_bias_shifted >= th;
            t2s = boot_bias_shifted <= -th;
            t1_V1s = boot_V1(t1s);
            t2_V1s = boot_V1(t2s);
            if ~isempty(t1_V1s) && ~isempty(t2_V1s)
                diff_tmp_shifted(i) = mean(t1_V1s, 'omitnan') - mean(t2_V1s, 'omitnan');
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    %%% Compute Statistics
    all_ripples_res.bias_diff_mean = mean(bias_diff_boot, 1, 'omitnan');
    all_ripples_res.bias_diff_CI = [prctile(bias_diff_boot, 2.5, 1); prctile(bias_diff_boot, 97.5, 1)];
    all_ripples_res.prop_mean = mean(prop_events_boot, 1, 'omitnan');
    all_ripples_res.prop_CI = [prctile(prop_events_boot, 2.5, 1); prctile(prop_events_boot, 97.5, 1)];
    all_ripples_res.thresholds = thresholds;

    all_ripples_res.bias_diff_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    all_ripples_res.bias_diff_shifted_CI = [prctile(bias_diff_shifted_boot, 2.5, 1); prctile(bias_diff_shifted_boot, 97.5, 1)];
    all_ripples_res.prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    all_ripples_res.prop_shifted_CI = [prctile(prop_events_shifted_boot, 2.5, 1); prctile(prop_events_shifted_boot, 97.5, 1)];

    % Area Under Curve (AUC) calculations
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds) - min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds) - min(thresholds)))';

    all_ripples_res.AUC_mean = mean(auc_boot, 'omitnan');
    all_ripples_res.AUC_CI = prctile(auc_boot, [2.5 97.5]);
    all_ripples_res.AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    all_ripples_res.AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % Save calculation results
    save(fullfile(analysis_folder,'processed_data','all_ripples_KDE_bias_difference_PRE.mat'), 'all_ripples_res');
else
    load(fullfile(analysis_folder,'processed_data','all_ripples_KDE_bias_difference_PRE.mat'), 'all_ripples_res');
end

%%% Figure 1: Multi-metric Diagnosis Panel (Single Row)
% fig1 = figure('Name', 'KDE bias difference in V1 - All Ripples 100ms window', 'Position', [640 100 1100 350]);
fig1 = figure('Name', 'KDE bias difference in V1 - All Ripples PRE', 'Position', [640 100 1100 350]);
tiledlayout(1, 3, 'TileSpacing', 'compact');

% Setup local variables for easy plotting syntax
b_mean  = all_ripples_res.bias_diff_mean;
b_CI_lo = all_ripples_res.bias_diff_CI(1,:);
b_CI_hi = all_ripples_res.bias_diff_CI(2,:);
p_mean  = all_ripples_res.prop_mean;
p_CI_lo = all_ripples_res.prop_CI(1,:);
p_CI_hi = all_ripples_res.prop_CI(2,:);

bs_mean  = all_ripples_res.bias_diff_shifted_mean;
bs_CI_lo = all_ripples_res.bias_diff_shifted_CI(1,:);
bs_CI_hi = all_ripples_res.bias_diff_shifted_CI(2,:);
ps_mean  = all_ripples_res.prop_shifted_mean;
ps_CI_lo = all_ripples_res.prop_shifted_CI(1,:);
ps_CI_hi = all_ripples_res.prop_shifted_CI(2,:);

% ---- Panel A: Bias difference vs. Threshold ----
clear p
nexttile; hold on;
fill([thresholds, fliplr(thresholds)], [b_CI_lo, fliplr(b_CI_hi)], main_color, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
p(1) = plot(thresholds, b_mean, 'Color', main_color, 'LineWidth', 2);
fill([thresholds, fliplr(thresholds)], [bs_CI_lo, fliplr(bs_CI_hi)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
p(2) = plot(thresholds, bs_mean, 'k-', 'LineWidth', 1.5);
ylim([-0.15 0.35]); xlim([0 1.4]); yline(0, '--r');
xlabel('HPC bias threshold'); ylabel('V1 bias diff (T1 - T2)');
title('All Ripples Combined');
set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);
legend([p(1:2)],'Shuffle','All events','Box', 'off');

% ---- Panel B: Proportion vs. Bias Difference (Shaded on X-axis) ----
nexttile; hold on;
valid_idx = isfinite(b_mean) & isfinite(p_mean);
fill([b_CI_lo(valid_idx), fliplr(b_CI_hi(valid_idx))], [p_mean(valid_idx), fliplr(p_mean(valid_idx))], main_color, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(b_mean(valid_idx), p_mean(valid_idx), '-', 'Color', main_color, 'LineWidth', 2);

fill([bs_CI_lo(valid_idx), fliplr(bs_CI_hi(valid_idx))], [ps_mean(valid_idx), fliplr(ps_mean(valid_idx))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot(bs_mean(valid_idx), ps_mean(valid_idx), 'k-', 'LineWidth', 1.5);
xlim([-0.1 0.35]); xline(0, '--r');
xlabel('V1 bias diff (T1 - T2)'); ylabel('Proportion of events detected');
title('Event Proportion vs. Bias Difference');
set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);

% ---- Panel C: Proportion vs. V1 Bias Difference (Shaded on Y-axis) ----
nexttile; hold on;
fill([b_mean(valid_idx), fliplr(b_mean(valid_idx))], [p_CI_lo(valid_idx), fliplr(p_CI_hi(valid_idx))], main_color, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(b_mean(valid_idx), p_mean(valid_idx), '-', 'Color', main_color, 'LineWidth', 2);

fill([bs_mean(valid_idx), fliplr(bs_mean(valid_idx))], [ps_CI_lo(valid_idx), fliplr(ps_CI_hi(valid_idx))], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot(bs_mean(valid_idx), ps_mean(valid_idx), 'k-', 'LineWidth', 1.5);
xlim([-0.1 0.35]); xline(0, '--r');
xlabel('V1 bias diff (T1 - T2)'); ylabel('Proportion of events detected');
title('Proportion vs. V1 Bias Difference');
set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);


%%% Figure 2: AUC Single Group Bar Chart
fig2 = figure('Name', 'KDE bias V1 AUC All Ripples PRE', 'Position', [640 100 240 325]);
% fig2 = figure('Name', 'KDE bias V1 AUC All Ripples 100ms window', 'Position', [640 100 240 325]);

hold on;

bar_width = 0.4;
% Shuffled Bar
y_shuf = all_ripples_res.AUC_mean_shuffled;
neg_err_shuf = y_shuf - all_ripples_res.AUC_CI_shuffled(1);
pos_err_shuf = all_ripples_res.AUC_CI_shuffled(2) - y_shuf;

h_shuf = bar(1, y_shuf, bar_width, 'FaceColor', 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
errorbar(1, y_shuf, neg_err_shuf, pos_err_shuf, 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');

% Real Bar
y_real = all_ripples_res.AUC_mean;
neg_err_real = y_real - all_ripples_res.AUC_CI(1);
pos_err_real = all_ripples_res.AUC_CI(2) - y_real;

h_real = bar(2, y_real, bar_width, 'FaceColor', main_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
errorbar(2, y_real, neg_err_real, pos_err_real, 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');

hold off;

% Configuration
set(gca, 'XTick', [1 2], 'XTickLabel', {'Shuffled', 'All Ripples'});
xlim([0.3, 2.7]);
ylim([-0.02 0.15])
ylabel('V1 bias AUC');
set(gca, 'TickDir', 'out', 'box', 'off', 'Color', 'none', 'FontSize', 12);


save_all_figures(fullfile(analysis_folder,'processed_data'),[])
