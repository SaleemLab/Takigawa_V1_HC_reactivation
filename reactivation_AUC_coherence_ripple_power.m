%% V1-HC reactivation AUC coherence analysis (Ripple power)
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




%% Temporal log odds AUC with different ripple power

fig = figure('Name','Ripple power distribution');
fig = histogram(ripple_info.ripple_power,5:0.2:18,'Normalization','probability');
xlabel('Ripple Power at ripple peak')
ylabel('Proportion of events')
set(gca,'TickDir','out','Box','off','FontSize',12);


% Ripple power binning
power_thresholds = prctile(ripple_info.ripple_power, 0:99.9/4:99.9);
nBins = length(power_thresholds) - 1;

% Time window parameters
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 20 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC selection window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;   % lightest
    226, 132, 187;   % mid-light
    212,  78, 156;   % mid-dark
    231,  41, 138    % darkest
    ] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

if isfile(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_ripple_power.mat'))==0;
    for t = 1:nTime
        t0 = time_bins(t)-win_size/2;
        t1 = time_bins(t) + win_size/2;

        % Sliding V1 window
        bins_to_select = bin_centers >= t0 & bin_centers < t1;

        fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

        for npower = 1:nBins
            % Ripple power bin index
            power_index = ripple_info.ripple_power > power_thresholds(npower) & ...
                ripple_info.ripple_power <= power_thresholds(npower+1);
            event_index = power_index > 0;
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
            AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
            AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
            AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
            AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
        end
    end

    save(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_ripple_power.mat'),'AUC')

else

    load(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_ripple_power.mat'))

end

% Plot temporal AUC traces
fig = figure('Name','Temporal V1 log-odds AUC different ripple powers','Position',[640 100 1100/3 900]);
tiledlayout(nBins, 1, 'TileSpacing','compact');

for npower = 1:nBins
    nexttile; hold on;
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

    % Shifted (black)
    fill([tvec fliplr(tvec)], ...
        [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
        [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0, '--r');
    xlabel('Time (s relative to ripple onset)');
    ylabel('V1 bias AUC');
    title(sprintf('Ripple power bin %d (%.2f–%.2f)', ...
        npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.2])
end



% Plot temporal AUC traces
fig = figure('Name','Temporal V1 log-odds AUC low vs high ripple powers','Position',[640 100 1100/3 900/4]);
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
    title('Ripple power low vs high');
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.2])
end

% Shifted (black)
fill([tvec fliplr(tvec)], ...
    [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

save_all_figures(fullfile(analysis_folder,'processed_data'),[])


%% Ripple power
power_thresholds = prctile(ripple_info.ripple_power,0:99.9/4:99.9);
nBins = length(power_thresholds) - 1;
bins_to_use = bin_centers>0 & bin_centers<0.1;
% bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
nBoot = 1000;

% Storage structure
ripple_power_KDE_bias_difference = struct;

% colour_lines = [158,202,225;33,113,181]/256;% two blue
% colour_lines = [158,202,225;107,174,214;66,146,198;33,113,181]/256;% two blue

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

if ~isfile(fullfile(analysis_folder,'processed_data','KDE_temporal_bias_ripple_power.mat'));

    % Plot layout
    fig = figure;
    fig.Position = [640 100 1100 650*2]
    fig.Name = 'KDE bias difference in V1 with different ripple powers';
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        % Ripple power bin index
        power_index = ripple_info.ripple_power > power_thresholds(npower) & ...
            ripple_info.ripple_power <= power_thresholds(npower+1);

        event_index = power_index >0;
        mean_bias = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
        mean_bias_V1 = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
        total_events = length(mean_bias);

        % Thresholds for bias
        thresholds = prctile(abs(mean_bias),0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = length(thresholds);

        % Bootstrap storage
        bias_diff_boot = NaN(nBoot, nThresh);
        prop_events_boot = NaN(nBoot, nThresh);
        % bias_diff_shifted_boot = NaN(nBoot, nThresh);
        % prop_events_shifted_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);


            boot_bias_shifted = mean_bias;
            boot_bias = mean_bias(idx);
            boot_V1 = mean_bias_V1(idx);

            diff_tmp = NaN(1, nThresh);
            prop_tmp = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);
                t1 = boot_bias >= th;
                t2 = boot_bias <= -th;

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);

                if ~isempty(t1_V1) && ~isempty(t2_V1)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            prop_events_boot(iBoot, :) = prop_tmp;
            %
            diff_tmp_shifted = NaN(1, nThresh);
            prop_tmp_shifted = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);
                t1 = boot_bias_shifted >= th;
                t2 = boot_bias_shifted <= -th;

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);

                if ~isempty(t1_V1) && ~isempty(t2_V1)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end

                prop_tmp_shifted(i) = (sum(t1) + sum(t2)) / total_events;
            end

            bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
            prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
        end

        % Compute stats
        bias_mean = mean(bias_diff_boot, 1, 'omitnan');
        bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
        bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);

        prop_mean = mean(prop_events_boot, 1, 'omitnan');
        prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
        prop_CI_hi = prctile(prop_events_boot, 97.5, 1);

        % Store results
        ripple_power_KDE_bias_difference(npower).power_range = [power_thresholds(npower), power_thresholds(npower+1)];
        ripple_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        ripple_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        ripple_power_KDE_bias_difference(npower).prop_mean = prop_mean;
        ripple_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        ripple_power_KDE_bias_difference(npower).thresholds = thresholds;


        % % Compute stats for shuffled (shifted) bias
        bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
        bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
        bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);

        prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
        prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
        prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

        % Store shifted (shuffled) results
        ripple_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        ripple_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        ripple_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

        % store AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

        ripple_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        ripple_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        ripple_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        ripple_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);
        % ---- Plot A: Bias difference vs. threshold ----
        nexttile((npower-1)*3 + 1);
        hold on;

        % Real
        x2 = [thresholds, fliplr(thresholds)];
        y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
        fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
        %
        % Time-shifted
        y_shift_lo = bias_shifted_CI_lo;
        y_shift_hi = bias_shifted_CI_hi;
        y2s = [y_shift_lo, fliplr(y_shift_hi)];
        fill(x2, y2s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);

        ylim([-0.15 0.35])
        xlim([0 1.4])
        yline(0,'--r')
        xlabel('HPC bias threshold');
        ylabel('V1 bias diff (T1 - T2)');
        title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
        %     grid on;
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        % ---- Plot B: Proportion vs. Bias Difference (X = bias, Y = proportion) ----
        nexttile((npower-1)*3 + 2);
        hold on;
        valid_idx = isfinite(bias_mean) & isfinite(prop_mean);
        x_vals = bias_mean(valid_idx);
        y_vals = prop_mean(valid_idx);
        x_lo = bias_CI_lo(valid_idx);
        x_hi = bias_CI_hi(valid_idx);
        x_shade = [x_lo, fliplr(x_hi)];
        y_shade = [y_vals, fliplr(y_vals)];

        % Real
        fill(x_shade, y_shade, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(x_vals, y_vals, '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);

        % Time-shifted
        x_vals_shift = bias_shifted_mean(valid_idx);
        y_vals_shift = prop_shifted_mean(valid_idx);
        x_lo_s = bias_shifted_CI_lo(valid_idx);
        x_hi_s = bias_shifted_CI_hi(valid_idx);
        x_shade_s = [x_lo_s, fliplr(x_hi_s)];
        y_shade_s = [y_vals_shift, fliplr(y_vals_shift)];

        fill(x_shade_s, y_shade_s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(x_vals_shift, y_vals_shift, 'k-', 'LineWidth', 1.5);

        xlim([-0.1 0.35])
        xline(0,'--r')
        xlabel('V1 bias diff (T1 - T2)');
        ylabel('Proportion of events detected');
        title('Event Proportion vs. Bias Difference');
        %     grid on;
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        % ---- Plot C: V1 Bias Difference vs. Proportion, shaded CI on Y ----
        nexttile((npower-1)*3 + 3);
        hold on;

        % Real
        valid_idx = isfinite(bias_mean) & isfinite(prop_mean);
        x_vals = bias_mean(valid_idx);
        y_vals = prop_mean(valid_idx);
        y_lo = prop_CI_lo(valid_idx);
        y_hi = prop_CI_hi(valid_idx);
        x_shade = [x_vals, fliplr(x_vals)];
        y_shade = [y_lo, fliplr(y_hi)];

        fill(x_shade, y_shade, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(x_vals, y_vals, '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);

        % Time-shifted
        x_vals_shift = bias_shifted_mean(valid_idx);
        y_vals_shift = prop_shifted_mean(valid_idx);
        y_lo_s = prop_shifted_CI_lo(valid_idx);
        y_hi_s = prop_shifted_CI_hi(valid_idx);
        x_shade_s = [x_vals_shift, fliplr(x_vals_shift)];
        y_shade_s = [y_lo_s, fliplr(y_hi_s)];

        fill(x_shade_s, y_shade_s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(x_vals_shift, y_vals_shift, 'k-', 'LineWidth', 1.5);

        xlim([-0.1 0.35])
        xline(0,'--r')
        xlabel('V1 bias diff (T1 - T2)');
        ylabel('Proportion of events detected');
        title('Proportion vs. V1 Bias Difference');
        %     grid on;
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    end

    % Save results
    save(fullfile(analysis_folder,'processed_data','ripple_power_KDE_bias_difference.mat'), 'ripple_power_KDE_bias_difference');
else
    load(fullfile(analysis_folder,'processed_data','ripple_power_KDE_bias_difference.mat'), 'ripple_power_KDE_bias_difference');
end



clear Fill
% Plot layout
fig = figure;
fig.Position = [640 100 2*1100/3 650/2]
fig.Name = 'KDE bias difference in V1 low vs high ripples';
% tiledlayout(nBins, 3, 'TileSpacing', 'compact');
% colour_lines = [158,202,225;33,113,181]/256;% two blue
nexttile
for npower = [1 4]
    bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_CI(2,:)

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color',colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    %     title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    ylim([-0.1 0.35])
    %     grid on;
end

bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);



yline(0,'--r')
legend(Fill([1 4 5]) ,{'Low ripple power','High ripple power','Shuffled'},'box','off')

nexttile
for npower = [1 4]
    bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_CI(2,:)
    prop_mean = ripple_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color',colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('Proportion of events detected');
    %     title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.1 0.35])
    %     grid on;
end

bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = ripple_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);

xline(0,'--r')
legend(Fill([1 4 5]) ,{'Low ripple power','High ripple power','Shuffled'},'box','off')


%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC low vs high ripples';
data = ripple_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.3;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:nBins
    if i == 1
        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(1).AUC_CI_shuffled(1);
        pos_err_shuf = data(1).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar

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
xlim([0.5, nBins + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','0-25','25-50','50-75','75-100'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'processed_data'),[]);


%% Ripple power (PRE)
power_thresholds = prctile(ripple_info.ripple_power,0:99.9/4:99.9);
nBins = length(power_thresholds) - 1;
bins_to_use = bin_centers>0 & bin_centers<0.1;
% bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
nBoot = 1000;

% Storage structure
ripple_power_bias_difference = struct;

% colour_lines = [158,202,225;33,113,181]/256;% two blue
% colour_lines = [158,202,225;107,174,214;66,146,198;33,113,181]/256;% two blue

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
    ] / 256;

if ~isfile(fullfile(analysis_folder,'processed_data','ripple_power_KDE_bias_difference_PRE.mat'))
    % Plot layout
    fig = figure;
    fig.Position = [640 100 1100 650*2]
    fig.Name = 'KDE bias difference in V1 with different ripple powers (PRE)';
    tiledlayout(nBins, 3, 'TileSpacing', 'compact');

    for npower = 1:nBins
        % Ripple power bin index
        power_index = ripple_info.ripple_power > power_thresholds(npower) & ...
            ripple_info.ripple_power <= power_thresholds(npower+1);

        event_index = power_index >0;
        mean_bias = mean(track_bias_HC(bins_to_use, event_index), 'omitnan');
        mean_bias_V1 = mean(track_bias_V1(bins_to_select, event_index), 'omitnan');
        total_events = length(mean_bias);

        % Thresholds for bias
        thresholds = prctile(abs(mean_bias),0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = length(thresholds);

        % Bootstrap storage
        bias_diff_boot = NaN(nBoot, nThresh);
        prop_events_boot = NaN(nBoot, nThresh);
        % bias_diff_shifted_boot = NaN(nBoot, nThresh);
        % prop_events_shifted_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1);


            boot_bias_shifted = mean_bias;
            boot_bias = mean_bias(idx);
            boot_V1 = mean_bias_V1(idx);

            diff_tmp = NaN(1, nThresh);
            prop_tmp = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);
                t1 = boot_bias >= th;
                t2 = boot_bias <= -th;

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);

                if ~isempty(t1_V1) && ~isempty(t2_V1)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            prop_events_boot(iBoot, :) = prop_tmp;
            %
            diff_tmp_shifted = NaN(1, nThresh);
            prop_tmp_shifted = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);
                t1 = boot_bias_shifted >= th;
                t2 = boot_bias_shifted <= -th;

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);

                if ~isempty(t1_V1) && ~isempty(t2_V1)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end

                prop_tmp_shifted(i) = (sum(t1) + sum(t2)) / total_events;
            end

            bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
            prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
        end

        % Compute stats
        bias_mean = mean(bias_diff_boot, 1, 'omitnan');
        bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
        bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);

        prop_mean = mean(prop_events_boot, 1, 'omitnan');
        prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
        prop_CI_hi = prctile(prop_events_boot, 97.5, 1);

        % Store results
        ripple_power_KDE_bias_difference(npower).power_range = [power_thresholds(npower), power_thresholds(npower+1)];
        ripple_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
        ripple_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
        ripple_power_KDE_bias_difference(npower).prop_mean = prop_mean;
        ripple_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
        ripple_power_KDE_bias_difference(npower).thresholds = thresholds;


        % % Compute stats for shuffled (shifted) bias
        bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
        bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
        bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);

        prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
        prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
        prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

        % Store shifted (shuffled) results
        ripple_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
        ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
        ripple_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
        ripple_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

        % store AUC
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

        ripple_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
        ripple_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
        ripple_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
        ripple_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

        % ---- Plot A: Bias difference vs. threshold ----
        nexttile((npower-1)*3 + 1);
        hold on;

        % Real
        x2 = [thresholds, fliplr(thresholds)];
        y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
        fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
        %
        % Time-shifted
        y_shift_lo = bias_shifted_CI_lo;
        y_shift_hi = bias_shifted_CI_hi;
        y2s = [y_shift_lo, fliplr(y_shift_hi)];
        fill(x2, y2s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);

        ylim([-0.15 0.35])
        xlim([0 1.4])
        yline(0,'--r')
        xlabel('HPC bias threshold');
        ylabel('V1 bias diff (T1 - T2)');
        title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
        %     grid on;
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        % ---- Plot B: Proportion vs. Bias Difference (X = bias, Y = proportion) ----
        nexttile((npower-1)*3 + 2);
        hold on;
        valid_idx = isfinite(bias_mean) & isfinite(prop_mean);
        x_vals = bias_mean(valid_idx);
        y_vals = prop_mean(valid_idx);
        x_lo = bias_CI_lo(valid_idx);
        x_hi = bias_CI_hi(valid_idx);
        x_shade = [x_lo, fliplr(x_hi)];
        y_shade = [y_vals, fliplr(y_vals)];

        % Real
        fill(x_shade, y_shade, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(x_vals, y_vals, '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);

        % Time-shifted
        x_vals_shift = bias_shifted_mean(valid_idx);
        y_vals_shift = prop_shifted_mean(valid_idx);
        x_lo_s = bias_shifted_CI_lo(valid_idx);
        x_hi_s = bias_shifted_CI_hi(valid_idx);
        x_shade_s = [x_lo_s, fliplr(x_hi_s)];
        y_shade_s = [y_vals_shift, fliplr(y_vals_shift)];

        fill(x_shade_s, y_shade_s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(x_vals_shift, y_vals_shift, 'k-', 'LineWidth', 1.5);

        xlim([-0.1 0.35])
        xline(0,'--r')
        xlabel('V1 bias diff (T1 - T2)');
        ylabel('Proportion of events detected');
        title('Event Proportion vs. Bias Difference');
        %     grid on;
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        % ---- Plot C: V1 Bias Difference vs. Proportion, shaded CI on Y ----
        nexttile((npower-1)*3 + 3);
        hold on;

        % Real
        valid_idx = isfinite(bias_mean) & isfinite(prop_mean);
        x_vals = bias_mean(valid_idx);
        y_vals = prop_mean(valid_idx);
        y_lo = prop_CI_lo(valid_idx);
        y_hi = prop_CI_hi(valid_idx);
        x_shade = [x_vals, fliplr(x_vals)];
        y_shade = [y_lo, fliplr(y_hi)];

        fill(x_shade, y_shade, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(x_vals, y_vals, '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);

        % Time-shifted
        x_vals_shift = bias_shifted_mean(valid_idx);
        y_vals_shift = prop_shifted_mean(valid_idx);
        y_lo_s = prop_shifted_CI_lo(valid_idx);
        y_hi_s = prop_shifted_CI_hi(valid_idx);
        x_shade_s = [x_vals_shift, fliplr(x_vals_shift)];
        y_shade_s = [y_lo_s, fliplr(y_hi_s)];

        fill(x_shade_s, y_shade_s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(x_vals_shift, y_vals_shift, 'k-', 'LineWidth', 1.5);

        xlim([-0.1 0.35])
        xline(0,'--r')
        xlabel('V1 bias diff (T1 - T2)');
        ylabel('Proportion of events detected');
        title('Proportion vs. V1 Bias Difference');
        %     grid on;
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
    save(fullfile(analysis_folder,'processed_data','ripple_power_KDE_bias_difference_PRE.mat'),'ripple_power_KDE_bias_difference')
else
    load(fullfile(analysis_folder,'processed_data','ripple_power_KDE_bias_difference_PRE.mat'),'ripple_power_KDE_bias_difference')
end

clear Fill
% Plot layout
fig = figure;
fig.Position = [640 100 2*1100/3 650/2]
fig.Name = 'KDE bias difference in V1 low vs high ripples (PRE)';
% tiledlayout(nBins, 3, 'TileSpacing', 'compact');
% colour_lines = [158,202,225;33,113,181]/256;% two blue
nexttile
for npower = [1 4]
    bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_CI(2,:)

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color',colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    %     title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    ylim([-0.1 0.35])
    %     grid on;
end

bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);



yline(0,'--r')
legend(Fill([1 4 5]) ,{'Low ripple power','High ripple power','Shuffled'},'box','off')

nexttile
for npower = [1 4]
    bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_CI(2,:)
    prop_mean = ripple_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color',colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('Proportion of events detected');
    %     title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    xlim([-0.1 0.35])
    %     grid on;
end


bias_mean = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = ripple_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = ripple_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0,'--r')
legend(Fill([1 4 5]) ,{'Low ripple power','High ripple power','Shuffled'},'box','off')


%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC low vs high ripples (PRE)';
data = ripple_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.3;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:nBins
    if i == 1
        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(1).AUC_CI_shuffled(1);
        pos_err_shuf = data(1).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');

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
xlim([0.5, nBins + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','0-25','25-50','50-75','75-100'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

save_all_figures(fullfile(analysis_folder,'processed_data'),[])
