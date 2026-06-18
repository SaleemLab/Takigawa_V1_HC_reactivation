%% Plotting ROC AUC of RUN 2 track discrimination 
analysis_folder = pwd;

%%% Load PLS ROC RUN with 10 fold validation using only bins with above-chance
%%% discriminiability
% load(fullfile(analysis_folder,'processed_data','PLS_RUN_validations.mat'),'PLS_RUN_validations_all')

%%% Or Load PLS ROC RUN with 10 fold validation using all bins
load(fullfile(analysis_folder,'processed_data','PLS_RUN_validations_all_bins.mat'),'PLS_RUN_validations_all')

%%%%% ROC track discrimination based on PLS
timebins = [20,100];
% title_text = 'PLS contextual discrimination HPC RUN1 scatter (95th percentile)';
title_text = 'PLS contextual discrimination HPC RUN1 scatter (all bins 95th percentile)';

fig = figure('Name', title_text, 'Position', [200 100 640 580]); hold on;

line_w = 0.2;
%%% n = 1 for 20ms timebin and n = 2 for 100ms timebin ROC
for n = 1:2
    nexttile
    % Extract base data
    fpr = PLS_RUN_validations_all.FPR(1,:);
    % These should be (sessions x fpr_points)
    data_real = squeeze(PLS_RUN_validations_all.HPC_TPR(n,:,:));
    data_shuf = squeeze(PLS_RUN_validations_all.HPC_TPR_shuffled(n,:,:));
    
    % Calculate Means
    TPR_real_mean = mean(data_real, 1, 'omitnan');
    TPR_shuf_mean = mean(data_shuf, 1, 'omitnan');
    
    hold on
    % 1. Plot the diagonal chance line
    plot([0 1],[0 1],'--k', 'HandleVisibility', 'off')
    
    % 2. Plot individual session lines (Low Alpha)
    % Color definitions
    real_col = [231,41,138]/255;
    shuf_col = [0, 0, 0]/255;
    
    % Plot all individual Real lines
    p_indiv_real = plot(fpr, data_real', 'Color', [real_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % Plot all individual Shuffled lines
    p_indiv_shuf = plot(fpr, data_shuf', 'Color', [shuf_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % 3. Plot Mean lines (Thick)
    PLOT(1) = plot(fpr, TPR_real_mean, 'Color', real_col, 'LineWidth', 2.5);
    PLOT(2) = plot(fpr, TPR_shuf_mean, 'Color', shuf_col, 'LineWidth', 2.5);
    
    % Formatting
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    title(sprintf('ROC curve HPC %i ms bins', timebins(n)));
    
    legend(PLOT(1:2), {'Real Mean', 'Shuffle Mean'}, 'Location', 'southeast', 'Box', 'off');
    
    xlim([0 1]); ylim([0 1]);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end



for n = 1:2
    nexttile
    % Extract data for the current condition
    real_data = PLS_RUN_validations_all.HPC_AUC(n,:);
    shuffled_data = PLS_RUN_validations_all.HPC_AUC_shuffled(n,:);
    
    % Combine for easier mean calculation
    % AUC = [real_data; shuffled_data];
    mean_real = mean(real_data, 'omitnan');
    mean_shuf = mean(shuffled_data, 'omitnan');
    bar_colors = [231,41,138; 0, 0, 0]/255;
    x_pos = [1 2];

    hold on

    % 1. Plot the Bars (Mean only, no error bars)
    plot([x_pos(1)-line_w, x_pos(1)+line_w], [mean_real, mean_real], 'Color', bar_colors(1,:), 'LineWidth', 4);
    plot([x_pos(2)-line_w, x_pos(2)+line_w], [mean_shuf, mean_shuf], 'Color', bar_colors(2,:), 'LineWidth', 4);

    % 2. Plot Raw Data Points and Connections with Jitter
    rng(1); % Seed for consistent jitter appearance
    jitter_strength = 0.2; % Adjusted for better visibility of distributions
    
    for s = 1:length(real_data)
        if ~isnan(real_data(s)) && ~isnan(shuffled_data(s))
            % Calculate jittered x-positions
            % This spreads the dots randomly around the x_pos center
            x1 = x_pos(1) + (rand - 0.5) * jitter_strength;
            x2 = x_pos(2) + (rand - 0.5) * jitter_strength;
            
            % Plot the connection line first (so it's behind the dots)
            plot([x1, x2], [real_data(s), shuffled_data(s)], 'Color', [0.8 0.8 0.8 0.4], 'LineWidth', 0.5);
            
            % Plot the individual dots
            scatter(x1, real_data(s), 50, bar_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.3);
            scatter(x2, shuffled_data(s), 50, bar_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.3);
        end
    end

    % Statistics and Formatting
    [p, h] = signrank(real_data, shuffled_data, 'tail', 'right');
    
    xlim([-0.5 3.5])
    xticks([1 2]);
    xticklabels({'Real', 'Shuffled'});
    ylabel('AUC');
    ylim([0 1.1]);
    yticks([0:0.25:1])
    title(['Condition ' num2str(n) ' (p=' num2str(p, '%.4e') ')']);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end


%%%%%%%%% V1
timebins = [20,100];
 % title_text = 'PLS contextual discrimination V1 RUN1 scatter (95th percentile)';
title_text = 'PLS contextual discrimination V1 RUN1 scatter (all bins 95th percentile)';
fig = figure('Name', title_text, 'Position', [200 100 640 580]); hold on;

for n = 1:2
    nexttile
    % Extract base data
    fpr = PLS_RUN_validations_all.FPR(1,:);
    % These should be (sessions x fpr_points)
    data_real = squeeze(PLS_RUN_validations_all.V1_TPR(n,:,:));
    data_shuf = squeeze(PLS_RUN_validations_all.V1_TPR_shuffled(n,:,:));
    
    % Calculate Means
    TPR_real_mean = mean(data_real, 1, 'omitnan');
    TPR_shuf_mean = mean(data_shuf, 1, 'omitnan');
    
    hold on
    % 1. Plot the diagonal chance line
    plot([0 1],[0 1],'--k', 'HandleVisibility', 'off')
    
    % 2. Plot individual session lines (Low Alpha)
    % Color definitions
    real_col = [231,41,138]/255;
    shuf_col = [0, 0, 0]/255;
    
    % Plot all individual Real lines
    p_indiv_real = plot(fpr, data_real', 'Color', [real_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % Plot all individual Shuffled lines
    p_indiv_shuf = plot(fpr, data_shuf', 'Color', [shuf_col, 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % 3. Plot Mean lines (Thick)
    PLOT(1) = plot(fpr, TPR_real_mean, 'Color', real_col, 'LineWidth', 2.5);
    PLOT(2) = plot(fpr, TPR_shuf_mean, 'Color', shuf_col, 'LineWidth', 2.5);
    
    % Formatting
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    title(sprintf('ROC curve V1 %i ms bins', timebins(n)));
    
    legend(PLOT(1:2), {'Real Mean', 'Shuffle Mean'}, 'Location', 'southeast', 'Box', 'off');
    
    xlim([0 1]); ylim([0 1]);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end



for n = 1:2
    nexttile
    % Extract data for the current condition
    real_data = PLS_RUN_validations_all.V1_AUC(n,:);
    shuffled_data = PLS_RUN_validations_all.V1_AUC_shuffled(n,:);
    
    % Combine for easier mean calculation
    % AUC = [real_data; shuffled_data];
    mean_real = mean(real_data, 'omitnan');
    mean_shuf = mean(shuffled_data, 'omitnan');
    bar_colors = [231,41,138; 0, 0, 0]/255;
    x_pos = [1 2];
    
    hold on

    % 1. Plot the Bars (Mean only, no error bars)
    plot([x_pos(1)-line_w, x_pos(1)+line_w], [mean_real, mean_real], 'Color', bar_colors(1,:), 'LineWidth', 4);
    plot([x_pos(2)-line_w, x_pos(2)+line_w], [mean_shuf, mean_shuf], 'Color', bar_colors(2,:), 'LineWidth', 4);

    % 2. Plot Raw Data Points and Connections with Jitter
    rng(1); % Seed for consistent jitter appearance
    jitter_strength = 0.2; % Adjusted for better visibility of distributions
    
    for s = 1:length(real_data)
        if ~isnan(real_data(s)) && ~isnan(shuffled_data(s))
            % Calculate jittered x-positions
            % This spreads the dots randomly around the x_pos center
            x1 = x_pos(1) + (rand - 0.5) * jitter_strength;
            x2 = x_pos(2) + (rand - 0.5) * jitter_strength;
            
            % Plot the connection line first (so it's behind the dots)
            plot([x1, x2], [real_data(s), shuffled_data(s)], 'Color', [0.8 0.8 0.8 0.4], 'LineWidth', 0.5);
            
            % Plot the individual dots
            scatter(x1, real_data(s), 50, bar_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.3);
            scatter(x2, shuffled_data(s), 50, bar_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.3);
        end
    end

    % Statistics and Formatting
    [p, h] = signrank(real_data, shuffled_data, 'tail', 'right');
    
    xlim([-0.5 3.5])
    xticks([1 2]);
    xticklabels({'Real', 'Shuffled'});
    ylabel('AUC');
    ylim([0 1.1]);
    yticks([0:0.25:1])
    title(['Condition ' num2str(n) ' (p=' num2str(p, '%.4e') ')']);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end
