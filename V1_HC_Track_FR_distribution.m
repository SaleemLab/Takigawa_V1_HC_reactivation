%% Code to 
analysis_folder = pwd;
load(fullfile(analysis_folder,'processed_data','spatial_map_all.mat'),'spatial_map_all');
% Initialize containers to gather data across all sessions
num_sessions =22;

%% Spatial tuning across all sessions
num_sessions = length(spatial_map_all.spatial_response);

% Initialize containers to gather data across all sessions
all_maps_trackL_odd  = [];
all_maps_trackL_even = [];
all_maps_trackR_odd  = [];
all_maps_trackR_even = [];
all_cell_labels      = {};

for nsession = 1:num_sessions
    % Fetch data structures for this session
    spatial_data = spatial_map_all.spatial_response{nsession};
    labels_data  = spatial_map_all.region{nsession};
    
    num_cells = size(spatial_data, 1);
    if num_cells == 0, continue; end
    
    % Store brain region labels for this session's cells
    all_cell_labels = [all_cell_labels; labels_data(:, 1)];
    
    % Temporary holders for this session's calculated averages
    sess_L_odd  = [];
    sess_L_even = [];
    sess_R_odd  = [];
    sess_R_even = [];
    
    for ncell = 1:num_cells
        % Track L (Column 1) & Track R (Column 2)
        map_L = spatial_data{ncell, 1}; 
        map_R = spatial_data{ncell, 2}; 
        
        % Track L Mean Maps (Odd vs. Even Laps)
        if ~isempty(map_L)
            sess_L_odd(ncell, :)  = mean(map_L(1:2:end, :), 1, 'omitnan');
            sess_L_even(ncell, :) = mean(map_L(2:2:end, :), 1, 'omitnan');
        else
            sess_L_odd(ncell, :)  = zeros(1, 140); 
            sess_L_even(ncell, :) = zeros(1, 140);
        end
        
        % Track R Mean Maps (Odd vs. Even Laps)
        if ~isempty(map_R)
            sess_R_odd(ncell, :)  = mean(map_R(1:2:end, :), 1, 'omitnan');
            sess_R_even(ncell, :) = mean(map_R(2:2:end, :), 1, 'omitnan');
        else
            sess_R_odd(ncell, :)  = zeros(1, 140);
            sess_R_even(ncell, :) = zeros(1, 140);
        end
    end
    
    % Append current session data to global matrix
    all_maps_trackL_odd  = [all_maps_trackL_odd;  sess_L_odd];
    all_maps_trackL_even = [all_maps_trackL_even; sess_L_even];
    all_maps_trackR_odd  = [all_maps_trackR_odd;  sess_R_odd];
    all_maps_trackR_even = [all_maps_trackR_even; sess_R_even];
end

% Replace remaining NaNs with 0
all_maps_trackL_odd(isnan(all_maps_trackL_odd))   = 0;
all_maps_trackL_even(isnan(all_maps_trackL_even)) = 0;
all_maps_trackR_odd(isnan(all_maps_trackR_odd))   = 0;
all_maps_trackR_even(isnan(all_maps_trackR_even)) = 0;

num_bins = size(all_maps_trackL_odd, 2); % Usually 140


% Loop Through the Target Regions and Generate Plots
% Define the list of regions you want to cycle through
target_regions = {'HPC', 'V1_L', 'V1_R'};

for r = 1:length(target_regions)
    target_region = target_regions{r};
    
    % Find cells matching the current region keyword
    target_indices = find(contains(all_cell_labels, target_region));
    
    if isempty(target_indices)
        warning('No cells found matching the region filter: %s. Skipping...', target_region);
        continue;
    end
    
    % Extract target sub-matrices for this region
    region_L_odd  = all_maps_trackL_odd(target_indices, :);
    region_L_even = all_maps_trackL_even(target_indices, :);
    region_R_odd  = all_maps_trackR_odd(target_indices, :);
    region_R_even = all_maps_trackR_even(target_indices, :);
    
    % --- Cross-Track Normalization & Sorting (ODD Laps) ---
    combined_odd = [region_L_odd, region_R_odd];
    norm_combined_odd = normalize(combined_odd, 2, 'range');
    
    norm_L_odd = norm_combined_odd(:, 1:num_bins);
    norm_R_odd = norm_combined_odd(:, num_bins+1:end);
    
    % Find independent peak locations for sorting reference
    [~, peak_idx_L] = max(norm_L_odd, [], 2);
    [~, peak_idx_R] = max(norm_R_odd, [], 2);
    
    [~, sort_order_L] = sort(peak_idx_L);
    [~, sort_order_R] = sort(peak_idx_R);
    
    % --- Cross-Track Normalization (EVEN Laps for Plotting) ---
    combined_even = [region_L_even, region_R_even];
    norm_combined_even = normalize(combined_even, 2, 'range');
    
    norm_L_even = norm_combined_even(:, 1:num_bins);
    norm_R_even = norm_combined_even(:, num_bins+1:end);
    
    % --- Generate Figures (2x2 Layout) ---
    fig = figure('Position', [50 + (r*40), 50 + (r*40), 700, 600]);
    sgtitle(sprintf('Population Maps for %s Cells (Cross-Track Normalized)', target_region), ...
        'FontSize', 15, 'FontWeight', 'bold');
    
    % 1. Track L (Even) sorted by Track L (Odd)
    subplot(2, 2, 1);
    imagesc(norm_L_even(sort_order_L, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track L (Even Laps) | Sorted by Track L (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
    
    % 2. Track L (Even) sorted by Track R (Odd)
    subplot(2, 2, 2);
    imagesc(norm_L_even(sort_order_R, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track L (Even Laps) | Sorted by Track R (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
    
    % 3. Track R (Even) sorted by Track L (Odd)
    subplot(2, 2, 3);
    imagesc(norm_R_even(sort_order_L, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track R (Even Laps) | Sorted by Track L (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
    
    % 4. Track R (Even) sorted by Track R (Odd)
    subplot(2, 2, 4);
    imagesc(norm_R_even(sort_order_R, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track R (Even Laps) | Sorted by Track R (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
end


%%%%
%%%%
%%%%
%% visualise Mean FR diff between Track L and Track R for V1 and HC
figure
subplot(2,2,1)
hold on
V1_L = [];V1_L_z=[];
V1_R = [];V1_R_z=[];
for nsession = 1:22
    V1_R_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'V1_R'),:);
    V1_L_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'V1_L'),:);
    scatter(V1_R_this(:,1),V1_R_this(:,2),'b','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency
    scatter(V1_L_this(:,1),V1_L_this(:,2), 'r','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency

    V1_L = [V1_L; V1_L_this];
    V1_R = [V1_R; V1_R_this];

    V1_R_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'V1_R'),:);
    V1_L_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'V1_L'),:);
    V1_L_z = [V1_L_z; V1_L_this];
    V1_R_z = [V1_R_z; V1_R_this];
end
xlim([0 50]);ylim([0 50])
plot([0 50],[0 50],'k--')
xlabel('Track L mean firing rate (Hz)')
ylabel('Track R mean firing rate (Hz)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(2,2,2)
histogram(V1_L(:,1)-V1_L(:,2),-20:1:20,'Normalization','probability')
hold on;
histogram(V1_R(:,1)-V1_R(:,2),-20:1:20,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (Hz)')
ylabel('Proportion of cells')
legend('V1 L','V1 R','box','off')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(2,2,3)
histogram(V1_L_z(:,1)-V1_L_z(:,2),-2:0.1:2,'Normalization','probability')
hold on;
histogram(V1_R_z(:,1)-V1_R_z(:,2),-2:0.1:2,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (z)')
ylabel('Proportion of cells')
sgtitle('V1')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

figure
subplot(2,2,1)
hold on
HC_L = [];HC_L_z=[];
HC_R = [];HC_R_z=[];
for nsession = 1:22
    HC_R_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'HPC_R'),:);
    HC_L_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'HPC_L'),:);
    scatter(HC_R_this(:,1),HC_R_this(:,2),'b','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency
    scatter(HC_L_this(:,1),HC_L_this(:,2), 'r','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency

    HC_L = [HC_L; HC_L_this];
    HC_R = [HC_R; HC_R_this];

    HC_R_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'HPC_R'),:);
    HC_L_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'HPC_L'),:);
    HC_L_z = [HC_L_z; HC_L_this];
    HC_R_z = [HC_R_z; HC_R_this];
end
xlim([0 50]);ylim([0 50])
plot([0 50],[0 50],'k--')
xlabel('Track L mean firing rate (Hz)')
ylabel('Track R mean firing rate (Hz)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(2,2,2)
histogram(HC_L(:,1)-HC_L(:,2),-20:1:20,'Normalization','probability')
hold on;
histogram(HC_R(:,1)-HC_R(:,2),-20:1:20,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (Hz)')
ylabel('Proportion of cells')
legend('HC L','HC R','box','off')

subplot(2,2,3)
histogram(HC_L_z(:,1)-HC_L_z(:,2),-2:0.1:2,'Normalization','probability')
hold on;
histogram(HC_R_z(:,1)-HC_R_z(:,2),-2:0.1:2,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (z)')
ylabel('Proportion of cells')
sgtitle('HC')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);