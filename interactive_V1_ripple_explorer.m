%% interactive_V1_ripple_explorer.m
function interactive_V1_ripple_explorer()
    %% Dynamically Identify the Current Repository Folder
    root_folder = pwd;
    if isempty(root_folder)
        root_folder = pwd; 
    end
    
    % Redirect paths to point to the processed_data folder inside the repository
    data_folder = fullfile(root_folder, 'processed_data');
    lme_folder  = fullfile(data_folder, 'processed_data');
    out_fldr    = data_folder; % Saves figures cleanly directly inside the processed_data folder
    
    %% Load and combine spatial and ripple PSTH Data
    files = {'V1_ripple_psth.mat', 'V1_low_ripple_psth.mat', 'V1_high_ripple_psth.mat','V1_spatial_psth.mat'};

    for f = 1:length(files)
        path = fullfile(data_folder, files{f});
        if ~exist(path, 'file'), error('Missing: %s', path); end
        tmp = load(path);
        vnames = fieldnames(tmp);
        curr_struct = tmp.(vnames{1}); % Extract main struct variable

        if f == 1
            V1_cells_data = curr_struct;
        else
            % Dynamically merge all fields from the other 3 files
            fields = fieldnames(curr_struct);
            for i = 1:length(V1_cells_data)
                for g = 1:length(fields)
                    V1_cells_data(i).(fields{g}) = curr_struct(i).(fields{g});
                end
            end
        end
    end
    n_cells = length(V1_cells_data);
    x_bins = V1_cells_data.x_bins;

        
     %% Load LME Regression Results
    lme_path = fullfile(lme_folder, 'V1_ripple_modulation_all.mat');
    if exist(lme_path, 'file')
        lme_data = load(lme_path, 'ripple_modulation_lme', 'ripple_modulation_lme_low', 'ripple_modulation_lme_high');
    else
        warning('LME result file not found. Save function will lack LME overlay.');
        lme_data = [];
    end
    
    % Precompute X (Track Diff) and Y (Ripple Diff) for UI
    windows = {'ALL', 'PRE', 'POST'};
    win_ranges = {[-0.2, 0.2], [-0.2, 0], [0, 0.2]};
    conditions = {'All ripples', 'Low ripples', 'High ripples'};
    
    % Load exact CSV data from the modulation subfolder
    f_all = fullfile(lme_folder, 'ripple_track_difference_ALL.csv');
    f_low = fullfile(lme_folder, 'ripple_track_difference_LOW.csv');
    f_high= fullfile(lme_folder, 'ripple_track_difference_HIGH.csv');
    
    session_ids = [V1_cells_data.session_id]';
    cell_ids = [V1_cells_data.cell_id]';
    x_track = zeros(n_cells, 1);
    z_track_diff = [V1_cells_data.z_track_diff];
    
    y_data = struct();
    y_data.All = zeros(n_cells, 3);
    y_data.Low = zeros(n_cells, 3);
    y_data.High= zeros(n_cells, 3);
    
    if exist(f_all, 'file') && exist(f_low, 'file') && exist(f_high, 'file')
        csv_all = readtable(f_all);
        csv_low = readtable(f_low);
        csv_high= readtable(f_high);
        
        x_track = csv_all.track_difference;
        y_data.All = [csv_all.ripple_difference_ALL, csv_all.ripple_difference_PRE, csv_all.ripple_difference_POST];
        y_data.Low = [csv_low.ripple_difference_ALL, csv_low.ripple_difference_PRE, csv_low.ripple_difference_POST];
        y_data.High= [csv_high.ripple_difference_ALL, csv_high.ripple_difference_PRE, csv_high.ripple_difference_POST];
    else
        % Fallback manual calculation if CSVs missing
        for i = 1:n_cells
            d = V1_cells_data(i);
            x_track(i) = d.z_track_diff;
            for w = 1:3
                idx = x_bins >= win_ranges{w}(1) & x_bins <= win_ranges{w}(2);
                y_data.All(i, w) = mean(d.z_all_T1_m(idx)) - mean(d.z_all_T2_m(idx));
                y_data.High(i, w) = mean(d.z_high_T1_m(idx)) - mean(d.z_high_T2_m(idx));
                y_data.Low(i, w) = mean(d.z_low_T1_m(idx)) - mean(d.z_low_T2_m(idx));
            end
        end
    end
    
    %% Setup UI
    fig = uifigure('Name', 'Interactive V1 Ripple Explorer', 'Position', [100 100 1400 700]);
    gl = uigridlayout(fig, [1, 2]);
    gl.ColumnWidth = {'1x', '1.5x'};
    
    % Left Panel: Controls + Scatter
    left_panel = uipanel(gl);
    left_gl = uigridlayout(left_panel, [2, 1]);
    left_gl.RowHeight = {110, '1x'};
    
    % Controls
    control_pnl = uipanel(left_gl, 'BorderType', 'none');
    
    uilabel(control_pnl, 'Text', 'Condition:', 'Position', [10 80 60 22]);
    cond_dd = uidropdown(control_pnl, 'Items', conditions, 'Position', [75 80 120 22], 'ValueChangedFcn', @(s, e) update_scatter());
    
    uilabel(control_pnl, 'Text', 'Window:', 'Position', [210 80 50 22]);
    win_dd = uidropdown(control_pnl, 'Items', windows, 'Position', [265 80 80 22], 'ValueChangedFcn', @(s, e) update_scatter());
    
    uilabel(control_pnl, 'Text', 'Value type:', 'Position', [10 50 65 22]);
    val_dd = uidropdown(control_pnl, 'Items', {'Raw (Hz)', 'Z-Scored'}, 'Position', [75 50 120 22], 'ValueChangedFcn', @(s, e) redraw_psths());
    
    % Manual Entry
    uilabel(control_pnl, 'Text', 'Session:', 'Position', [10 20 50 22]);
    sess_edit = uieditfield(control_pnl, 'numeric', 'Position', [60 20 40 22]);
    uilabel(control_pnl, 'Text', 'Cell:', 'Position', [110 20 40 22]);
    cell_edit = uieditfield(control_pnl, 'numeric', 'Position', [140 20 40 22]);
    uibutton(control_pnl, 'Text', 'Add Cell', 'Position', [190 20 70 22], 'ButtonPushedFcn', @(s, e) add_manual_cell());
    uibutton(control_pnl, 'Text', 'Clear Cells', 'Position', [270 20 80 22], 'ButtonPushedFcn', @(s, e) clear_cells());
    
    % Save Button
    uibutton(control_pnl, 'Text', 'Save Figures', 'Position', [360 80 100 40], ...
             'BackgroundColor', [0.8 0.9 1], 'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) save_figures());
             
    ax_scatter = uiaxes(left_gl);
    hold(ax_scatter, 'on');
    
    % Right Panel: Redesigned Row-Based Layout (Row 1 = Cell 1, Row 2 = Cell 2)
    right_panel = uipanel(gl, 'Title', 'Select up to 2 dots to view PSTHs');
    right_gl = uigridlayout(right_panel, [2, 2]);
    
    % Row 1: Cell 1 Panels
    ax_spat_c1 = uiaxes(right_gl); title(ax_spat_c1, 'Cell 1: Spatial PSTH'); hold(ax_spat_c1, 'on');
    ax_rip_c1  = uiaxes(right_gl); title(ax_rip_c1, 'Cell 1: Ripple PSTH'); hold(ax_rip_c1, 'on');
    
    % Row 2: Cell 2 Panels
    ax_spat_c2 = uiaxes(right_gl); title(ax_spat_c2, 'Cell 2: Spatial PSTH'); hold(ax_spat_c2, 'on');
    ax_rip_c2  = uiaxes(right_gl); title(ax_rip_c2, 'Cell 2: Ripple PSTH'); hold(ax_rip_c2, 'on');
    
    % State Variables
    selected_cells = []; 
    colors = [175, 141, 195; 127, 191, 123] / 256; 
    hi_pts = gobjects(1, 2); 
    
    cur_xv = []; cur_yv = []; cur_orig_idx = [];
    
    %% Nested Functions
    function update_scatter()
        c_idx = find(strcmp(conditions, cond_dd.Value));
        w_idx = find(strcmp(windows, win_dd.Value));
        
        cla(ax_scatter);
        
        x = x_track;
        if c_idx == 1; y = y_data.All(:, w_idx);
        elseif c_idx == 2; y = y_data.Low(:, w_idx);
        else; y = y_data.High(:, w_idx); end
        
        valid = ~isnan(x) & ~isnan(y);
        cur_xv = x(valid); cur_yv = y(valid);
        cur_orig_idx = find(valid);
        
        if length(cur_xv) > 1
            fitted = polyfit(cur_xv, cur_yv, 1);
            y_fit = polyval(fitted, [min(cur_xv) max(cur_xv)]);
            
            p_text = '';
            if ~isempty(lme_data)
                if c_idx == 1; lme_struc = lme_data.ripple_modulation_lme;
                elseif c_idx == 2; lme_struc = lme_data.ripple_modulation_lme_low;
                else; lme_struc = lme_data.ripple_modulation_lme_high; end
                beta = lme_struc(w_idx).model.fixedEffects;
                if sign(beta(2)) ~= sign(diff(y_fit))
                    y_fit(2) = y_fit(1) - diff(y_fit);
                end
                f_pval = lme_struc(w_idx).p;
                f_stat = lme_struc(w_idx).t_stat^2;
                p_text = sprintf('p = %.3e\nF = %.3f', f_pval, f_stat);
            end
            
            plot(ax_scatter, [min(cur_xv) max(cur_xv)], y_fit, 'r-', 'LineWidth', 2, 'HitTest', 'off');
            
            if ~isempty(p_text)
                text(ax_scatter, -2, 0.1, p_text, 'FontSize', 12, 'HitTest', 'off');
            end
        end
        xline(ax_scatter, 0, '--', 'HitTest', 'off'); yline(ax_scatter, 0, '--', 'HitTest', 'off');
        
        h_scatter = scatter(ax_scatter, cur_xv, cur_yv, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2);
        
        ylim(ax_scatter, [-0.2 0.2]); xlim(ax_scatter, [-2.5 2.5]);
        xlabel(ax_scatter, 'Track difference (z)'); ylabel(ax_scatter, 'Ripple difference (z)');
        title(ax_scatter, sprintf('%s (%s)', cond_dd.Value, win_dd.Value));
        ax_scatter.TickDir = 'out'; ax_scatter.Box = 'off';
        
        h_scatter.ButtonDownFcn = @(src, event) scatter_clicked(event);
        
        redraw_highlights();
        redraw_psths();
    end
    function scatter_clicked(event)
        click_x = event.IntersectionPoint(1);
        click_y = event.IntersectionPoint(2);

        % Calculate raw squared Euclidean distance without warping the axes
        distances = (cur_xv - click_x).^2 + (cur_yv - click_y).^2;
        [~, min_idx] = min(distances);

        add_to_selection(cur_orig_idx(min_idx));
    end
    function add_manual_cell()
        s_val = sess_edit.Value;
        c_val = cell_edit.Value;
        idx = find(session_ids == s_val & cell_ids == c_val, 1);
        if ~isempty(idx)
            add_to_selection(idx);
        else
            uialert(fig, 'Cell not found!', 'Error');
        end
    end
    function clear_cells()
        selected_cells = [];
        redraw_highlights();
        redraw_psths();
        right_panel.Title = 'Select up to 2 dots to view PSTHs';
    end
    function add_to_selection(idx)
        if length(selected_cells) >= 2
            selected_cells(2) = idx; 
        else
            selected_cells = [selected_cells, idx];
        end
        redraw_highlights();
        redraw_psths();
    end
    function redraw_highlights()
        delete(hi_pts(isgraphics(hi_pts)));

        for k = 1:length(selected_cells)
            cell_idx = selected_cells(k);

            % Get current X and Y directly from our data arrays using the cell index
            c_idx = find(strcmp(conditions, cond_dd.Value));
            w_idx = find(strcmp(windows, win_dd.Value));

            x_val = x_track(cell_idx);
            if c_idx == 1; y_val = y_data.All(cell_idx, w_idx);
            elseif c_idx == 2; y_val = y_data.Low(cell_idx, w_idx);
            else; y_val = y_data.High(cell_idx, w_idx); end

            % Plot directly at the correct coordinates
            hi_pts(k) = scatter(ax_scatter, x_val, y_val, 100, colors(k,:), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HitTest', 'off');
        end

        % Update title panel string logic
        if isempty(selected_cells)
            right_panel.Title = 'Select up to 2 dots to view PSTHs';
        else
            str = '';
            for k = 1:length(selected_cells)
                [~, index] = min(abs(x_track(selected_cells(k)) - z_track_diff));
                d = V1_cells_data(index);

                if c_idx == 1; y_val = y_data.All(selected_cells(k), w_idx);
                elseif c_idx == 2; y_val = y_data.Low(selected_cells(k), w_idx);
                else; y_val = y_data.High(selected_cells(k), w_idx); end

                str = [str, sprintf('Cell %d (S%d) X:%.2f Y:%.2f | ', d.cell_id, d.session_id, d.z_track_diff, y_val)];
            end
            right_panel.Title = str(1:end-3);
        end
    end


    function redraw_psths()
        % Clear all four targeted UI axes panels
        cla(ax_spat_c1); cla(ax_rip_c1);
        cla(ax_spat_c2); cla(ax_rip_c2);

        use_z = strcmp(val_dd.Value, 'Z-Scored');
        pos_bins = 1:1:140;

        % Iterate and direct data streams to dedicated rows instead of overlaying
        for k = 1:length(selected_cells)
            [~, index] = min(abs(x_track(selected_cells(k)) - z_track_diff));
            d = V1_cells_data(index);

            % d = V1_cells_data(selected_cells(k));

            % Set color profile targets based on row assignments
            col1 = [0 0 1];  if k==2; col1 = [0.5 0.5 1]; end
            col2 = [1 0 0];  if k==2; col2 = [1 0.5 0.5]; end

            % Pick row axes targets
            if k == 1
                ax_s = ax_spat_c1; ax_r = ax_rip_c1;
                title(ax_s, sprintf('Cell %d (S%d): Spatial', d.cell_id, d.session_id));
                title(ax_r, sprintf('Cell %d (S%d): Ripple', d.cell_id, d.session_id));
            else
                ax_s = ax_spat_c2; ax_r = ax_rip_c2;
                title(ax_s, sprintf('Cell %d (S%d): Spatial', d.cell_id, d.session_id));
                title(ax_r, sprintf('Cell %d (S%d): Ripple', d.cell_id, d.session_id));
            end

            spatial_all = [d.spatial_1_m, d.spatial_2_m];
            sd_sp = std(spatial_all); if sd_sp==0; sd_sp=1; end
            mu_sp = mean(spatial_all);

            if use_z
                plot_psth(ax_s, pos_bins, (d.spatial_1_m-mu_sp)/sd_sp, d.spatial_1_se/sd_sp, col1);
                plot_psth(ax_s, pos_bins, (d.spatial_2_m-mu_sp)/sd_sp, d.spatial_2_se/sd_sp, col2);

                plot_psth(ax_r, x_bins, d.z_all_T1_m, d.z_all_T1_se, col1);
                plot_psth(ax_r, x_bins, d.z_all_T2_m, d.z_all_T2_se, col2);
            else
                plot_psth(ax_s, pos_bins, d.spatial_1_m, d.spatial_1_se, col1);
                plot_psth(ax_s, pos_bins, d.spatial_2_m, d.spatial_2_se, col2);

                plot_psth(ax_r, x_bins, d.raw_all_T1_m, d.raw_all_T1_se, col1);
                plot_psth(ax_r, x_bins, d.raw_all_T2_m, d.raw_all_T2_se, col2);
            end
        end

        % Enforce baseline limits and indicator decorations across active views
        xlim(ax_spat_c1, [0 140]); xlim(ax_spat_c2, [0 140]);
        xlim(ax_rip_c1, [-0.5 0.5]); xline(ax_rip_c1, 0, '--');
        xlim(ax_rip_c2, [-0.5 0.5]); xline(ax_rip_c2, 0, '--');
    end

    function plot_psth(ax, x_vals, m, se, color)
        hold(ax, 'on');
        if isempty(m) || all(isnan(m(:))), return; end
        m = m(:)'; se = se(:)'; x_vals = x_vals(:)';
        valid = ~isnan(m) & ~isnan(se);
        if ~any(valid), return; end
        m = m(valid); se = se(valid); x_vals = x_vals(valid);
        fill(ax, [x_vals fliplr(x_vals)], [m-se fliplr(m+se)], color, 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HitTest', 'off');
        plot(ax, x_vals, m, 'Color', color, 'LineWidth', 1.2, 'HitTest', 'off');
    end
    %% Save Figures Function
    function save_figures()
        if isempty(selected_cells)
            uialert(fig, 'Please select at least 1 cell to save plots.', 'Info'); return;
        end
        
        if ~exist(out_fldr, 'dir'), mkdir(out_fldr); end
        
        % 1. Regression Figure (3x3 Grid)
        fReg = figure('Name', 'regression_3x3_scatter', 'Position', [60 60 1200 900]);
        csv_files = {'ripple_track_difference_ALL.csv', 'ripple_track_difference_LOW.csv', 'ripple_track_difference_HIGH.csv'};
        c_names  = {'All', 'Low', 'High'};
        w_names  = {'All', 'PRE', 'POST'};
        plot_idx = 1;
        
        for c = 1:3
            csv_path = fullfile(lme_folder, csv_files{c});
            if exist(csv_path, 'file') && ~isempty(lme_data)
                tbl = readtable(csv_path);
                x_full = tbl.track_difference;
                
                if c == 1; lme_struc = lme_data.ripple_modulation_lme;
                elseif c == 2; lme_struc = lme_data.ripple_modulation_lme_low;
                else; lme_struc = lme_data.ripple_modulation_lme_high; end
                
                for w = 1:3
                    subplot(3, 3, plot_idx);
                    
                    y = table2array(tbl(:, 1 + w));
                    x = x_full(~isnan(y)); y = y(~isnan(y));
                    
                    f_pval = lme_struc(w).p;
                    f_stat = lme_struc(w).t_stat^2;
                    beta = lme_struc(w).model.fixedEffects;
                    
                    fitted = polyfit(x, y, 1);
                    y_fit = polyval(fitted, [min(x) max(x)]);
                    
                    if sign(beta(2)) ~= sign(diff(y_fit))
                        y_fit(2) = y_fit(1) - diff(y_fit);
                    end
                    
                    scatter(x, y, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1); hold on;
                    plot([min(x) max(x)], y_fit, 'r-', 'LineWidth', 2);
                    xline(0, '--'); yline(0, '--');
                    text(-2, 0.1, sprintf('p = %.3e\nF = %.3f', f_pval, f_stat), 'FontSize', 10);
                    
                    % Overlay highlights
                    for k = 1:length(selected_cells)
                        
                        d = V1_cells_data(find(x_track(selected_cells(k)==z_track_diff)));
                        % d = V1_cells_data(selected_cells(k));
                        scatter(d.z_track_diff, y_data.(strrep(c_names{c}, ' ripples', ''))(selected_cells(k), w), 60, colors(k,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                    end
                    
                    ylim([-0.2 0.2]); xlim([-2.5 2.5]);
                    xlabel('Track difference (z)'); ylabel('Ripple difference (z)');
                    title(sprintf('%s %s', w_names{w}, lower(c_names{c})));
                    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
                    
                    plot_idx = plot_idx + 1;
                end
            else
                plot_idx = plot_idx + 3; 
            end
        end
        
        if ~exist('save_all_figures', 'file')
            saveas(fReg, fullfile(out_fldr, 'regression_3x3_scatter.pdf'));
        end
        
        % 2. PSTH Figures
        for k = 1:length(selected_cells)
             % d = V1_cells_data(find(x_track(selected_cells(k))==z_track_diff));
            % d = V1_cells_data(selected_cells(k));
            [~, index] = min(abs(x_track(cell_indx) - z_track_diff));
            d = V1_cells_data(index);
            fPsth = figure('Name', sprintf('PSTH_Sess_%d_Cell_%d', d.session_id, d.cell_id), 'Position', [100 100 1200 600]);
            
            % Raw
            subplot(2, 4, 1);
            hold on;
            plot_psth(gca, 1:1:140, d.spatial_1_m, d.spatial_1_se, [0 0 1]); hold on; plot_psth(gca, 1:1:140, d.spatial_2_m, d.spatial_2_se, [1 0 0]);
            xlim([0 140]); title(sprintf('Cell %d (S %d)\nSpatial (Hz)', d.cell_id, d.session_id));
            subplot(2, 4, 2);
            plot_psth(gca, x_bins, d.raw_all_T1_m, d.raw_all_T1_se, [0 0 1]); hold on; plot_psth(gca, x_bins, d.raw_all_T2_m, d.raw_all_T2_se, [1 0 0]);
            xlim([-0.5 0.5]); xline(0,'--'); title('All (Hz)');
            subplot(2, 4, 3);
            plot_psth(gca, x_bins, d.raw_high_T1_m, d.raw_high_T1_se, [0 0 1]); hold on; plot_psth(gca, x_bins, d.raw_high_T2_m, d.raw_high_T2_se, [1 0 0]);
            xlim([-0.5 0.5]); xline(0,'--'); title('High (Hz)');
            subplot(2, 4, 4);
            plot_psth(gca, x_bins, d.raw_low_T1_m, d.raw_low_T1_se, [0 0 1]); hold on; plot_psth(gca, x_bins, d.raw_low_T2_m, d.raw_low_T2_se, [1 0 0]);
            xlim([-0.5 0.5]); xline(0,'--'); title('Low (Hz)');
            
            % Spatial Z-score configuration
            spatial_all = [d.spatial_1_m, d.spatial_2_m];
            sd_sp = std(spatial_all); if sd_sp==0; sd_sp=1; end
            mu_sp = mean(spatial_all);
            
            % ZScored
            subplot(2, 4, 5); 
            plot_psth(gca, 1:1:140, (d.spatial_1_m-mu_sp)/sd_sp, d.spatial_1_se/sd_sp, [0 0 1]); hold on; 
            plot_psth(gca, 1:1:140, (d.spatial_2_m-mu_sp)/sd_sp, d.spatial_2_se/sd_sp, [1 0 0]);
            xlim([0 140]); title('Spatial (Z)'); 
            
            subplot(2, 4, 6);
            plot_psth(gca, x_bins, d.z_all_T1_m, d.z_all_T1_se, [0 0 1]); hold on; plot_psth(gca, x_bins, d.z_all_T2_m, d.z_all_T2_se, [1 0 0]);
            xlim([-0.5 0.5]); xline(0,'--'); title('All (Z)');
            subplot(2, 4, 7);
            plot_psth(gca, x_bins, d.z_high_T1_m, d.z_high_T1_se, [0 0 1]); hold on; plot_psth(gca, x_bins, d.z_high_T2_m, d.z_high_T2_se, [1 0 0]);
            xlim([-0.5 0.5]); xline(0,'--'); title('High (Z)');
            subplot(2, 4, 8);
            plot_psth(gca, x_bins, d.z_low_T1_m, d.z_low_T1_se, [0 0 1]); hold on; plot_psth(gca, x_bins, d.z_low_T2_m, d.z_low_T2_se, [1 0 0]);
            xlim([-0.5 0.5]); xline(0,'--'); title('Low (Z)');
            
            for ax_idx = 1:8
                ax = subplot(2,4,ax_idx); ax.TickDir = 'out'; ax.Box = 'off';
            end
            
            if ~exist('save_all_figures', 'file')
                saveas(fPsth, fullfile(out_fldr, sprintf('PSTH_S%d_C%d.pdf', d.session_id, d.cell_id)));
            end
        end
        
        if exist('save_all_figures', 'file')
            save_all_figures(out_fldr, []);
        end
        
        uialert(fig, ['Figures saved successfully inside your repository under: ' out_fldr], 'Saved');
    end
    % Initialize
    update_scatter();
end