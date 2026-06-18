function ripple_modulation = ripple_modulation_analysis(spike_times,spike_id,window,psthBinSize,varargin)

% Default values
p = inputParser;
addParameter(p,'event_times',[],@isnumeric);
addParameter(p,'event_id',[],@isnumeric);
addParameter(p,'shuffle_option',1,@isnumeric);
addParameter(p,'unit_id',unique(spike_id),@isnumeric);
addParameter(p,'saving_PSTH',1,@isnumeric);
addParameter(p,'smooth_option',1,@isnumeric);


parse(p,varargin{:});
event_times = p.Results.event_times;
event_id =  p.Results.event_id;
unit_id = p.Results.unit_id;
shuffle_option = p.Results.shuffle_option;
saving_PSTH = p.Results.saving_PSTH;
smooth_option = p.Results.smooth_option;

event_conditions = unique(event_id(:))';
n_conditions = numel(event_conditions);
no_events = size(event_times,1);
time_edges = window(1):psthBinSize:window(2);
no_cluster = length(unit_id);

timevec = [min(spike_times) max(spike_times)];
timevec_edge = (timevec(1)-(psthBinSize)/2 : psthBinSize : timevec(end)+(psthBinSize)/2)';

% Gaussian smoothing window
gaussianWindow = gausswin(0.05*1/psthBinSize);
gaussianWindow = gaussianWindow / sum(gaussianWindow);

spike_times_events = spike_times;

for iCell = 1:no_cluster
    cluster_spike_id = spike_id == unit_id(iCell);
    y = histcounts(spike_times_events(cluster_spike_id), timevec_edge)';
%     y = conv(y, gaussianWindow, 'same');

    binnedArrays = cell(1, n_conditions);
    psth = cell(1, n_conditions);

    for ncond = 1:n_conditions
        cond_id = event_conditions(ncond);
        [~, bins, ~, ~, ~, binnedArray] = psthAndBA(spike_times_events(cluster_spike_id), ...
            event_times(event_id == cond_id), window, psthBinSize);

        if size(binnedArray,1) == 1
            smoothed = conv(binnedArray, gaussianWindow, 'same');
        else
            smoothed = zeros(size(binnedArray));
            for nevent = 1:size(binnedArray,1)
                smoothed(nevent,:) = conv(binnedArray(nevent,:), gaussianWindow, 'same');
            end
        end

        if smooth_option==1
            binnedArrays{ncond} = smoothed;
            psth{ncond} = mean(smoothed, 'omitnan');
        else
            binnedArrays{ncond} = binnedArray;
            psth{ncond} = mean(binnedArray, 'omitnan');
        end
    end

    if shuffle_option == 1
        PSTH_shuffled = cell(1, n_conditions);
        for ncond = 1:n_conditions
            this_binned = binnedArrays{ncond};
            n_trials = size(this_binned, 1);
            n_timepoints = size(this_binned, 2);
            % PSTH_shuffled{ncond} = zeros(1000, n_timepoints);
            temp = zeros(1000, n_timepoints);

            parfor nshuffle = 1:1000
                s = RandStream('mrg32k3a','Seed',nshuffle+iCell*1000+ncond*10000);
                shift_values = round(n_timepoints * rand(s, 1, n_trials))';
                shifted = arrayfun(@(nevent) circshift(this_binned(nevent,:), shift_values(nevent), 2), ...
                    1:n_trials, 'UniformOutput', false);
                temp(nshuffle,:) = mean(cell2mat(shifted'), 'omitnan');
            end
            PSTH_shuffled{ncond} = temp;
        end
    end

    for ncond = 1:n_conditions
        ripple_modulation(ncond).bins = bins;

        if saving_PSTH==1
            ripple_modulation(ncond).PSTH(iCell,:,:) = single(binnedArrays{ncond}) / psthBinSize;  % normalize by condition 1
            ripple_modulation(ncond).PSTH_zscored(iCell,:,:) = single((binnedArrays{ncond} - mean(y)) ./ std(y));
            ripple_modulation(ncond).PSTH_zscored(iCell,:,:) = single((binnedArrays{ncond} - mean(y)) ./ std(y));
        end

        if shuffle_option == 1
            twin = bins > 0 & bins < 0.2;
            real_psth = psth{ncond};
            base_psth = mean(PSTH_shuffled{ncond}(:, twin), 1);
            MSD_real = mean((real_psth(twin) - base_psth).^2);

            FR_diff = zeros(1, 1000);
            MSD_shuffled = zeros(1, 1000);
            for nshuffle = 1:1000
                shuffled = PSTH_shuffled{ncond}(nshuffle,twin);
                FR_diff(nshuffle) = (max(shuffled) - min(shuffled)) / mean(shuffled);
                MSD_shuffled(nshuffle) = mean((shuffled - base_psth).^2);
            end

            FR_mod = (max(real_psth(twin)) - min(real_psth(twin))) / mean(real_psth(twin));
            ripple_modulation(ncond).ripple_modulation(iCell) = FR_mod;
            ripple_modulation(ncond).ripple_modulation_peak_percentile(iCell) = ...
                sum(abs(FR_mod) > abs(FR_diff)) / numel(FR_diff);
            ripple_modulation(ncond).ripple_modulation_percentile(iCell) = ...
                mean(MSD_shuffled < MSD_real);

            % PRE window
            twin_pre = bins > -0.2 & bins < 0;
            real_psth_pre = psth{ncond};
            base_psth_pre = mean(PSTH_shuffled{ncond}(:, twin_pre), 1);
            MSD_real_pre = mean((real_psth_pre(twin_pre) - base_psth_pre).^2);

            FR_diff_pre = zeros(1, 1000);
            MSD_shuffled_pre = zeros(1, 1000);
            for nshuffle = 1:1000
                shuffled = PSTH_shuffled{ncond}(nshuffle, twin_pre);
                FR_diff_pre(nshuffle) = (max(shuffled) - min(shuffled)) / mean(shuffled);
                MSD_shuffled_pre(nshuffle) = mean((shuffled - base_psth_pre).^2);
            end

            FR_mod_pre = (max(real_psth_pre(twin_pre)) - min(real_psth_pre(twin_pre))) / mean(real_psth_pre(twin_pre));
            ripple_modulation(ncond).peak_percentile_PRE(iCell) = ...
                sum(abs(FR_mod_pre) > abs(FR_diff_pre)) / numel(FR_diff_pre);
            ripple_modulation(ncond).modulation_percentile_PRE(iCell) = ...
                mean(MSD_shuffled_pre < MSD_real_pre);
        end
    end
end
end
