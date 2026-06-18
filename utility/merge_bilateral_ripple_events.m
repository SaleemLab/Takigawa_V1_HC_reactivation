function [event_ids_first,event_ids_second] = merge_bilateral_ripple_events(hemisphere_id,event_times,max_diff);
%%%% Merge bilateral ripple events

% Separate indices
idx_hemi1 = find(hemisphere_id == 1);
idx_hemi2 = find(hemisphere_id == 2);

% Sort event times for both hemispheres
[pt1, sort1] = sort(event_times(idx_hemi1));
[pt2, sort2] = sort(event_times(idx_hemi2));
idx1_sorted = idx_hemi1(sort1);
idx2_sorted = idx_hemi2(sort2);

% Initialize matched flags
used1 = false(size(idx1_sorted));
used2 = false(size(idx2_sorted));

% Store synchronised pairs
pair_first = [];
pair_second = [];

i = 1; j = 1;
while i <= length(pt1) && j <= length(pt2)
    dt = pt1(i) - pt2(j);
    if abs(dt) <= max_diff
        if ~used1(i) && ~used2(j)
            if dt <= 0
                pair_first(end+1) = idx1_sorted(i);
                pair_second(end+1) = idx2_sorted(j);
            else
                pair_first(end+1) = idx2_sorted(j);
                pair_second(end+1) = idx1_sorted(i);
            end
            used1(i) = true;
            used2(j) = true;
            i = i + 1;
            j = j + 1;
        else
            if used1(i), i = i + 1; end
            if used2(j), j = j + 1; end
        end
    elseif dt < -max_diff
        i = i + 1;
    else
        j = j + 1;
    end
end

% Combine all synchronised ripple indices
paired_idx = [pair_first, pair_second];

% Unilateral events: not in any pair
all_idx = 1:numel(event_times);
unilateral_idx = setdiff(all_idx, paired_idx);

% Define the two output event ID sets
[event_ids_first sort_idx]  = sort([unilateral_idx, pair_first]);
event_ids_second = [unilateral_idx, pair_second];
event_ids_second = event_ids_second(sort_idx);