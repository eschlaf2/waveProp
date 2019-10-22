function out = keep_random_channels(mea, n)
% Keep n random channels from mea

out = false;  % Initialize action indicator
channels = 1:size(mea.Data,2);  % Get all recorded channels
channels(mea.BadChannels) = [];  % Exclude bad channels
nC = numel(channels);  % convenience

if n >= numel(channels), disp('No channels removed.'), return, end

out = true;  % Indicate changes being made
bad_ch = channels(randperm(nC, nC - n));  % Randomly choose additional bad channels
mea.BadChannels = unique([bad_ch mea.BadChannels]);  % Add them to mea
mea.excluded_channels = bad_ch;  % Store which channels were removed
mea = rmfield(mea, {'event_inds', 'firingRate'});  % These should be recomputed later if needed

assignin('base', 'mea', mea)

end