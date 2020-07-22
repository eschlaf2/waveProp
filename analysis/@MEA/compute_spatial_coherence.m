function coh_sp = compute_spatial_coherence(mea)
% coh_sp = self.compute_spatial_coherence(mea)
% returns spatial coherence at each time point
if strcmpi(mea.patient, 'c7')
    data_i = WaveFlow.interpolate_nans( ...
        mea.Data(mea.Time > 0 & mea.Time < 30, :), mea.Position);
else
    data_i = WaveFlow.interpolate_nans( ...
        mea.Data(mea.Active, :), mea.Position);
end
params.tapers = [3 5]; 
cx = arrayfun(@(tt) coherencyc(data_i(:, 5*ones(1, 9), tt), data_i(:, [1:4 6:10], tt), params), 1:length(data_i), 'uni', 0);

cy = arrayfun(@(tt) coherencyc(data_i(5*ones(1, 9), :, tt)', data_i([1:4 6:10], :, tt)', params), 1:length(data_i), 'uni', 0);

cx_mn = cellfun(@(x) mean(x, [1 2]), cx);
cy_mn = cellfun(@(x) mean(x, [1 2]), cy);

coh_sp = mean([cx_mn.^2; cy_mn.^2]);

% save(['coh_sp/' strrep(names{ii}, ' ', '_')], 'coh_sp', 'cx', 'cy'); 

end
