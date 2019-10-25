function metrics = rename_metrics(metrics)

map = {'maxdescent', 'M'; ...
	'events', 'E'; ...
	'delays_T01_fband1_50', 'D'; ...
	'delays_T0p2_fband0_50', 'Ds'; ...
	'delays_T10_fband1_13', 'D10'};

if ischar(metrics)
	metrics = map{strcmp(metrics, map(:, 1)), 2};
else
	for ii = 1:length(metrics)
		metrics{ii} = map{strcmp(metrics{ii}, map(:, 1)), 2};
	end
end

end