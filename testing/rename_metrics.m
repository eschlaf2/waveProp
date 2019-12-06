function metrics = rename_metrics(metrics)

map = {'maxdescent', 'M'; ...
	'events', 'E'; ...
	'delays_T01_fband1_50', 'D1w'; ...
	'delays_T0p2_fband0_50', 'Ds'; ...
	'delays_T10_fband1_13', 'D10'; ...
	'groundtruth', 'G'; ...
	'delays_T01_fband1_13', 'D1'};

if ischar(metrics)
	if ismember(metrics, map(:, 1))
		metrics = map{strcmp(metrics, map(:, 1)), 2};
	end
else
	for ii = 1:numel(metrics)
		if ismember(metrics{ii}, map)
			metrics{ii} = map{strcmp(metrics{ii}, map(:, 1)), 2};
		end
	end
end

end