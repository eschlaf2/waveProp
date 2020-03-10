function res = compile_wave_prop(varargin)
% Create <res>. Compiles computed wave properties from the files produced
% after analyzing wave directions.
% Inputs are in the form of name-value pairs
%     res = compile_wave_prop('files', files);  % files is a struct
%     res = compile_wave_prop('pat', pats);  % pats is a char array or
%		cell of strings
%     res = compile_wave_prop('pat', '*', 'seizure', '*');
%       % this will pull pat/seizure pairs from seizures2.txt

args = parse_inputs(varargin);

files = args.files;
metrics = args.metrics;
sig = args.sig;


%% Make res

nF = numel(files);
res(nF) = struct(...
	'name', [], ...
	'data', [], ...
	'Z', [], ...
	'time', [], ...
	'p', [], ...
	'Vx', [], ...
	'Vy', []);

for ii = 1:nF
	finfo = strsplit(files(ii).name, {'_', '.'});  % Get the file name
	res(ii).name = sprintf('%s %s', finfo{1}, finfo{2}(8:end));  % get the patient and seizure number
	res(ii).data = load(fullfile(files(ii).folder, files(ii).name), metrics{:});  % store the data used to fit each wave
	res(ii).time = get_unique_wavetimes(res(ii).data);  % Get the times of waves detected for each method (these should be identical, but prior analyses ran delays only every 1 second)

	time = res(ii).time;
	fields = fieldnames(res(ii).data);
	[res(ii).Z, res(ii).Vx, res(ii).Vy] = ...
		deal(zeros(length(time), length(fields)));
	for jj = 1:numel(fields)  % for each metric
		for F = 'Zpxy'  % and each property
			
			[data, f] = ...  % match results to unique times
				interpolate_data(F, res(ii).data.(fields{jj}), time);
			res(ii).(f)(:, jj) = data;  
			mask = ...  % Remove values where fit is not significant
				res(ii).data.(fields{jj}).p >= sig;
			res(ii).(f)(mask, jj) = nan;
			mask = ...  % Remove values where slope is zero in both directions
				all(abs(res(ii).data.(fields{jj}).beta(1:2, :)) < eps);
			res(ii).(f)(mask, jj) = nan;
			
		end
	end
end
        
end

%% HELPERS

function [dataI, f] = interpolate_data(f, data, time)
% Match wave properties to nearest wave time and rename f. Originally this
% was to control for delay metric which only computed wave properties once
% per second. Now, wave properties are computed at each discharge so this
% should be identical for each metric.

switch f
	case {'Z', 'p'}
		dataI = interp1(...
			data.computeTimes / 1e3, ...
			data.(f), time, 'nearest');
	case 'x'
		dataI = interp1(...
			data.computeTimes / 1e3, ...
			data.V(1, :), time, 'nearest');
		f = 'Vx'; 
	case 'y'
		dataI = interp1(...
			data.computeTimes / 1e3, ...
			data.V(2, :), time, 'nearest');
		f = 'Vy'; 
end
			
end

function times = get_unique_wavetimes(data)
fields = fieldnames(data);
alltimes = cellfun(@(f) ...
	data.(f).computeTimes(:), ...
	fields, 'uni', 0);
alltimes = cat(1, alltimes{:});
times = unique(alltimes / 1e3, 'stable');
end

function args = parse_inputs(args)

P = inputParser;
p =@(varargin) addParameter(P, varargin{:});

% Defaults
metrics = {...
	'maxdescent', ...
	'events', ...
	'delays_T10_fband1_13', ...
	'delays_T01_fband1_13'}; 
allMetrics = {...
	'maxdescent', ...
	'events', ...
	'delays_T10_fband1_13', ...
	'delays_T01_fband1_13', ...
	'delays_T10_fband1_50', ...
	'delays_T01_fband1_50'}; 

p('pat', '*');
p('seizure', '*');
p('files', []);
p('metrics', metrics, @(c) all(contains(c, allMetrics)));
p('sig', 5e-2);

parse(P, args{:});

args = P.Results;

% Cleaning
if isnumeric(args.seizure), args.seizure = num2str(args.seizure); end
if isempty(args.files)
	if ischar(args.pat) && ischar(args.seizure)
		if strcmpi([args.pat args.seizure], '**')
			fid = fopen('seizures2.txt');
			A = textscan(fid, '%s %d'); 
			pats = A{1}; seizures = A{2};
			for ii = numel(pats):-1:1
				files(ii) = dir(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop.mat', pats{ii}, seizures(ii)));
			end
			fclose(fid);
			args.files = files;
		else
			args.files = dir(sprintf('%s_Seizure%s_Neuroport_10_10_wave_prop.mat', args.pat, args.seizure));
		end
	else
		
		ii = 1;
		for p = 1:numel(args.pat)
			for s = 1:numel(args.seizure)
				try 
					args.files{ii} = ...
						dir(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop.mat', args.pat{p}, args.seizure(s)));
					ii = ii + 1;
				catch ME
					if ~strcmpi(ME.identifier, 'MATLAB:matrix:singleSubscriptNumelMismatch')
						rethrow(ME)
					end
				end
			end
		end	
	end

end

end
