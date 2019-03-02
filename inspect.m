function [] = inspect(pat)
%%
% pat = 'CUCX3';
patPath = genpath(pat);
addpath(patPath);
SF = 100;
% try 
% 	map = neuroport_electrode_map();
% % 	map = metadata.map;
% catch MEmap
% 	md = jsondecode(fileread('metadata.json'));
% 	map = md.electrodemap;
% end

% try 
% 	start_times = [md.seizures.onset];
% catch MEonset
% 	metas = dir([pat filesep pat '*metadata.mat']);
% 	start_times = zeros(numel(metas), 1);
% 	for i = 1:numel(metas)
% 		md = load(metas(i).name);
% 		start_times(i) = md.metadata.seizure_start;
% 	end
% end

% seizures = struct('ns5', [], ...
% 	'epoch_type', [], ...
% 	'seizure_start', [], ...
% 	'seizure_end', [] ...
% 	);
micros = dir([pat filesep pat '_Neuroport' filesep '*.ns5']);
%%
for i = 1:numel(micros)
% 	seizures(i).ns5 = micros(i).name;

	nsx = openNSx(micros(i).name, 'skipfactor', SF);
%%
	samplingRate = nsx.MetaTags.SamplingFreq / SF;
	if nsx.RawData.PausedFile
% 		seizures.pausePoints = nsx.MetaTags.DataPoints;
		Data = cell(size(nsx.Data))';
		for di = 1:numel(nsx.Data)
			temp = single(nsx.Data{di})';
			temp = smoothdata(temp ./ max(abs(temp)), 'movmean', samplingRate / 10);
			Data{di} = zca_whitening(temp);
		end
		Data = cell2mat(Data);
	else
		Data = single(nsx.Data');
	end
	Data = diff(Data);
	[nSamp, nCh] = size(Data);
	if nCh > 96
		Data = Data(:, 1:96);
		nCh = 96;
	end
% 	Data = Data ./ max(abs(Data)) / 2 + int16(1 : nCh);
% 	Data = smoothdata(Data, 'gaussian', samplingRate / 10);
	Data = Data ./ max(abs(Data));
	T = (1:nSamp) / samplingRate;
	figure(998); fullwidth(); 
% 	plot((1:nSamp-1) / samplingRate, zscore(diff((Data))) + (1:nCh))
	yyaxis right
	plot(T, smoothdata(std(Data, [], 2), 'movmean', samplingRate), 'k', 'linewidth', 2); 
	yyaxis left
	set(gca, 'colororder', cool(nCh)); set(gca, 'nextplot', 'replacechildren')
	plot(T, (((Data))) + (1:nCh))
	axis tight;
	title([pat ' ' nsx.MetaTags.Filename]); xlabel('Time (s)'); ylabel('Channel')
	print(998, [pat filesep pat '_' nsx.MetaTags.Filename(1:end-4) '_full'], '-dpng')
	savefig(998, [pat filesep pat '_' nsx.MetaTags.Filename(1:end-4) '_full'])
% 	figure(999); fullwidth; 
% 	p1 = subplot(2, 2, [1 2]); set(p1, 'colororder', cool(nCh))
% 	p1.NextPlot = 'replaceChildren';  
% 	for ts = 1: round(30 * samplingRate) : nSamp
% 		S = 60;
% 		inds = (ts : ts + S * samplingRate);
% 		inds(inds > nSamp) = [];
% 		data_temp = Data(inds, :);
% 		plot(p1, inds / samplingRate, data_temp)
% 		p1.ColorOrder = cool(nCh);
% 		axis('tight');
% 		ylim([-1 1])
% 		title([pat ' ' nsx.MetaTags.Filename]);
% 		xticks(p1, round(inds(1:2*samplingRate:end) / samplingRate))
% 		subplot(223);
% 		[X, Y] = ndgrid(inds / samplingRate, 1:nCh);
% 		mesh(X, Y, data_temp); zlim([-1 1])
% 		subplot(224)
% 		plot(inds / samplingRate, median(data_temp, 2)); axis('tight'); ylim([-1 1])
% 		print(gcf, [pat filesep pat '_' nsx.MetaTags.Filename(1:end-4) '_' num2str(ts / samplingRate, '%04.0f')], '-dpng')
% % 		pause();
% % 		close(gcf)
% 	end
% 	figure(1); fullwidth(1); 
% 	plot((1:nSamp) / samplingRate, ...
% 		Data ./ max(abs(Data) / 2) + int16((1:nCh)));
% 	seizures(i).epoch_type = input('Epoch type (ictal/non): ', 's');
% 	if strcmpi(seizures(i).epoch_type, 'ictal')
% 		fprintf('Proposed starts:\n')
% 		if nsx.RawData.PausedFile
% 			fprintf('Pauses at ');
% 			fprintf('%0.4f  ', nsx.MetaTags.DataDurationSec(1:end-1));
% 			fprintf('\n');
% 		end
% 		fprintf('%0.4f  ', start_times); 
% 		fprintf('\n');
% 		seizures(i).seizure_start = input('Seizure start: ');
% 		seizures(i).seizure_end = input('Seizure end: ');
% 	end
% pause();
	
end
% metadata.patient = pat;
% metadata.mea = seizures;
% metadata.map = map;
% metadata.samplingRate = nsx.MetaTags.SamplingFreq;
% fid = fopen([pat '_metadata.json'], 'w');
% fwrite(fid, jsonencode(metadata));
% fclose(fid);

rmpath(patPath);