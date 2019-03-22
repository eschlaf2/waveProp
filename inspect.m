function [] = inspect(pat, files)
% Creates a PNG of an entire ns5 file for the indicated patient. If no
% specific files are identified, then all NS5 files in the given patient's
% folder are inspected. The raw data are whitened and smoothed and then the
% first derivative is taken. This is then sent to plot_Neuroport() to be
% plotted. 

patPath = genpath(pat);
addpath(patPath);
SF = 100;

micros = dir([pat filesep pat '_Neuroport' filesep '*.ns5']);

if ~exist('files', 'var') || isempty(files)
	files = 1:numel(micros);
else
	files(files > numel(micros)) = [];
	files(files < 1) = [];
end

%%
for i = files
	disp(micros(i).name);
	nsx = openNSx(micros(i).name, 'skipfactor', SF);
%%
	samplingRate = nsx.MetaTags.SamplingFreq / SF;
	if nsx.RawData.PausedFile
		nsx.Data = nsx.Data(logical(nsx.MetaTags.Timestamp));
		Data = cell(size(nsx.Data))';
		for di = 1:numel(nsx.Data)
			temp = single(nsx.Data{di})';
			temp = smoothdata(temp ./ max(abs(temp)), 'movmean', samplingRate / 10);
			Data{di} = zca_whitening(temp);
		end
		Data = cell2mat(Data);
	else
		Data = single(nsx.Data') + 1e-6;
		temp = smoothdata(Data ./ max(abs(Data)), 'movmean', samplingRate / 10);
		Data = zca_whitening(temp);
	end
	Data = diff(Data);
	plot_Neuroport(Data, samplingRate);
	title([pat ' ' nsx.MetaTags.Filename]); xlabel('Time (s)'); ylabel('Channel')
	print(998, [pat filesep pat '_' nsx.MetaTags.Filename(1:end-4) '_full'], '-dpng')
	try
		savefig(998, [pat filesep pat '_' nsx.MetaTags.Filename(1:end-4) '_full'])
	catch ME
		display(ME)
	end
	close 998;

	
end

rmpath(patPath);