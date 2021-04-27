function [files, names] = txt2files(fname)
if nargin < 1 || isempty(fname), fname = 'seizures2.txt'; end
fid = fopen(fname);
A = textscan(fid, '%s %d'); 
pats = A{1}; seizures = A{2};
for ii = numel(pats):-1:1
	files(ii) = dir(sprintf('%s_Seizure%d_fits.mat', ...
	    pats{ii}, seizures(ii)));
end
fclose(fid);
fname = [files(ii).folder filesep files(ii).name];
names = arrayfun(@(ii) [pats{ii} ' ' num2str(seizures(ii))], 1:numel(pats), 'uni', 0);
