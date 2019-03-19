pat = 'c7';
files = dir(['/Volumes/elements/' pat '/' pat '_seizure*metadata.mat']);
md.patient = pat;
i = 0;
for f = files'
	i = i + 1;
	load([f.folder filesep f.name]);
	metadata.seizure = str2double(f.name(11));
	if ~isfield(md, 'seizures')
		md.seizures = metadata;
	else
		md.seizures = [md.seizures; metadata];
	end
end

fid = fopen([f.folder filesep pat '_metadata_original.json'], 'w');
fprintf(fid, '%s', jsonencode(md));
fclose(fid);
md = rmfield(md, 'seizures');


% load(['/Volumes/Elements/' pat '/' pat '_seizure _metadata.mat', 'metadata')
% metadata.seizure = 3;
% % md = metadata;
% md = [md; metadata];