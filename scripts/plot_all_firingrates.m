fid = fopen('seizures.txt'); seizures = textscan(fid, '%s %d'); fclose(fid);
r = 3; c = 4;
fr = matfile('firingrates', 'Writable', true);

% for ii = 1:3, figure(ii); fullwidth(true); end

for ii = 1:31
%     pat = seizures{1}{ii};
%     seiz = seizures{2}(ii);
%     varname = sprintf('%s_%d', pat, seiz);
%     mea = load(sprintf('%s%s%s_Seizure%d_Neuroport_10_10.mat', ...
%         pat, filesep, pat, seiz)); 
%     temp = mean(mua_firing_rate(mea), 2); 
%     fr.(varname) = temp;
    name = strrep(fields{ii}, '_', ' ');
    figure(ceil(ii / 12));
    subplot(r, c, mod(ii - 1, 12) + 1); 
    plot(fr.(fields{ii})), 
    title(name); 
end


