% filename = m.Name;
% mea = m.Neuroport;
% pat = 'c5'; seizure = 3;
patpath = genpath(pat);
addpath(patpath);

fname = sprintf('%s_Seizure%d_Neuroport', pat, seizure);
if ~exist([fname '.mat'], 'file')
	if exist([fname '_10_10.mat'], 'file')
		fname = [fname '_10_10.mat'];
		fprintf('Using %s\n', fname);
	elseif exist(sprintf('%s_Seizure%d.mat', pat, seizure), 'file')
		disp('Creating epoch file ...')
		create_epoch(pat, seizure, 'padding', [10 10]);
		fprintf('Using %s\n', fname);
		fname = [fname '_10_10.mat'];
	end
end
mea = matfile(fname);
[~, name, ~] = fileparts(mea.Properties.Source);

disp('Computing wave directions from delays...')
[delays, mea] = wave_prop(mea, 'delays');
plot_wave_directions(mea, 
wave_prop(mea, 'bos', true);
save([name '_wave_prop'], 'delays')
disp('Done.')

% disp('Saving figures...')
% nn = strrep(mea.Name, 'Seizure', '');
% print(2, nn, '-dpng')
% print(3, [nn '_FF'], '-dpng')

% disp('ALL DONE! GO HOME!!!')
% print(4, [nn, '_coh'], '-dpng')

rmpath(patpath);
