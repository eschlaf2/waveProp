% filename = m.Name;
% mea = m.Neuroport;
% pat = 'c5'; seizure = 3;
% patpath = genpath(pat);
% addpath(patpath);
% 
% fname = sprintf('%s_Seizure%d_Neuroport', pat, seizure);
% if ~exist([fname '.mat'], 'file')
% 	if exist([fname '_10_10.mat'], 'file')
% 		fname = [fname '_10_10.mat'];
% 		fprintf('Using %s\n', fname);
% 	elseif exist(sprintf('%s_Seizure%d.mat', pat, seizure), 'file')
% 		disp('Creating epoch file ...')
% 		create_epoch(pat, seizure, 'padding', [10 10]);
% 		fprintf('Using %s\n', fname);
% 		fname = [fname '_10_10.mat'];
% 	end
% end
% mea = matfile(fname);
% [~, name, ~] = fileparts(mea.Properties.Source);

mea = load('SIM/seizing_cortical_field_sim.mat');
name = mea.Name;
outfile = matfile([name '_wave_prop'], 'writable', true);

% disp('Computing wave directions from delays ...')
% [delays, mea] = wave_prop(mea, 'delays');
% plot_wave_directions(mea, delays);
% print(gcf, delays.Name, '-dpng')
% outfile.delays = delays;

disp('Computing wave directions from events ...')
[events, mea] = wave_prop(mea, 'events');
plot_wave_directions(mea, events);
print(gcf, events.Name, '-dpng');
outfile.events = events;

disp('Computing wave directions from maxdescent ...')
[maxdescent, mea] = wave_prop(mea, 'maxdescent');
plot_wave_directions(mea, maxdescent);
print(gcf, maxdescent.Name, '-dpng');
outfile.maxdescent = maxdescent;

disp('Computing wave directions from deviance ...')
[dev, mea] = wave_prop(mea, 'deviance');
plot_wave_directions(mea, dev);
print(gcf, dev.Name, '-dpng');
outfile.dev = dev;

disp('Done.')

rmpath(patpath);
