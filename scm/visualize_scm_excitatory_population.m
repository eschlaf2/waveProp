function visualize_scm_excitatory_population(params)
%% Visualize the activity of the excitatory population. -----------
	%%%% Specifiy directory where simulation results were saved. --------------
	struct2var(params);
	%%%% Choose the map type, and specify total number of indices. ------------
	%%%% Visualize the excitatory population activity. ------------------------
	figure(); fullwidth(true);
	counter=1;
	for k=1:2:(sum(padding) + duration)                                     %For each 1 s interval, in steps of 2 s,
		fprintf(['Read in ' num2str(k*T0) '\n'])       %... print counter,
													%... load the saved data,
		load(sprintf('%s_%d_%03d.mat', basename, sim_num, k*T0), 'last')
		Q0 = last.Qe;                               %... get the excitatory activity at last moment of time,
	    subplot(10,10,counter)                      %... assign a subplot
		imagesc(Q0, [0 25])                         %... and image the excitatory activity with range [0 25] Hz.
		colormap bone;
% 		mov(counter) = getframe(gcf);
		axis off                                    %... exclude axes for visualization.
		counter=counter+1;                          %... augment the counter.
		if counter > 100, continue, end
	end

end