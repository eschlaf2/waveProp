function save_scm_ecog_data(params)
%% Create example 10 s of ECoG data to analyze wave dynamics.

	%%%% Specify directory where simulation results were saved ---------------
	struct2var(params);

	%%%% Choose the map type, and specify a starting index to load. -----------
	%%%% Load 10 s of simulated ECoG and create a new data variable -----------
	duration=90;                                           %# of intervals to load.
	T=5000;                                         %Time indices per interval.
	N=9;                                            %# electrodes.

	ECOG = zeros(duration*T, N);                           %Variable to hold ECoG.
	phi = zeros(duration, 3, 3);
	for k=1:duration                                       %For each time point,
		fprintf(['Read in ' num2str(k0+k) '\n'])    %... load in simulation, 
		load(sprintf('%s_%d_%03d.mat', basename, sim_num, k0 + k*T0), 'NP', 'last')
		ECOG(1+(k-1)*T:k*T,:)=reshape(NP.Qe, T, N);                %... and store the result.
		phi(k,:, :)=last.phi_ee(49:51, 49:51);
		figure(1); fullwidth(); subplot(121); imagesc(last.Qe); subplot(122); imagesc(squeeze(NP.Qe(1, :, :))); title(sprintf('k = %d', k+k0))
		drawnow()
		pause(.01)
	end

	%%%% Downsample the simulation data from to 500 Hz. -----------------------
	dec = 10;
	data = zeros(duration*T/10, N);
	for k=1:N                                       %For each channel,
		d0 = ECOG(:,k);                             %... get the ECoG,
		d0 = decimate(d0,dec);                      %... downsample it,
		data(:,k) = d0;                             %... and store the result.
	end

	dt = (time(10)-time(9))*dec;                    %Define sampling interval.
	fs = 1/dt;                                      %Define sampling frequency.

	%%%% Define the positions. ------------------------------------------------
	y = [-12, -12, -12,   0,  0,  0,  12, 12, 12];
	x = [-12,   0,  12, -12,  0, 12, -12,  0, 12];
	position = [x;y]';

	%%%% Save the results. ----------------------------------------------------
	% NOTE: This file can be loaded and analyzed in analysis/main_seizure_wave.m
	save(sprintf('%s_%d_ecogwaves', basename, sim_num), 'data', 'fs', 'position');

end
