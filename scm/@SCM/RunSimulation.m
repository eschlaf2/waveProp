function RunSimulation(params)

PM = params.meta;

% Extract variables from meta-parameters
BASENAME = PM.basename;
DURATION = PM.duration;
PADDING = PM.padding;
SAVE = PM.save;
SIM_NUM = PM.sim_num;
T_STEP = PM.t_step;
T0_START = PM.t0_start;

% Create stimulus map
if T0_START > -PADDING(1)
	% If not starting a fresh sim, use previously saved <last>
    try
        load(sprintf('%s_%d_%03d', BASENAME, SIM_NUM, T0_START - 1), 'last')
    catch ME
        if ~(strcmpi('MATLAB:load:couldNotReadFile', ME.identifier))
            rethrow(ME);
        else
            last = params.IC;
        end
    end
else
	% ... otherwise, start a fresh sim
	last = params.IC;  %Load the initial conditions to start.
end

K = PADDING(2) + DURATION;
fig = [];
for t0 = PM.t0_start:T_STEP:K-T_STEP  % For each step
	tic
	% Update time offset
	params.t0 = t0;
	
	% ... show progress, 
	fprintf('Running %d / %d .. ', t0+1, floor(K));  
		
	% ... run simulation for duration T_STEP,
	[NP, EC, time, last, fig] = ...  
		seizing_cortical_field('legacy variable', ...
            min(T_STEP, K - t0), last, fig, params);
	
	% ... save the results of this run,
	if SAVE
		fprintf('Saving .. ')
		fname = checkname(sprintf('%s_%d_%03d', BASENAME, SIM_NUM, t0*T_STEP));
		save(fname, 'NP','EC','time','last');
	end
	toc
	% ... update progress.
	fprintf('Done.\n')  
end

end
