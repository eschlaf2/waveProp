# Seizure simulation and analysis


## Creating simulations

1. Create a simulation using the default parameters:
```
params = SCMParams; 
scm_sim;
```
This creates a simulation with a 60 second seizure and 10 seconds of pre- and post-ictal activity. The result is a file for each second of the simulation (`./SCM/SCM/SCM_X_TTT.mat`) and a simulated MEA recording (`./SCM/SCM_SeizureX_10_10.mat`). The parameters used are saved to `./SCM/SCM/SCM_X_info.mat`. 

2. Visualize the simulation while running it:
```
params = SCMParams('visualization_rate', 10, 'padding', [0 0], 'duration', 3, 'save', 0);
scm_sim;
```
This creates a simulation with a 3 second seizure and no pre- or post-ictal activity. No files are saved in this case.

3. Visualize a different set of states and save the results:
```
P.visualization_rate = 10;
P.duration = 3;
P.padding = [0 0];
P.out_vars = {'Qe', 'map', 'K'};
P.sim_num = 99;
params = SCMParams(P);
scm_sim;
```

4. Continue or extend a simulation:
```
load('SCM/SCM/SCM_99_info', 'params')
params.duration = 5;
params.t0_start = 3;
scm_sim;
```

