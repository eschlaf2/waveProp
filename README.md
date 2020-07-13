# Seizure simulation and analysis


## Creating simulations

1. Create a simulation using the default parameters:
```
sim = SCM; 
sim.Run();
sim.Preview();
```
This creates a simulation with a 60 second seizure and 10 seconds of pre- and post-ictal activity. The result is a file for each second of the simulation (`./SCM/SCM/SCM_X_TTT.mat`) and a simulated MEA recording (`./SCM/SCM_SeizureX_10_10.mat`). The parameters used are saved to `./SCM/SCM/SCM_X_info.mat`. 

2. Visualize the simulation while running it:
```
sim = SCM('visualization_rate', 10, 'padding', [0 0], 'duration', 3, 'save', 0);
sim.Run();
```
This creates a simulation with a 3 second seizure and no pre- or post-ictal activity. No files are saved in this case.

3. Visualize a different set of states and save the results:
```
P = struct();
P.visualization_rate = 10;
P.duration = 3;
P.padding = [0 0];
P.out_vars = {'Qe', 'map', 'K'};
P.sim_num = 99;
sim = SCM(P);
sim.Run();
```

4. Continue or extend a simulation:
```
load('SCM/SCM/SCM_99_info', 'params')
params.duration = 5;
params.t0_start = 3;
sim = SCM(params);
sim.Run();
```

## Viewing a simulated micro electrode array (MEA) recording

1. Load a simulated MEA recording
```
mea = MEA('SCM/SCM_Seizure1_Neuroport_10_10.mat');
mea.preview(30:.01:31);
lfp = mea.lfp;
plot(mea.Time, mean(mea.lfp, 2));

```


Martinet, L. E., G. Fiddyment, J. R. Madsen, E. N. Eskandar, W. Truccolo, U. T. Eden, S. S. Cash, and M. A. Kramer. “Human Seizures Couple across Spatial Scales through Travelling Wave Dynamics.” Nature Communications 8 (2017). https://doi.org/10.1038/ncomms14896.

Steyn-Ross, Moira L, D A Steyn-Ross, and J W Sleigh. “Interacting Turing-Hopf Instabilities Drive Symmetry-Breaking Transitions in a Mean-Field Model of the Cortex: A Mechanism for the Slow Oscillation,” n.d. https://doi.org/10.1103/PhysRevX.3.021005.





