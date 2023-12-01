# 3D_CA_BoneRemodeling
A 3D stochastic cellular automaton (CA) model governed by evolutionary game theory rules to simulate local dynamics of the microenvironment that contains four different voxels states: Formation, Quiescence, Resorption and Environment.

Conceptualized by Anna-Dorothea Heller, Angelo Valleriani and Amaia Cipitria, written by Anna-Dorothea Heller. 

When used please reference:

    not published, yet


Tested for Matlab 2020b, 2021b, 2022b and 2023b on MacOS and GNU/Linux.

"Model" contains:

- main.m    >    to start single simulation
- main_batch.m    >    to start multiple simulations with different parameters

- model functions: create_setup.m, BR_model.m, choose_expansion_voxel.m, choose_occupied_voxel.m

- comb_para.m    >    to create parameter list for main_batch.m
- comb_pair.m    >    to create parameter list of pre-defined pairs for main_batch.m
