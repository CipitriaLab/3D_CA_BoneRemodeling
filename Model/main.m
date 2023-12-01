%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    3D stochastic EGT CA for Bone Remodeling             %
%                             Anna-Dorothea Heller                        %
%                                   Dez 2023                              %
%                                Matlab R2023b                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n       - size of domain (side length of cube)
% type    - different settings for initial domains, customizable in
%           create_setup.m
%           > 'quiet_resorp_X': quiescent cube with resorption patch
%                             positioned on layer X (X=2 -> surface,
%                             X>2 -> deeper within cube)
%           > 'quiet_form_X': quiescent cube with formation patch
%                             positioned on layer X (X=2 -> surface,
%                             X>2 -> deeper within cube)
% NB      - different neighborhood settings, customizable in
%           create_setup.m; possibly choose_expansion_voxel.m,
%           choose_occupied_voxel.m and BR_model.m have to be adapted as
%           well
%           > 'vN_1': von-Neumann neighborhood with radius 1 
%           > 'M_1': Mooore neighborhood with radius 1
%           > 'M_2': Mooore neighborhood with radius 2
% iter    - number of iterations
% model   - simulation modes
%           > 'regular': run model start to end with given initial
%                        conditions and iterations
%           > 'rm_cycle': 1. model resorption for number of iterations (iter1)
%                            > make sure to set 'quiet_resorp_X' as type
%                         2. seed number of Formation voxels (seed_F)
%                         3. continue simulation for number of iterations
%                            (iter2 = iter-iter1)
% further model variables described in BR_model.m

clear all;
tic
% If requested, define random seed by uncomment following statement

%rand_seed = 1;
%rng(rand_seed);

%% Initial Settings

% Initial Config
n = 30;
type = 'quiet_resorp_3';
NB = 'M_2';

g_EE = 0.;
g_RE = 0.;
g_FE = 0.;
g_QE = 0.;

g_ER = 0.5; 
g_RR = 0.6; 
g_FR = 0.;
g_QR = 0.;

g_EF = 0.;
g_RF = 0.;
g_FF = 0.;
g_QF = 0.;

g_EQ = 0.;
g_RQ = 0.5;
g_FQ = 0.;
g_QQ = 0.;

G = [g_EE g_ER g_EF g_EQ ; g_RE g_RR g_RF g_RQ ; g_FE g_FR g_FF g_FQ ; g_QE g_QR g_QF g_QQ];

% Iterations
iter = 5000;

% Simulation mode
model = 'regular';

% Create Set-up
[M_init, g_bar, stencil_ind, neighb, def_neighbourX,...
 def_neighbourY, def_neighbourZ] = create_setup(n, type, G, NB);

%% Run Model

switch model

    case 'regular'

        [M_RFQ, count, time_count, field_exp, time_field, iter, n, G, t, spreads, no_neigh_count] = ...
        BR_model(M_init, g_bar, G, stencil_ind, neighb, def_neighbourX,...
        def_neighbourY, def_neighbourZ, iter, n, NB);
        
        % save results

        formatOut = 'ddmmyy_HHMM';
        save(join(['res', datestr(now,formatOut)]), 'count', 'time_count', 'field_exp', 'time_field', 'iter', 'n', 'G', 't', 'spreads', 'no_neigh_count');

    case 'rm_cycle'

        iter1 = 500;
        iter2 = iter-iter1;
        seed_F = 9;

        % Resorption
        M_init1 = M_init;

        [M_RFQ1, count1, time_count1, field_exp1, time_field1, iter1, n, G, t1, spreads1, no_neigh_count1] = ...
            BR_model(M_init1, g_bar, G, stencil_ind, neighb, def_neighbourX,...
            def_neighbourY, def_neighbourZ, iter1, n, NB);
        
        % seed Formation
        
        void_voxel = M_RFQ1 == 0;
        
        inner_voxels = zeros(n, n, n);
        inner_voxels(2:end-1, 2:end-1, 2:end-1) = ones(n-2, n-2, n-2);
        
        [i_v_i, i_v_j, i_v_k] = ind2sub(size(inner_voxels), find(void_voxel & inner_voxels));
        inner_voids = [i_v_i, i_v_j, i_v_k];
        inner_voids_sort = sortrows(inner_voids, 2, 'descend');
        
        if length(inner_voids_sort)>seed_F
           fill_form = sub2ind(size(inner_voxels),inner_voids_sort(1:15,1), inner_voids_sort(1:15,2), inner_voids_sort(1:15,3));
        else
           fill_form = sub2ind(size(inner_voxels),inner_voids_sort(:,1), inner_voids_sort(:,2), inner_voids_sort(:,3));
        end
        
        M_init2 = M_RFQ1;
        M_init2(fill_form) = 2;
        
        
        % Formation
        
        [M_RFQ2, count2, time_count2, field_exp2, time_field2, iter2, n, G, t2, spreads2, no_neigh_count2] = ...
            BR_model(M_init2, g_bar, G, stencil_ind, neighb, def_neighbourX,...
            def_neighbourY, def_neighbourZ, iter2, n, NB);
        
        % save results

        % 1 stands for the Resorption (1st) round of the cycle
        % 2 stands for the Formation (2nd) round of the cycle
        % > depending on evaluation purposes, further varibales can be
        %   grouped together as shown done below and/or included to the
        %   saved file

        count = [count1, count2];
        time_count = [time_count1, time_count2];
        time_field = [time_field1, time_field2];

        formatOut = 'ddmmyy_HHMM';
        save(join(['res', datestr(now,formatOut)]), 'count', 'time_count', 'time_field', 'n', 'G',...
         'iter1', 'iter2', 'spreads1', 'spreads2', 'no_neigh_count1', 'no_neigh_count2',...
          'M_init1', 'M_init2', 'M_RFQ1', 'M_RFQ2', 'field_exp1', 'field_exp2');


end

toc


