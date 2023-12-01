%{
VARIABLES

    stencil_ind (21x1)          - indices of stencil_legit_neigh (needed
                                  for picking destinations
    neigh_setup (3x3x3)         - matrix to check eligible neighbors to
                                  spread to)
    no_neigh_count              - check to see how often there was no
                                  eligible voxel to spread to
    new_type                    - state of the voxel that was chosen to
                                  spread
    step                        - spreading step to the occupied voxel
    legit (3x3x3)               - bool matrix of eligible neighbors
    legit_x/_y/_z (0-6 each)    - (3x3x3) indices of eligible neighbors
    legit_steps (0-6x3)         - steps to eligible neighbors
    step (3)                    - randomly chosen step (from eligible ones)                  

%}
function step = choose_occupied_voxel(new_type, neigh_setup, stencil_ind, NB)

    global no_neigh_count
    global spreads
    
    if strcmp(NB, 'M_2')
        ind_corr = 3;
    else
        ind_corr = 2;
    end
    
    switch new_type
        
        case 0

            legit = double(neigh_setup == 0) + double(neigh_setup == 1);
            legit(stencil_ind) = 0;
            [legit_x,legit_y,legit_z] = ind2sub(size(legit), find(legit));
            if legit_x
                % transform indices (1,2,3) to steps (-1,0,1)
                legit_steps = [legit_x, legit_y, legit_z] - ind_corr;
                step = legit_steps(randi(size(legit_steps,1)),:);
                spreads(1,1) = spreads(1,1) + 1;
            else
                step = [0,0,0];
                no_neigh_count = no_neigh_count + 1;
            end

    
        case 1

            %legit = double(neigh_setup == 1) + double(neigh_setup == 3);
            legit = double(neigh_setup == 1) + double(neigh_setup == 2) + double(neigh_setup == 3);
            legit(stencil_ind) = 0;
            [legit_x,legit_y,legit_z] = ind2sub(size(legit), find(legit));
            if legit_x
                % transform indices (1,2,3) to steps (-1,0,1)
                legit_steps = [legit_x, legit_y, legit_z] - ind_corr;
                step = legit_steps(randi(size(legit_steps,1)),:);
                spreads(2,1) = spreads(2,1) + 1;
            else
                step = [0,0,0];
                no_neigh_count = no_neigh_count + 1;
            end

        case 2

            %legit = double(neigh_setup == 2) + double(neigh_setup == 0);
            legit = double(neigh_setup == 1) + double(neigh_setup == 2) + double(neigh_setup == 0);
            legit(stencil_ind) = 0;
            [legit_x,legit_y,legit_z] = ind2sub(size(legit), find(legit));
            if legit_x
                % transform indices (1,2,3) to steps (-1,0,1)
                legit_steps = [legit_x, legit_y, legit_z] - ind_corr;
                step = legit_steps(randi(size(legit_steps,1)),:);
                spreads(3,1) = spreads(3,1) + 1;
            else
                step = [0,0,0];
                no_neigh_count = no_neigh_count + 1;
            end

        case 3

            legit = double(neigh_setup == 2) + double(neigh_setup == 3);
            legit(stencil_ind) = 0;
            [legit_x,legit_y,legit_z] = ind2sub(size(legit), find(legit));
            if legit_x
                % transform indices (1,2,3) to steps (-1,0,1)
                legit_steps = [legit_x, legit_y, legit_z] - ind_corr;
                step = legit_steps(randi(size(legit_steps,1)),:);
                spreads(4,1) = spreads(4,1) + 1;
            else
                step = [0,0,0];
                no_neigh_count = no_neigh_count + 1;
            end
        
    end

end