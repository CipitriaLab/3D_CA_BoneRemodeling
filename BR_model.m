%{
VARIABLES:

    n                           - initial size of cube
    iter                        - number of all iterations
    t                           - iteration step 
    time                        - current simulation time depending on
                                  update steps
    time_count                  - actual time vector for cell count tracking
                                  depending on update steps
    time_field                  - actual time vector for field tracking
                                  depending on update steps
    def_neighbourX/Y/Z (1x3)    - possible steps per direction
    neighb (6x3)                - collected possible neighbours
    stencil_ind (21x1)          - indices of stencil_legit_neigh (needed
                                  for picking destinations
    G (4x4)                     - payoff matrix
    omega                       - scaling of payoff matrix impact
    g_bar (16x1)                - transformed payoff parameters
    M_init (nxnxn)              - initial configuration
    M_age (nxnxn)               - age property of voxels
    M_RFQ_B ((n+2)x(n+2)x(n+2)) - domain with 2-voxel environment layer to
                                  model space to spread in
    Fit (nxnxn)                 - matrix with expansion rate of each voxel   
    neigh_setup (3x3x3)         - matrix to check eligible neighbors to
                                  spread to)
    no_neigh_count              - check to see how often there was no
                                  eligible voxel to spread to
    new_type                    - state of the voxel that was chosen to
                                  spread
    step                        - spreading step to the occupied voxel
    old_type                    - state of the voxel that was chosen to
                                  be occupied
    min_val                     - time that update step takes
    update (3x3)                - x-, y- & z-indices of the domain to be updated 
    filled_voxels               - domain voxels, that are not void
    
%}


function [M_RFQ, count, time_count, field_exp, time_field, iter, n, G, t, spreads, no_neigh_count]...
    = BR_model(M_init, g_bar, G, stencil_ind, neighb,...
    def_neighbourX, def_neighbourY, def_neighbourZ, iter, n, NB)

    %% Intialize model variables

    M_RFQ = M_init;
    
    if strcmp(NB, 'M_2')       
        M_RFQ_B = cat(3, zeros(n+6, n+6, 3),...
                [zeros(3, n+6, n); zeros(n,3,n), M_RFQ, zeros(n,3,n); zeros(3, n+6, n)],...
                 zeros(n+6, n+6, 3));    
    else
        M_RFQ_B = cat(3, zeros(n+4, n+4, 2),...
                        [zeros(2, n+4, n); zeros(n,2,n), M_RFQ, zeros(n,2,n); zeros(2, n+4, n)],...
                         zeros(n+4, n+4, 2));
    end

    Fit = zeros(size(M_RFQ));
    M_age = zeros(size(M_init));
    
    update = true;
    
    update_x = [];
    update_y = [];
    update_z = [];

    time = 0;

    %% Initialize evaluation data

    [count, field_exp, time_count, time_field, step_count, step_field] = init_eval_data(iter);

    count(:,1) = [sum(sum(sum(M_RFQ(:,:,:)==1)));...
                 sum(sum(sum(M_RFQ(:,:,:)==2)));...
                 sum(sum(sum(M_RFQ(:,:,:)==3)))];
    time_count(1) = time;
    field_exp{1} = M_RFQ(:,:,:);
    time_field(1) = time;
    % how often did each of the strategy spread
    global spreads
    spreads = zeros(4,1);

    % check how often a state does not have a possible neighbor to move to
    global no_neigh_count
    no_neigh_count = 0;
    
    %% First Step
    
    t = 1;

    [min_val, min_x, min_y, min_z, Fit] = choose_expansion_voxel(n, t, neighb, M_RFQ_B, Fit, g_bar, update, NB);

    % Proliferation
    
    if strcmp(NB, 'M_2')
        
        new_type = M_RFQ_B(min_x+3, min_y+3, min_z+3);
        neigh_setup = M_RFQ_B(min_x+1:min_x+5, min_y+1:min_y+5, min_z+1:min_z+5);

    else
        
        new_type = M_RFQ_B(min_x+2, min_y+2, min_z+2);
        neigh_setup = M_RFQ_B(min_x+1:min_x+3, min_y+1:min_y+3, min_z+1:min_z+3);
    
    end
    
    step = choose_occupied_voxel(new_type, neigh_setup, stencil_ind, NB);

    % position in bone domain M_RFQ (without void border)
    new_pos = [min_x, min_y, min_z] + step;
    
    %defines, whether new positions need to be assigned or not   
    update = true;

    if any(new_pos < 1) || any(new_pos > n)
        
        if new_type == 0
            
            update = false;
            
        else
            
            [M_RFQ, M_RFQ_B, Fit, M_age, n, new_pos] = update_borders(M_RFQ, Fit, M_age, n, new_pos, NB);
            
        end
        
    end
    
    if update
        
        if strcmp(NB, 'M_2')
            
            old_type = M_RFQ_B(new_pos(1)+3, new_pos(2)+3, new_pos(3)+3);
            M_RFQ_B(new_pos(1)+3, new_pos(2)+3, new_pos(3)+3) = new_type;
            M_RFQ = M_RFQ_B(4:end-3,4:end-3,4:end-3);
            
        else
            
            old_type = M_RFQ_B(new_pos(1)+2, new_pos(2)+2, new_pos(3)+2);
            M_RFQ_B(new_pos(1)+2, new_pos(2)+2, new_pos(3)+2) = new_type;
            M_RFQ = M_RFQ_B(3:end-2,3:end-2,3:end-2);
        
        end

        M_age = M_age + min_val;
        time = time + min_val;
        M_age(new_pos(1), new_pos(2), new_pos(3)) = 0;
    end


    %% Results (only if iter small enough)

    if iter <= 1e5

            %filled_voxels = sum(sum(sum(M_RFQ(:,:,:)>0)));
            count(:,2) = [sum(sum(sum(M_RFQ(:,:,:)==1)));...
                        sum(sum(sum(M_RFQ(:,:,:)==2)));...
                        sum(sum(sum(M_RFQ(:,:,:)==3)))];
            time_count(2) = time;

    end

    if iter <= 2*1e3

           field_exp{2} = M_RFQ(:,:,:);
           time_field(2) = time;

    end


    %% Iterations
    
    for t = 3:iter
        
        if n > 46
            break
        end
        
        if update 
            update_x = def_neighbourX + new_pos(1);
            update_x(update_x <= 0) = [];
            update_x(update_x >= n+1) = [];

            update_y = def_neighbourY + new_pos(2);
            update_y(update_y <= 0) = [];
            update_y(update_y >= n+1) = [];

            update_z = def_neighbourZ + new_pos(3);
            update_z(update_z <= 0) = [];
            update_z(update_z >= n+1) = [];

            Fit(update_x,update_y,update_z) = 0;
        end

        [min_val, min_x, min_y, min_z, Fit] = choose_expansion_voxel(n, t, neighb, M_RFQ_B, Fit, g_bar, update, NB, update_x, update_y, update_z);
    
        % Proliferation
         
        if strcmp(NB, 'M_2')
            
            new_type = M_RFQ_B(min_x+3, min_y+3, min_z+3);
            neigh_setup = M_RFQ_B(min_x+1:min_x+5, min_y+1:min_y+5, min_z+1:min_z+5);
        
        else
    
            new_type = M_RFQ_B(min_x+2, min_y+2, min_z+2);
            neigh_setup = M_RFQ_B(min_x+1:min_x+3, min_y+1:min_y+3, min_z+1:min_z+3);

        end
        
        step = choose_occupied_voxel(new_type, neigh_setup, stencil_ind, NB);

        % position in bone domain M_RFQ (without void border)
        new_pos = [min_x, min_y, min_z] + step;
        
        update = true;

        if any(new_pos < 1) || any(new_pos > n)
            
            if new_type == 0
                
                update = false;
                
            else 
                
                [M_RFQ, M_RFQ_B, Fit, M_age, n, new_pos] = update_borders(M_RFQ, Fit, M_age, n, new_pos, NB);
            
            end

        end
        
        if update
            
            if strcmp(NB, 'M_2')
                
                old_type = M_RFQ_B(new_pos(1)+3, new_pos(2)+3, new_pos(3)+3);
                M_RFQ_B(new_pos(1)+3, new_pos(2)+3, new_pos(3)+3) = new_type;
                M_RFQ = M_RFQ_B(4:end-3,4:end-3,4:end-3);

            else
                
                old_type = M_RFQ_B(new_pos(1)+2, new_pos(2)+2, new_pos(3)+2);
                M_RFQ_B(new_pos(1)+2, new_pos(2)+2, new_pos(3)+2) = new_type;
                M_RFQ = M_RFQ_B(3:end-2,3:end-2,3:end-2);
            
            end

            M_age = M_age + min_val;
            time = time + min_val;
            M_age(new_pos(1), new_pos(2), new_pos(3)) = 0;
        end

        %% Results
        
        if strcmp(NB, 'M_2')
            M_RFQ = M_RFQ_B(4:end-3,4:end-3,4:end-3);
        else
            M_RFQ = M_RFQ_B(3:end-2,3:end-2,3:end-2);
        end

        if iter <= 1e5

            %filled_voxels = sum(sum(sum(M_RFQ(:,:,:)>0)));
            count(:,t) = [sum(sum(sum(M_RFQ(:,:,:)==1)));...
                        sum(sum(sum(M_RFQ(:,:,:)==2)));...
                        sum(sum(sum(M_RFQ(:,:,:)==3)))];
            time_count(t) = time;

        elseif mod(t,step_count) == 0

            %filled_voxels = sum(sum(sum(M_RFQ(:,:,:)>0)));
            count(:,fix(t/step_count)) = [sum(sum(sum(M_RFQ(:,:,:)==1)));...
                        sum(sum(sum(M_RFQ(:,:,:)==2)));...
                        sum(sum(sum(M_RFQ(:,:,:)==3)))];
            time_count(fix(t/step_count)) = time;

        end

        if iter <= 2*1e3

           field_exp{t} = M_RFQ(:,:,:);
           time_field(t) = time;

        elseif mod(t,step_field) == 0    

           field_exp{fix(t/step_field)} = M_RFQ(:,:,:);
           time_field(fix(t/step_field)) = time;

        end

    end 
    
    %% Functions

    function [count, field_exp, time_count, time_field, step_count, step_field] = init_eval_data(iter)


        if iter <= 1e5
            count = zeros(3,iter);
            time_count = zeros(1,iter);
        else
            count = zeros(3,1e5);
            time_count = zeros(1,1e5);
        end

        if iter <= 2*1e3
            time_field = zeros(1, iter);
        else
            time_field = zeros(1, 2*1e3);
        end

        field_exp = {};

        step_count = iter./1e5;

        step_field = iter./(2*1e3);

    end

    function [M_RFQ, M_RFQ_B, Fit, M_age, n, new_pos] = update_borders(M_RFQ, Fit, M_age, n, new_pos, NB)

        if strcmp(NB, 'M_2')
            
            if any(new_pos == -1) || any(new_pos == n+2)
                
                1
                
                M_RFQ_B = cat(3, zeros(n+10,n+10,5),...
                                [zeros(5, n+10, n); zeros(n,5,n), M_RFQ, zeros(n,5,n); zeros(5, n+10, n)],...
                                     zeros(n+10, n+10, 5));
                M_RFQ = M_RFQ_B(4:end-3,4:end-3,4:end-3);
                n = n+4;

                if n ~= size(M_RFQ,1)
                    msg = 'Border update dimensions do not match.';
                    error(msg)     
                end

                Fit = cat(3, zeros(n,n,2),...
                                [zeros(2, n, n-4); zeros(n-4,2,n-4), Fit, zeros(n-4,2,n-4); zeros(2, n, n-4)],...
                                     zeros(n, n, 2));
                M_age = cat(3, zeros(n,n,2),...
                                [zeros(2, n, n-4); zeros(n-4,2,n-4), M_age, zeros(n-4,2,n-4); zeros(2, n, n-4)],...
                                     zeros(n, n, 2));
                new_pos = new_pos+2; 

            else
                
                2

                M_RFQ_B = cat(3, zeros(n+8,n+8,4),...
                                [zeros(4, n+8, n); zeros(n,4,n), M_RFQ, zeros(n,4,n); zeros(4, n+8, n)],...
                                     zeros(n+8, n+8, 4));
                M_RFQ = M_RFQ_B(4:end-3,4:end-3,4:end-3);
                n = n+2;

                if n ~= size(M_RFQ,1)
                    msg = 'Border update dimensions do not match.';
                    error(msg)     
                end

                Fit = cat(3, zeros(n,n,1),...
                                [zeros(1, n, n-2); zeros(n-2,1,n-2), Fit, zeros(n-2,1,n-2); zeros(1, n, n-2)],...
                                     zeros(n, n, 1));
                M_age = cat(3, zeros(n,n,1),...
                                [zeros(1, n, n-2); zeros(n-2,1,n-2), M_age, zeros(n-2,1,n-2); zeros(1, n, n-2)],...
                                     zeros(n, n, 1));
                new_pos = new_pos+1;

            end         
            
        else

            M_RFQ_B = cat(3, zeros(n+6,n+6,3),...
                            [zeros(3, n+6, n); zeros(n,3,n), M_RFQ, zeros(n,3,n); zeros(3, n+6, n)],...
                                 zeros(n+6, n+6, 3));
            M_RFQ = M_RFQ_B(3:end-2,3:end-2,3:end-2);
            n = n+2;

            if n ~= size(M_RFQ,1)
                msg = 'Border update dimensions do not match.';
                error(msg)     
            end

            Fit = cat(3, zeros(n,n,1),...
                            [zeros(1, n, n-2); zeros(n-2,1,n-2), Fit, zeros(n-2,1,n-2); zeros(1, n, n-2)],...
                                 zeros(n, n, 1));
            M_age = cat(3, zeros(n,n,1),...
                            [zeros(1, n, n-2); zeros(n-2,1,n-2), M_age, zeros(n-2,1,n-2); zeros(1, n, n-2)],...
                                 zeros(n, n, 1));
            new_pos = new_pos+1;

        end


    end


end
