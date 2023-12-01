function [M_init, g_bar, stencil_ind, neighb, def_neighbourX,...
          def_neighbourY, def_neighbourZ] = create_setup(n, type, G, NB)
    
    switch NB
        
        case 'vN_1'
            
            def_neighbourX = [-1 0 1];
            def_neighbourY = [-1 0 1];
            def_neighbourZ = [-1 0 1];

            neighb = [0 0 1;  0 1 0; 1 0 0; 0 0 -1; 0 -1 0; -1 0 0];

            stencil_legit_neigh = ones(3,3,3);
            stencil_legit_neigh(2,2,1:3) = 0;
            stencil_legit_neigh(2,1:3,2) = 0;
            stencil_legit_neigh(1:3,2,2) = 0;
            stencil_legit_neigh(2,2,2) = 1;
            stencil_ind = find(stencil_legit_neigh);
            
        case 'M_1'
            
            def_neighbourX = [-1 0 1];
            def_neighbourY = [-1 0 1];
            def_neighbourZ = [-1 0 1];
            
            neighb = comb_parameter({def_neighbourX; def_neighbourY; def_neighbourZ});

            stencil_legit_neigh = zeros(3,3,3);
            stencil_legit_neigh(2,2,2) = 1;
            stencil_ind = find(stencil_legit_neigh);         
            
        case 'M_2'
            
            def_neighbourX = [-2 -1 0 1 2];
            def_neighbourY = [-2 -1 0 1 2];
            def_neighbourZ = [-2 -1 0 1 2];

            neighb = comb_parameter({def_neighbourX; def_neighbourY; def_neighbourZ});

            stencil_legit_neigh = zeros(5,5,5);
            stencil_legit_neigh(3,3,3) = 1;
            stencil_ind = find(stencil_legit_neigh);
            
    end
    
    g_bar = [flip(G(1,:)');flip(G(2,:)');flip(G(3,:)');flip(G(4,:)')];
    
    M_init = create_initial_config(type, n);
    
    %% Functions

    function M = create_initial_config(type, n)

        switch type
                
            case 'quiet_resorp_2'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2,14:16) = ones(3,1,3);
                
            case 'quiet_resorp_3'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2,14:16) = zeros(3,1,3);
                M(14:16,3,14:16) = ones(3,1,3);
                
            case 'quiet_resorp_4'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2:3,14:16) = zeros(3,2,3);
                M(14:16,4,14:16) = ones(3,1,3);
                
            case 'quiet_resorp_5'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2:4,14:16) = zeros(3,3,3);
                M(14:16,5,14:16) = ones(3,1,3);
                
            case 'quiet_form_2'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2,14:16) = 2.*ones(3,1,3);
                
            case 'quiet_form_3'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2,14:16) = zeros(3,1,3);
                M(14:16,3,14:16) = 2.*ones(3,1,3);
                
            case 'quiet_form_4'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2:3,14:16) = zeros(3,2,3);
                M(14:16,4,14:16) = 2.*ones(3,1,3);
                
            case 'quiet_form_5'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,2:4,14:16) = zeros(3,3,3);
                M(14:16,5,14:16) = 2.*ones(3,1,3);
                        
            case 'BMU'
                
                M = zeros(n,n,n);
                M(2:n-1, 2:n-1, 2:n-1) = 3;
                
                M(14:16,5,14:16) = ones(3,1,3);
                M(14:16,4,14:16) = 2.*ones(3,1,3);
                M(14:16,2:3,14:16) = zeros(3,2,3);
                               
        end
    end

    function comb_set = comb_parameter(para_set)

        num_comb = prod(cellfun('length',para_set));
        num_pos = length(para_set);

        comb_set = zeros(num_comb, num_pos);
        step = num_comb;

        for i = 1:num_pos

            num_para = length(para_set{i});

            step = step/length(para_set{i});

            pattern = zeros(step*num_para,1);

            for j = 0:num_para-1

                pattern(j*step+1:(j+1)*step) = para_set{i}(j+1);

            end

            len_pat = length(pattern);  
            rep_pat = num_comb/len_pat;

            for j = 0:rep_pat-1

                comb_set(j*len_pat+1:(j+1)*len_pat,i)=pattern;

            end

            comb_set(find(sum(abs(comb_set),2)==0),:)=[];

        end

    end

end