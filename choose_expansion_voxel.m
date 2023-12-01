function [min_val, min_x, min_y, min_z, Fit] = ...
    choose_expansion_voxel(n, t, neighb, M_RFQ_B, Fit, g_bar, ...
    update, NB, update_x, update_y, update_z)

    %% Calculate expansion rate
    
    if update
        
        % 1st step
        if t == 1

            for k = 1:length(neighb)
                
                if strcmp(NB, 'M_2')
                    
                    step = neighb(k,:);
                    inter = 4 .* M_RFQ_B(4:n+3,4:n+3,4:n+3) ...
                               - M_RFQ_B(4+step(1):n+3+step(1), 4+step(2):n+3+step(2), 4+step(3):n+3+step(3))...
                               + 4;
                    Fit = Fit + g_bar(inter);
                    
                else

                    step = neighb(k,:);
                    inter = 4 .* M_RFQ_B(3:n+2,3:n+2,3:n+2) ...
                               - M_RFQ_B(3+step(1):n+2+step(1), 3+step(2):n+2+step(2), 3+step(3):n+2+step(3))...
                               + 4;
                    Fit = Fit + g_bar(inter);
                    
                end

            end

        % 2nd-nth step
        else

            for k = 1:length(neighb)
                
                if strcmp(NB, 'M_2')
                    
                    step = neighb(k,:);
                    inter = 4 .* M_RFQ_B(update_x+3,update_y+3,update_z+3)...
                               - M_RFQ_B(update_x+3+step(1), update_y+3+step(2), update_z+3+step(3))...
                               + 4;
                    Fit(update_x,update_y,update_z) = Fit(update_x,update_y,update_z) + g_bar(inter);

                else

                    step = neighb(k,:);
                    inter = 4 .* M_RFQ_B(update_x+2,update_y+2,update_z+2)...
                               - M_RFQ_B(update_x+2+step(1), update_y+2+step(2), update_z+2+step(3))...
                               + 4;
                    Fit(update_x,update_y,update_z) = Fit(update_x,update_y,update_z) + g_bar(inter);

                end
                
            end

        end
        
    end

    %% Calculate expansion time and choose voxel
        
    % Exponential distribution
    M_exp = exprnd(1./Fit);

    % Expansion

    [min_val, min_idx] = min(M_exp(:));

    [min_x, min_y, min_z] = ind2sub(size(M_exp),min_idx);

end