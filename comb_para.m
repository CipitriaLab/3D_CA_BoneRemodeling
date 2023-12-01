% Use this function to create a prameter list containing all combinations
% of pre-defined parmeter values.
% > if you want to keep a parameter constant, add one value
% > if you want to vary a parameter, add a list of values

clear all;

g_EE = [0.];
g_RE = [0.];
g_FE = [0.];
g_QE = [0.];

g_ER = [0.5];
g_RR = [-1.:0.05:1.];
g_FR = [0.];
g_QR = [0.];

g_EF = [0.];
g_RF = [0.];
g_FF = [0.];
g_QF = [0.];

g_EQ = [0.];
g_RQ = [-1.:0.05:1.];
g_FQ = [0.];
g_QQ = [0.];

para_set = {g_EE ; g_RE ; g_FE ; g_QE ; g_ER ; g_RR ; g_FR ; g_QR ;...
                g_EF ; g_RF ; g_FF ; g_QF ; g_EQ ; g_RQ ; g_FQ ; g_QQ};

para_list = comb_parameter(para_set);

save('initial_para_list', 'para_list');
% optimally use significant file name, e.g.:
% save('singlePatchRes_gER_05', 'para_list')

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



    end

end