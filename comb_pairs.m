% Use this function to create a parameter list containing all combinations
% of pre-defined parmeter pairs. This function differs from comb_para.m by
% keeping the pre-defined pair(s) together in each combination,
% > if you want to keep a parameter constant, add one value
% > if you want to add a parameter pair:
%      1. Comment the definition of the parameter, which are covered by
%         your pairs
%      2. Define number k of parameters each pair contains (e.g. 2 - g_RR and g_RQ are a pair)
%         > this function offers only combination of pairs with the same
%           number of parameters
%      3. Create a n x k array of parameter pairs, where n is the number of
%         different parameter pairs:
%         > example: [g_RR_value_1, g_RQ_value_1;
%                     g_RR_value_2, g_RQ_value_2;
%                     ...
%                     g_RR_value_n, g_RQ_value_n]
%         > you can also load lists from previous evaluations (see example
%           in code)
%      4. Include all parameter-pair-arrays (example: g_Res, g_For) in para_set

clear all;

g_VV = 0.;
g_RV = 0.;
g_FV = 0.;
g_QV = 0.;

g_VR = 0.8;
%g_RR = [0.9];
g_FR = 0;
g_QR = 0.;

g_VF = 0.;
g_RF = 0;
%g_FF = -0.85;
g_QF = 0.3;

g_VQ = 0.;
%g_RQ = [-0.6];
%g_FQ = [0.4];
g_QQ = 0.;

%% example

length_pair = 2;

% to be substituted with your lists
load('paralists_eval_Res.mat');
load('paralists_eval_For.mat');

g_Res = [g_RR, g_RQ];
g_For = [g_FF, g_FQ];

%%

para_set = {g_Res ; g_For};

pair_list = comb_parameter(para_set, length_pair);

para_list = [repmat([g_VV, g_RV, g_FV, g_QV, g_VR],length(pair_list),1),...
    pair_list(:,1), repmat([g_FR, g_QR, g_VF, g_RF],length(pair_list),1),...
    pair_list(:,3), repmat([g_QF, g_VQ],length(pair_list),1),...
    pair_list(:,2), pair_list(:,4), repmat([g_QQ],length(pair_list),1)];

save('initial_para_list', 'para_list');
% optimally use significant file name, e.g.:
% save('combSim_gRRgRQ_gFFgFQ', 'para_list')

function comb_set = comb_parameter(para_set, length_pair)

    num_comb = prod(cellfun('length',para_set));
    num_pos = length(para_set)*length_pair;

    comb_set = zeros(num_comb, num_pos);
    step = num_comb;

    for i = 1:(num_pos/length_pair)

        num_para = length(para_set{i});

        step = step/length(para_set{i});

        pattern = zeros(step*num_para,length_pair);

        for j = 0:num_para-1

            pattern(j*step+1:(j+1)*step,:) = repmat(para_set{i}(j+1,:),step,1);

        end

        len_pat = length(pattern);  
        rep_pat = num_comb/len_pat;

        for j = 0:rep_pat-1

            comb_set(j*len_pat+1:(j+1)*len_pat,(i-1)*2+1:i*2)=pattern;

        end

    end

end