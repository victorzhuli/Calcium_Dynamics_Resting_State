%%%% 07/04/2017
%%%% based on the separated NA/WA scalegrams, conduct permutation test for
%%%% NW vs. WA across subjects

function [wc_permut_pval, wc_permut_t_orig, wc_permut_crit_t, wc_permut_est_alpha] = lz_permutation_NW_WA_2(pVal)

%%%% stack the mean WC matrix across subjects for NW and WA respectively
wc_NW_mn_stack = [];
wc_WA_mn_stack = [];

exp_session_array = {'am', 'pm'}; 
for mouse_id = 1:6
    for iexp_session = 1:2
        exp_session = exp_session_array{iexp_session};  
        
        % exclude mouse 3 pm because of sparse WA
        if mouse_id == 3 && strcmp(exp_session, 'pm')
            continue
        end
        
        % load wc and labels
        [folder_name, mouse_name, loadName_RS_ind, loadName_wc, loadName_WA_ind, ...
            loadName_wc_RS, loadName_wc_WA, loadName_RS_ind_pool, loadName_WA_ind_pool, ...
            loadName_wc_stat, loadName_wc_separate] ...
            = lz_build_folder_name(mouse_id, exp_session);
          
        if ispc
            !hostname > hostname.txt
            if strcmp(textread('hostname.txt','%s'), 'Victor-THINK')
                load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoh_stat\',loadName_wc_stat]);
            else
                load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoh_stat\',loadName_wc_stat]);
            end
        elseif ismac
            load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoh_stat/',loadName_wc_stat]);
        elseif isunix
            load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoh_stat/',loadName_wc_stat]);
        end
                
        wc_NW_mn_stack = cat(1, wc_NW_mn_stack, wc_NW_mean);
        wc_WA_mn_stack = cat(1, wc_WA_mn_stack, wc_WA_mean);
    end
end

fprintf('\nPermutation test for wavelet coherence across subject for NW vs. WA ......\n');
[nFreq, nChPair] = size(wc_NW_mean);

%%%% initiate
wc_permut_pval = [];
wc_permut_t_orig = [];
wc_permut_crit_t = [];
wc_permut_est_alpha = [];
for iChPair = 1:nChPair
    data_NW = reshape( wc_NW_mn_stack(:,iChPair), nFreq, 11)';
    data_WA = reshape( wc_WA_mn_stack(:,iChPair), nFreq, 11)';
    dif = data_NW - data_WA;
    
    [wc_permut_pval(:, iChPair), wc_permut_t_orig(:, iChPair), wc_permut_crit_t(:, iChPair), ...
        wc_permut_est_alpha(:, iChPair)] ...
        = mult_comp_perm_t1(dif,50000, 0, pVal,0,0);
end

    
    




