%%%% 06/26/2017
%%%% compute stationary mean and std for RS

function [wc_NW_mean_within_subj, wc_NW_std_within_subj] = lz_mean_std_WC_NA(mouse_id, exp_session)

% load wc and labels
[folder_name, mouse_name, loadName_RS_ind, loadName_wc, loadName_WA_ind, loadName_wc_RS, loadName_wc_WA] ...
    = lz_build_folder_name(mouse_id, exp_session);

fprintf('\nLoading computed RS wavelet coherent for subject %d, ''%s'' ......', mouse_id, exp_session);
if ispc
    !hostname > hostname.txt
    if strcmp(textread('hostname.txt','%s'), 'Victor-THINK')
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',loadName_wc_RS]);
    else
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',loadName_wc_RS]);
    end
elseif ismac
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc_RS]);
elseif isunix
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc_RS]);
end

wc_NW_mean_within_subj = mean(wc_NW, 2);
wc_NW_std_within_subj = std(wc_NW, 0, 2);

wc_WA_mean_within_subj = mean(wc_WA, 2);
wc_WA_std_within_subj = std(wc_WA, 0, 2);
