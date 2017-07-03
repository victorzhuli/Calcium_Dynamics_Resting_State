%%%% 07/02/2017
%%%% grab result from lz_separate_WC_NA_WA_ind.m
%%%% 06/26/2017
%%%% compute stationary mean and std for RS

function lz_mean_std_t_NW_WA_2(mouse_id, exp_session)

nCh = 30;
% load wc and labels
[folder_name, mouse_name, loadName_RS_ind, loadName_wc, loadName_WA_ind, loadName_wc_RS, loadName_wc_WA, loadName_RS_ind_pool, loadName_WA_ind_pool] ...
    = lz_build_folder_name(mouse_id, exp_session);

fprintf('\nLoading computed wavelet coherent for subject %d, ''%s'' ......\n', mouse_id, exp_session);
if ispc
    !hostname > hostname.txt
    if strcmp(textread('hostname.txt','%s'), 'Victor-THINK')
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_RS_ind_pool]);
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_WA_ind_pool]);
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',loadName_wc]);
    else
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_RS_ind_pool]);
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_WA_ind_pool]);
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',loadName_wc]);
    end
elseif ismac
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind_pool]);
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_WA_ind_pool]);
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
elseif isunix
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind_pool]);
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_WA_ind_pool]);
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
end

if (mouse_id == 2 && strcmp(exp_session, 'am')) || (mouse_id == 3 && strcmp(exp_session, 'pm'))
    nTr = 15;
else
    nTr = 16;
end
nFreq = size(wcoh,1);
nChPair = nCh*(nCh-1)/2;

%%%% initial
wc_NW_mean = nan(nFreq, nChPair);
wc_NW_std  = nan(nFreq, nChPair);
wc_NW_num  = nan(nFreq, nChPair);
wc_WA_mean = nan(nFreq, nChPair);
wc_WA_std  = nan(nFreq, nChPair);
wc_WA_num  = nan(nFreq, nChPair);

wc_t       = nan(nFreq, nChPair);

for iFreq = 1:size(wcoh,1)
    fprintf('\nCompute mean and std for subject %d, session %s, Freq bin %03d ...\n', mouse_id, exp_session, iFreq);

    for iCh_pair = 1 : nChPair
        wcoh_freq_chPair = squeeze(wcoh(iFreq,:,iCh_pair,:));
        wc_NW_pool_chPair_freq = [];
        wc_WA_pool_chPair_freq = [];
        
        for iTr = 1:nTr
            if ~isempty(wc_NW_ind{iTr})
                wc_NW_pool_chPair_freq = [wc_NW_pool_chPair_freq, wcoh_freq_chPair(wc_NW_ind{iTr}, iTr)'];
            end
            if ~isempty(wc_WA_ind{iTr})
                wc_WA_pool_chPair_freq = [wc_WA_pool_chPair_freq, wcoh_freq_chPair(wc_WA_ind{iTr}, iTr)'];
            end
        end
        
        wc_NW_mean(iFreq, iCh_pair) = mean(wc_NW_pool_chPair_freq);
        wc_NW_std(iFreq, iCh_pair)  = std(wc_NW_pool_chPair_freq);
        wc_NW_num(iFreq, iCh_pair)  = length(wc_NW_pool_chPair_freq);
        
        wc_WA_mean(iFreq, iCh_pair) = mean(wc_WA_pool_chPair_freq);
        wc_WA_std(iFreq, iCh_pair)  = std(wc_WA_pool_chPair_freq);
        wc_WA_num(iFreq, iCh_pair)  = length(wc_WA_pool_chPair_freq);
        
        wc_t(iFreq, iCh_pair) = lz_ttest2(wc_NW_pool_chPair_freq, wc_WA_pool_chPair_freq);
        
    end
end

saveName = sprintf('GC6f_emx_%02d_%s_wc_stat_NW_WA', mouse_id, exp_session);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoh_stat\',saveName], ...
        'wc_NW_mean','wc_NW_std','wc_NW_num','wc_WA_mean','wc_WA_std','wc_WA_num','wc_t', 'f');
elseif ismac
    save(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoh_stat/',saveName], ...
        'wc_NW_mean','wc_NW_std','wc_NW_num','wc_WA_mean','wc_WA_std','wc_WA_num','wc_t', 'f');
elseif isunix
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoh_stat/',saveName], ...
        'wc_NW_mean','wc_NW_std','wc_NW_num','wc_WA_mean','wc_WA_std','wc_WA_num','wc_t', 'f');
end













