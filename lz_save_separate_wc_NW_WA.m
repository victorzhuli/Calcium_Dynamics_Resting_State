%%%% 07/02/2017
%%%% grab result from lz_separate_WC_NA_WA_ind.m
%%%% save NA and WA wavelet coherence separately

function lz_save_separate_wc_NW_WA(mouse_id, exp_session)

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

%%%% length of NW and WA respectively?
len_NW_tmp = 0;
len_WA_tmp = 0;
for iTr = 1:nTr
    len_NW_tmp = len_NW_tmp + length(wc_NW_ind{iTr});
    len_WA_tmp = len_WA_tmp + length(wc_WA_ind{iTr});
end
    
%%%% initial
wc_NW = nan(nFreq, len_NW_tmp, nChPair);
wc_WA = nan(nFreq, len_WA_tmp, nChPair);

for iChPair = 1: nChPair
    for iFreq = 1:size(wcoh,1)

        fprintf('\nSeparating for subject %d, session %s, channel pair %03d, Freq bin %03d ...\n', mouse_id, exp_session, iChPair, iFreq);

        wcoh_freq_chPair = squeeze(wcoh(iFreq,:,iChPair,:));
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
        wc_NW(iFreq, :, iChPair) = wc_NW_pool_chPair_freq;
        wc_WA(iFreq, :, iChPair) = wc_WA_pool_chPair_freq;
        
    end
end

saveName = sprintf('GC6f_emx_%02d_%s_wc_separate_NW_WA', mouse_id, exp_session);
fprintf('\nSaving separated wavelet coherent for subject %d, ''%s'' ......\n', mouse_id, exp_session);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoh_separate\',saveName], ...
        'wc_NW','wc_WA', '-v7.3');
elseif ismac
    save(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoh_separate/',saveName], ...
        'wc_NW','wc_WA', '-v7.3');
elseif isunix
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoh_separate/',saveName], ...
        'wc_NW','wc_WA', '-v7.3');
end













