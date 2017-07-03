%%%% 07/02/2017
%%%% developed from lz_separate_WC_NA_WA (06/25/2017) - pool index of wcoh
%%%% for NW and WA conditions, instead of pooling wcoh data.
%%%% based on the computed wavelet coherence and NA/WA borders, separate
%%%% NA/WA scalegrams.

function [wc_NW_ind, wc_WA_ind] = lz_separate_WC_NA_WA_ind(mouse_id, exp_session)

% load wc and labels
[folder_name, mouse_name, loadName_RS_ind, loadName_wc, loadName_WA_ind] = lz_build_folder_name(mouse_id, exp_session);

fprintf('\nLoading computed wavelet coherent for subject %d, ''%s'' ......', mouse_id, exp_session);
if ispc
    !hostname > hostname.txt
    if strcmp(textread('hostname.txt','%s'), 'Victor-THINK')
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_RS_ind]);
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_WA_ind]);
    else
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_RS_ind]);
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_WA_ind]);
    end
elseif ismac
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind]);
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_WA_ind]);
elseif isunix
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind]);
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_WA_ind]);
end

if (mouse_id == 2 && strcmp(exp_session, 'am')) || (mouse_id == 3 && strcmp(exp_session, 'pm'))
    nTr = 15;
else
    nTr = 16;
end

% initiate concatenated wc associated WA and NW conditions
wc_NW_ind = cell(nTr,1); 
wc_WA_ind = cell(nTr,1);

for iTr = 1:nTr
    tic
    fprintf('\n\nStacking index for trial %02d of subject %d, ''%s'' ...... \n', iTr, mouse_id, exp_session);
    % resting state
    if ~isempty(RS_eff.ind_start_cal{iTr})
        for iSeg_RS = 1: length(RS_eff.ind_start_cal{iTr})
            wc_NW_ind{iTr} = [wc_NW_ind{iTr}, RS_eff.ind_start_cal{iTr}(iSeg_RS):RS_eff.ind_end_cal{iTr}(iSeg_RS)];
        end
    end

    % whisking active
    if ~isempty(WA_eff.ind_start_cal{iTr})
        for iSeg_WA = 1: length(WA_eff.ind_start_cal{iTr})
            wc_WA_ind{iTr} = [wc_WA_ind{iTr}, WA_eff.ind_start_cal{iTr}(iSeg_WA):WA_eff.ind_end_cal{iTr}(iSeg_WA)];
        end
    end
    toc
end

saveName_RS = sprintf('GC6f_emx_%02d_%2s_spont_wc_RS_ind', mouse_id, exp_session);
saveName_WA = sprintf('GC6f_emx_%02d_%2s_spont_wc_WA_ind', mouse_id, exp_session);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',saveName_RS], 'wc_NW_ind');
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',saveName_WA], 'wc_WA_ind');
elseif ismac
    save(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',saveName_RS], 'wc_NW_ind');
    save(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',saveName_WA], 'wc_WA_ind');
elseif isunix
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',saveName_RS], 'wc_NW_ind');
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',saveName_WA], 'wc_WA_ind');
end

