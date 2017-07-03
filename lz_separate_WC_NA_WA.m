
%%%% 06?25?2017 
%%%% based on the computed wavelet coherence and NA/WA borders, separate
%%%% NA/WA scalegrams.

function [wc_NW, wc_WA] = lz_separate_WC_NA_WA(mouse_id, exp_session)

nCh = 30;
% load wc and labels
[folder_name, mouse_name, loadName_RS_ind, loadName_wc, loadName_WA_ind, loadName_RS_ind_pool, loadName_WA_ind_pool] ...
    = lz_build_folder_name(mouse_id, exp_session);

fprintf('\nLoading computed wavelet coherent for subject %d, ''%s'' ......', mouse_id, exp_session);
if ispc
    !hostname > hostname.txt
    if strcmp(textread('hostname.txt','%s'), 'Victor-THINK')
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_RS_ind]);
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_WA_ind]);
        load(['D:\Victor\Dropbox\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',loadName_wc]);
    else
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_RS_ind]);
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',loadName_WA_ind]);
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',loadName_wc]);
    end
elseif ismac
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind]);
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_WA_ind]);
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
elseif isunix
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind]);
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_WA_ind]);
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
end

if (mouse_id == 2 && strcmp(exp_session, 'am')) || (mouse_id == 3 && strcmp(exp_session, 'pm'))
    nTr = 15;
else
    nTr = 16;
end

% initiate concatenated wc associated WA and NW conditions
wc_NW = nan(109,0,nCh*(nCh-1)/2); 
wc_WA = nan(109,0,nCh*(nCh-1)/2);

for iTr = 1:nTr
    tic
    fprintf('\n\nStacking on trial %02d of subject %d, ''%s'' ...... \n', iTr, mouse_id, exp_session);
    % resting state
    if ~isempty(RS_eff.ind_start_cal{iTr})
        for iSeg_RS = 1: length(RS_eff.ind_start_cal{iTr})
            for iChan_pair = 1: nCh*(nCh-1)/2
                wc_NW(:,:,iChan_pair) = cat(2, wcoh(:, RS_eff.ind_start_cal{iTr}(iSeg_RS):RS_eff.ind_end_cal{iTr}(iSeg_RS),iChan_pair,iTr), wc_NW(:,:,iChan_pair));
            end
        end
    end

    % whisking active
    if ~isempty(WA_eff.ind_start_cal{iTr})
        for iSeg_WA = 1: length(WA_eff.ind_start_cal{iTr})
            for iChan_pair = 1: nCh*(nCh-1)/2
                wc_WA(:,:,iChan_pair) = cat(2, wcoh(:, WA_eff.ind_start_cal{iTr}(iSeg_WA):WA_eff.ind_end_cal{iTr}(iSeg_WA),iChan_pair,iTr), wc_WA(:,:,iChan_pair));
            end
        end
    end
    toc
end

saveName_RS = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_wc_RS', mouse_id, exp_session);
saveName_WA = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_wc_WA', mouse_id, exp_session);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',saveName_RS], 'wc_NW', 'f', 'coi', '-v7.3');
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',saveName_WA], 'wc_WA', 'f', 'coi', '-v7.3');
elseif ismac
    save(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',saveName_RS], 'wc_NW', 'f', 'coi', '-v7.3');
    save(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',saveName_WA], 'wc_WA', 'f', 'coi', '-v7.3');
elseif isunix
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',saveName_RS], 'wc_NW', 'f', 'coi', '-v7.3');
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',saveName_WA], 'wc_WA', 'f', 'coi', '-v7.3');
end