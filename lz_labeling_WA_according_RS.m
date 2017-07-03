%%%% 06/24/2017
%%%% look for labels for WA according to NA labels
function WA_eff = lz_labeling_WA_according_RS(mouse_id, exp_session)

% load RS_eff
[folder_name, mouse_name, loadName_RS_ind, loadName_wc] = lz_build_folder_name(mouse_id, exp_session);
if ispc
elseif ismac
    load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind]);
elseif isunix
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',loadName_RS_ind]);
end
        
len = 2000;
nTr = length(RS_eff.ind_start_cal);

% initiate
WA_eff.ind_start_cal = cell(1, nTr);
WA_eff.ind_end_cal = cell(1, nTr);

for iTr = 1: nTr
    nSeg_RS = length(RS_eff.ind_start_cal{iTr});
    
    % if whisking active for all durating (no resting state exist) for this trial
    if nSeg_RS == 0
        WA_eff.ind_start_cal{iTr} = 1;
        WA_eff.ind_end_cal{iTr} = len;
        
    % if there is resting state exist    
    else
        for iSeg_RS = 1:nSeg_RS
            
            % first RS seg
            if iSeg_RS == 1
                if RS_eff.ind_start_cal{iTr}(iSeg_RS) ~= 1
                    WA_eff.ind_start_cal{iTr} = 1;
                    WA_eff.ind_end_cal{iTr}   = RS_eff.ind_start_cal{iTr}(iSeg_RS) - 1;
                end
                
                % following RS seg
            else
                WA_eff.ind_start_cal{iTr} = [WA_eff.ind_start_cal{iTr}, RS_eff.ind_end_cal{iTr}(iSeg_RS-1) + 1];
                WA_eff.ind_end_cal{iTr}   = [WA_eff.ind_end_cal{iTr}, RS_eff.ind_start_cal{iTr}(iSeg_RS) - 1];
            end
            
        end
        
        % check the end of each trial
        if RS_eff.ind_end_cal{iTr}(nSeg_RS) ~= len
            WA_eff.ind_start_cal{iTr} = [WA_eff.ind_start_cal{iTr}, RS_eff.ind_end_cal{iTr}(nSeg_RS)+1];
            WA_eff.ind_end_cal{iTr}   = [WA_eff.ind_end_cal{iTr}, len];
        end
    end
end
     
% save results
saveName = sprintf('GC6f_emx_%02d_%2s_spont_whisking_active_index', mouse_id, exp_session);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_WA_index\',saveName], 'WA_eff');
elseif ismac
    save(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',saveName], 'WA_eff');
elseif isunix
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_WA_index/',saveName], 'WA_eff');
end
    