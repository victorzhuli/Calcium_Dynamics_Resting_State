%%%% 06/23/2017
%%%% %% combine 'am' and 'pm' for both NA label and WC pattern
function combine_RS_ind_wc

for mouse_id = 1:6
    fprintf('Combining mouse %d', mouse_id);
    
    % load exp_session specified files
    exp_session = 'am';
    [folder_name, mouse_name, loadName_RS_ind, loadName_wc] = lz_build_folder_name(mouse_id, exp_session);
    if ispc
        RS_ind_am = load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',loadName_RS_ind]);
        wc_am = load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
    elseif isunix
        RS_ind_am = load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',loadName_RS_ind]);
        wc_am = load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
    end
    
    exp_session = 'pm';
    [folder_name, mouse_name, loadName_RS_ind, loadName_wc] = lz_build_folder_name(mouse_id, exp_session);
    if ispc
        RS_ind_pm = load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',loadName_RS_ind]);
        wc_pm = load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
    elseif isunix
        RS_ind_pm = load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',loadName_RS_ind]);
        wc_pm = load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
    end
    
    % combine files across exp_sessions
    RS_ind_all.ind_start_cal = horzcat(RS_ind_am.ind_start_cal, RS_ind_pm.ind_start_cal); 
    RS_ind_all.ind_end_cal = horzcat(RS_ind_am.ind_end_cal, RS_ind_pm.ind_end_cal);
    
    wc_all = cat(4, wc_am.wcoh, wc_pm.wcoh);
    f = wc_am.f;
    coi = wc_am.f;
    
    saveName_RS_ind = sprintf('GC6f_emx_%02d_spont_resting_state_index', mouse_id);
    saveName_wc = sprintf('GC6f_emx_%02d_spont_resting_state_wcoherence', mouse_id);
    if ispc
        save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_index\',saveName_RS_ind], 'RS_ind_all');
        save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',saveName_wc], 'wc_all', 'f', 'coi', '-v7.3');
    elseif ismac
        
    elseif isunix
        save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',saveName_RS_ind], 'RS_ind_all');
        save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',saveName_wc], 'wc_all', 'f', 'coi', '-v7.3');
    end
end