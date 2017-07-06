%%%% 06/22/2017
%%%% construct folder name and mouse name for loading spontaneous mouse data
function [folder_name, mouse_name, loadName_RS_ind, loadName_wc, loadName_WA_ind, ...
    loadName_wc_RS, loadName_wc_WA, loadName_RS_ind_pool, loadName_WA_ind_pool, ...
    loadName_wc_stat, loadName_wc_separate] ...
    = lz_build_folder_name(mouse_id, exp_session)

mouse_name = sprintf('GC6f_emx_%02d',mouse_id);
loadName_RS_ind_pool = sprintf('GC6f_emx_%02d_%2s_spont_wc_RS_ind', mouse_id, exp_session);
loadName_RS_ind = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_index', mouse_id, exp_session);
loadName_wc = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_wcoherence', mouse_id, exp_session);
loadName_WA_ind_pool = sprintf('GC6f_emx_%02d_%2s_spont_wc_WA_ind', mouse_id, exp_session);
loadName_WA_ind = sprintf('GC6f_emx_%02d_%2s_spont_whisking_active_index', mouse_id, exp_session);
loadName_wc_RS = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_wc_RS', mouse_id, exp_session);
loadName_wc_WA = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_wc_WA', mouse_id, exp_session);
loadName_wc_stat = sprintf('GC6f_emx_%02d_%s_wc_stat_NW_WA', mouse_id, exp_session);
loadName_wc_separate = sprintf('GC6f_emx_%02d_%s_wc_separate_NW_WA', mouse_id, exp_session);

if (mouse_id == 01 || mouse_id == 02 || mouse_id == 03) &&  strcmp(exp_session, 'am')
    folder_name = '150421am GC6-emx 1-3 spont';
elseif (mouse_id == 01 || mouse_id == 02 || mouse_id == 03) &&  strcmp(exp_session, 'pm')
    folder_name = '150421pm GC6-emx 1-3 spont';
elseif (mouse_id == 04 || mouse_id == 05 || mouse_id == 06) &&  strcmp(exp_session, 'am')
    folder_name = '150609am GC6-emx 4-6 spont';
else
    folder_name = '150609pm GC6-emx 4-6 spont';
end