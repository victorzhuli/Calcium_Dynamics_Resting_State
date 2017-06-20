%%%% Study the dynamics of resting state, prepared for SPIE Conference
%%%% 2017 in San Diego.
%%
mouse_id = 1;
exp_session = 'am';
%% Load HbO and HbR signals
mouse_name = sprintf('GC6f_emx_%02d',mouse_id);

if (mouse_id == 01 || mouse_id == 02 || mouse_id == 03) &&  strcmp(exp_session, 'am')
    folder_name = '150421am GC6-emx 1-3 spont';
elseif (mouse_id == 01 || mouse_id == 02 || mouse_id == 03) &&  strcmp(exp_session, 'pm')
    folder_name = '150421pm GC6-emx 1-3 spont';
elseif (mouse_id == 04 || mouse_id == 05 || mouse_id == 06) &&  strcmp(exp_session, 'am')
    folder_name = '150609am GC6-emx 4-6 spont';
else
    folder_name = '150609pm GC6-emx 4-6 spont';
end

%--------------------------------------------------------------------------
if ispc
    
elseif ismac
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/Ca.mat']);
elseif isunix
    
end
%--------------------------------------------------------------------------

