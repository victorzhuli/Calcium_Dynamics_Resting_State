%%%% Study the dynamics of resting state, prepared for SPIE Conference
%%%% 2017 in San Diego.
%%
mouse_id = 1;
exp_session = 'am';
% Load HbO and HbR signals
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
    cd('/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/Ca.mat']); % Calcium signals
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/anglekeeper.mat']); % whisker signals
elseif isunix
    cd('/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
    load(['/home/lz206/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/Ca.mat']); % Calcium signals
    load(['/home/lz206/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/anglekeeper.mat']); % whisker signals
end
%--------------------------------------------------------------------------

% change cell to matrix: channel X time X trial
Cal = reshape(cell2mat(Ca.Ch0), 30, 2047, size(Ca.Ch0,2));
if (mouse_id == 03 && strcmp(exp_session, 'pm'))
    Cal = Cal(:, 1:2000, [1:11,13:16]);
else
    Cal = Cal(:,1:2000,:);
end
%%
RS_eff = lz_labeling_nonWhisking(mouse_id, exp_session);
%--------------------------------------------------------------------------
%% Check point
ch2plot = 8; tr2plot = 16;
figure(1);clf;
subplot(411);plot(Cal(ch2plot,:,1),'b'); title(['Channel ',num2str(ch2plot),', Trial ',num2str(tr2plot)]);
hold on; 
for iRS = 1: length(RS_eff.ind_start_cal{tr2plot})
    h1 = vline(RS_eff.ind_start_cal{tr2plot}(iRS),'k'); set(h1, 'linew', 2);
    h2 = vline(RS_eff.ind_end_cal{tr2plot}(iRS),'r'); set(h2, 'linew', 2);
end
subplot(412);plot(anglekeeper(tr2plot,:)); xlim([0 length(anglekeeper(tr2plot,:))]);
subplot(413);stem(whisker_conv_win(tr2plot,:),'r');  ylim([-.5, 1.5]);
subplot(414);stem(diff_whisker_flag(tr2plot,:),'r'); ylim([-1.5, 1.5]);
%--------------------------------------------------------------------------





