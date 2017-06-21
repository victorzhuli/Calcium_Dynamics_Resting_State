%%%% Study the dynamics of resting state, prepared for SPIE Conference
%%%% 2017 in San Diego.
%%
exp_session_array = {'am', 'pm'};

for mouse_id = 1:6
    for iexp_session = 1:2
       
        exp_session = exp_session_array{iexp_session};
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
        
        if ispc
            
        elseif ismac
            cd('/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
            load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/', ...
                folder_name,'/',mouse_name,'/Ca.mat']); % Calcium signals
        elseif isunix
            cd('/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
            load(['/home/lz206/Dropbox/GCaMP6f spont and tone reward/', ...
                folder_name,'/',mouse_name,'/Ca.mat']); % Calcium signals
        end
        
        % change cell to matrix: channel X time X trial
        Cal = reshape(cell2mat(Ca.Ch0), 30, 2047, size(Ca.Ch0,2));
        if (mouse_id == 03 && strcmp(exp_session, 'pm'))
            Cal = Cal(:, 1:2000, [1:11,13:16]);
        else
            Cal = Cal(:,1:2000,:);
        end
        %--------------------------------------------------------------------------
        % call function
        RS_eff = lz_labeling_nonWhisking(mouse_id, exp_session);
        %--------------------------------------------------------------------------
        
        saveName = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_index', mouse_id, exp_session);
        
        if ispc
            save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_index\',saveName], 'RS_eff');
        elseif ismac
            
        elseif isunix
            save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',saveName], 'RS_eff');
        end
        
    end
end


%--------------------------------------------------------------------------
%% Check point
mouse_id = 1; exp_session = 'am'; ch2plot = 8; tr2plot = 6;
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

loadName = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_index', mouse_id, exp_session);
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\RS_index\',loadName]);
elseif ismac
    
elseif isunix
    load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',loadName]);
    load(['/home/lz206/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/anglekeeper.mat']); % whisker signals
    load(['/home/lz206/Dropbox/GCaMP6f spont and tone reward/', ...
                folder_name,'/',mouse_name,'/Ca.mat']); % Calcium signals
end

figure(1);clf;
subplot(211);plot(linspace(0,20,2000), Cal(ch2plot,:,tr2plot),'k'); 
title(['Calcium Signal, Channel ',num2str(ch2plot),', Trial ',num2str(tr2plot)]);
hold on; 
for iRS = 1: length(RS_eff.ind_start_cal{tr2plot})
    h1 = vline(RS_eff.ind_start_cal{tr2plot}(iRS)/100,'b'); set(h1, 'linew', 2);
    h2 = vline(RS_eff.ind_end_cal{tr2plot}(iRS)/100,'r'); set(h2, 'linew', 2);
end
subplot(212);
plot(linspace(0,20,length(anglekeeper(tr2plot,:))), anglekeeper(tr2plot,:), 'k'); 
title(['Whisker Activitly, Trial ',num2str(tr2plot)])
%--------------------------------------------------------------------------





