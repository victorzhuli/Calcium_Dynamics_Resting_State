%%%% 06/22/2017
%%%% compute wavelet coherence
function [wcoh,f,coi] = lz_calcium_wavelet_coherence(mouse_id, exp_session)

srate_cal = 100;

% load calcium signal
[folder_name, mouse_name] = lz_build_folder_name(mouse_id, exp_session);
if ispc
elseif ismac
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/',folder_name,'/',mouse_name,'/Ca.mat']);
elseif isunix
    load(['/home/lz206/Dropbox/GCaMP6f spont and tone reward/',folder_name,'/',mouse_name,'/Ca.mat']);
end

% change cell to matrix: channel X time X trial
Cal = reshape(cell2mat(Ca.Ch0), 30, 2047, size(Ca.Ch0,2));
if (mouse_id == 03 && strcmp(exp_session, 'pm'))
    Cal = Cal(:, 1:2000, [1:11,13:16]);
else
    Cal = Cal(:,1:2000,:);
end
nTr = size(Cal,3);

% initiate
wcoh = nan(109,2000,30*29/2,nTr); % 109 is the number of frequency bins
for iTr = 1:nTr
    iChan_pair = 1;
    for iCh = 1:30
        for jCh = iCh+1 : 30
            [wcoh(:,:,iChan_pair,iTr),~,f,coi] = wcoherence(Cal(iCh,:,iTr),Cal(jCh,:,iTr),srate_cal);
            iChan_pair = iChan_pair + 1;
        end
    end
end

saveName = sprintf('GC6f_emx_%02d_%2s_spont_resting_state_wcoherence', mouse_id, exp_session);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State\wcoherence\',saveName], 'wcoh', 'f', 'coi', '-v7.3');
elseif ismac
    
elseif isunix
    save(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',saveName], 'wcoh', 'f', 'coi', '-v7.3');
end