%%%% Study the dynamics of resting state, prepared for SPIE Conference
%%%% 2017 in San Diego.
if ispc
    cd('C:\Users\Li_Lab537\Dropbox\projects\calcium\Calcium_Dynamics_Resting_State')
elseif ismac
    cd('/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
elseif isunix
    cd('/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
end
%% calculate the duration of boarder of non-whisking period (Resting-state) 
%%
exp_session_array = {'am', 'pm'};

for mouse_id = 1:6
    for iexp_session = 1:2
        exp_session = exp_session_array{iexp_session};        
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
        % call function
        RS_eff = lz_labeling_nonWhisking(mouse_id, exp_session);
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------      
    end
end
%% compute continuous wavelet coherence for each trial
%%
for mouse_id = 1:6
    for iexp_session = 1:2
        exp_session = exp_session_array{iexp_session};       
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
        % call function
        fprintf('\n\nComputing wavelet coherence on mouse %d, ''%s'' ......\n', mouse_id, exp_session);
        [wcoh,f,coi] = lz_calcium_wavelet_coherence(mouse_id, exp_session);              
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
    end
end

%% combine 'am' and 'pm' for both NA label and WC pattern
%%
%%%% unfortunately: out of memory!!
lz_combine_RS_ind_wc            

%% 
%% compute stationary mean and std
%%
exp_session_array = {'am', 'pm'}; nCh = 30;

% for mouse_id = 1%:6
%     for exp_session_array = 1%:2
        mouse_id = 1; iexp_session = 1;
        exp_session = exp_session_array{iexp_session};  
        
        % load wc and labels
        [folder_name, mouse_name, loadName_RS_ind, loadName_wc] = lz_build_folder_name(mouse_id, exp_session);
        if ispc
        elseif ismac
            load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',loadName_RS_ind]);
            load(['/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
        elseif isunix
            load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/RS_index/',loadName_RS_ind]);
            load(['/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State/wcoherence/',loadName_wc]);
        end
        
        % initiate concatenated wc associated WA and NW conditions
        wc_WA = [];  wc_NW = [];
        if (mouse_id == 2 && strcmp(exp_session, 'am')) || (mouse_id == 3 && strcmp(exp_session, 'pm'))
            nTr = 15;
        else
            nTr = 16;
        end
        
        tic
        for iTr = 1:nTr
            for iSeg = 1: length(RS_eff.ind_start_cal{iTr})
                for iChan_pair = 1: nCh*(nCh-1)/2
                    wc_WA = cat(2, wcoh(:, RS_eff.ind_start_cal{iSeg}:RS_eff.ind_end_cal{iSeg},iChan_pair), wc_WA);
%                     wc_NA =
                end
            end
        end
        toc
        
        
        
        % initiate results
        wc_WA_mean_within = nan(109,nCh*(nCh-1)/2);
%         wc_NW_mean_within = nan(109,nCh*(nCh-1)/2);
        
        
    end
end
        



%%=========================================================================
%% Check point for wavelet coherence --------------------------------------
%% plot two raw calcium signals and their wavelet coherence
mouse_id = 1; exp_session = 'am'; iTr = 5; iCh = 1; jCh = iCh+18;
lz_checkpoint_raw_calcium_pair_wc(mouse_id, exp_session, iTr, iCh, jCh);

% Check point for non-whisking labeling ----------------------------------
lz_checkpoint_nw_labeling(mouse_id, exp_session, iCh, iTr);





