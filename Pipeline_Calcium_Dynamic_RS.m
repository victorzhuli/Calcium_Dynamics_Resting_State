%%%% Study the dynamics of resting state, prepared for SPIE Conference
%%%% 2017 in San Diego.
if ispc
    
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

%% compute stationary mean and std
%%
exp_session_array = {'am', 'pm'};

for mouse_id = 1:6
    for exp_session_array = 1:2
        exp_session = exp_session_array{iexp_session};  
        
        



%%=========================================================================
%% Check point for wavelet coherence --------------------------------------
%% plot two raw calcium signals and their wavelet coherence
mouse_id = 1; exp_session = 'am'; iTr = 5; iCh = 1; jCh = iCh+18;
lz_checkpoint_raw_calcium_pair_wc(mouse_id, exp_session, iTr, iCh, jCh);

% Check point for non-whisking labeling ----------------------------------
lz_checkpoint_nw_labeling(mouse_id, exp_session, iCh, iTr);





