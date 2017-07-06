%%%% Study the dynamics of resting state, prepared for SPIE Conference
%%%% 2017 in San Diego.
if ispc
    
elseif ismac
    cd('/Users/lizhu/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
elseif isunix
    cd('/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
end
%% calculate the boarder of non-whisking period (Resting-state) 
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
%% calculate the boarder of whisking period based on RS-eff
%%
exp_session_array = {'am', 'pm'};
for mouse_id = 1:6
    for iexp_session = 1:2
        exp_session = exp_session_array{iexp_session};  
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
        % call function
        WA_eff = lz_labeling_WA_according_RS(mouse_id, exp_session);
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
    end
end
%% compute continuous wavelet coherence for each trial
%%
exp_session_array = {'am', 'pm'};
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
%% unfortunately: out of memory!!
lz_combine_RS_ind_wc            

%% based on the computed wavelet coherence and NA/WA borders, pool index for NA/WA scalegrams
%%
exp_session_array = {'am', 'pm'}; 
for mouse_id = 1:6
    for iexp_session = 1:2
        exp_session = exp_session_array{iexp_session};  
        [wc_NW, wc_WA] = lz_separate_WC_NA_WA_ind(mouse_id, exp_session);
    end
end

%% save NA and WA wavelet coherence separately
exp_session_array = {'am', 'pm'}; 
for mouse_id = 1:6
    for iexp_session = 1:2
        exp_session = exp_session_array{iexp_session};  
        lz_save_separate_wc_NW_WA(mouse_id, exp_session)
    end
end

%% based on the pool index for NA/WA scalegrams, conduct permutation test across subjects for NW vs. WA
%% exclude mouse 3 pm !!!!!!! becasue of sparse WA
pVal = 0.05;
[wc_permut_pval, wc_permut_t_orig, wc_permut_crit_t, wc_permut_est_alpha] ...
    = lz_permutation_NW_WA_2(pVal);

%%%% take entire frequency bins, decide significant changes 
nCh = 30; nChPair = nCh*(nCh-1)/2;
% initial
pVal_th = ones(nChPair, 1);
for iChPair = 1:nChPair
    if min(wc_permut_t_orig(:,iChPair)>0.05)
        pVal_th(iChPair) = 0;
    end
end
        
        
load('f');
figure(12); clf; imagesc(1:nChPair, f, wc_permut_t_orig); 
xlabel('Channel Pair'); ylabel('Frequency (Hz)'); title('Permutation Test Result (t-val)')
set(gca, 'fontsize', 20, 'linew',2);
colorbar

%% 


%% based on the pool index for NA/WA scalegrams, compute mean, std, and t-values (using t-test and permutation test) across NW and WA
exp_session_array = {'am', 'pm'}; 
for mouse_id = 1:6
    for iexp_session = 1:2
        exp_session = exp_session_array{iexp_session};  
        lz_mean_std_t_NW_WA_2(mouse_id, exp_session);
    end
end


        



%%=========================================================================
%% Check point for wavelet coherence --------------------------------------
%% plot two raw calcium signals and their wavelet coherence
mouse_id = 3; exp_session = 'pm'; iTr = 5; iCh = 1; jCh = iCh+5;
lz_checkpoint_raw_calcium_pair_wc(mouse_id, exp_session, iTr, iCh, jCh);

% Check point for non-whisking labeling ----------------------------------
lz_checkpoint_nw_labeling(mouse_id, exp_session, iCh, iTr);

%% check point for identifying boarders of whisking active periods
mouse_id = 6; exp_session = 'am'; iTr = 4;
lz_checkpoint_WA_labeling(mouse_id, exp_session, iTr);

%% Check point for mean, std, and t-values for wavelet coherence across NW/WA
mouse_id = 2; exp_session = 'am';
lz_checkpoint_tval_wc_NW_WA(mouse_id, exp_session);




