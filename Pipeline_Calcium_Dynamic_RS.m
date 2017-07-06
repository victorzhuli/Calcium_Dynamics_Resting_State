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
pVal = 0.01;
[wc_permut_pval, wc_permut_t_orig, wc_permut_crit_t, wc_permut_est_alpha] ...
    = lz_permutation_NW_WA_2(pVal);
% plot resulted p values for all freq bins and channel pairs
load('f'); nCh = 30; nChPair = nCh*(nCh-1)/2;
figure(12); clf; %imagesc(1:nChPair, f, wc_permut_t_orig); 

surf(1:nChPair, f, zeros(size(wc_permut_t_orig)),'CData',wc_permut_t_orig, 'Linestyle','none')
view(0,90)
set(gca,'yscale','log')
xlim([1 nChPair]); ylim([0 f(1)])
xlabel('Channel Pair'); ylabel('Frequency (Hz)'); title('Permutation Test Result (t-val)')
set(gca, 'fontsize', 20, 'linew',2);
colorbar

%%%% take entire frequency bins, decide significant changes 
% initial
t_orig_pos = zeros(nChPair, 1);
t_orig_neg = zeros(nChPair, 1);
t_orig_pos_freq = cell(nChPair, 1); % where does the peak t value appear
t_orig_neg_freq = cell(nChPair, 1);
for iChPair = 1:nChPair
    if max(wc_permut_t_orig(:,iChPair)) > wc_permut_crit_t(2,iChPair)
        t_orig_pos(iChPair) = 1;
        [~,ind] = max(wc_permut_t_orig(:,iChPair));
        t_orig_pos_freq(iChPair) = f(ind);
    end
    if min(wc_permut_t_orig(:,iChPair)) < wc_permut_crit_t(1,iChPair)
        t_orig_neg(iChPair) = 1;
        [~,ind] = min(wc_permut_t_orig(:,iChPair));
        t_orig_neg_freq(iChPair) = f(ind);
    end
end
t_orig_pos_freq = find(wc_permut_t_orig > 
% plot channel connection matrix
diff_mat_pos = nan(nCh, nCh);
diff_mat_neg = nan(nCh, nCh);
diff_mat_pos_freq  = nan(nCh, nCh);
diff_mat_neg_freq  = nan(nCh, nCh);
for iCh = 1: nCh
    diff_mat_pos(iCh, iCh) = 0;
    diff_mat_neg(iCh, iCh) = 0;
    diff_mat_pos_freq(iCh, iCh) = 0;
    diff_mat_neg_freq(iCh, iCh) = 0;
    for jCh = iCh + 1 : nCh
        ind_wcoh = lz_ind_loc_wcoh(iCh, jCh, nCh);      
        diff_mat_pos(iCh, jCh) = t_orig_pos(ind_wcoh);
        diff_mat_pos(jCh, iCh) = diff_mat_pos(iCh, jCh);
        diff_mat_pos_freq(iCh, jCh) = t_orig_pos_freq(ind_wcoh);
        diff_mat_pos_freq(jCh, iCh) = diff_mat_pos_freq(iCh, jCh);
        diff_mat_neg(iCh, jCh) = t_orig_neg(ind_wcoh);
        diff_mat_neg(jCh, iCh) = diff_mat_neg(iCh, jCh);
        diff_mat_neg_freq(iCh, jCh) = t_orig_neg_freq(ind_wcoh);
        diff_mat_neg_freq(jCh, iCh) = diff_mat_neg_freq(iCh, jCh);
    end
end
figure(14); clf; 
subplot(121);imagesc(diff_mat_pos); title('FC pattern (NW > WA)')
colormap('gray');axis square; colorbar;
subplot(122);imagesc(diff_mat_neg); title('FC pattern (NW < WA)')
colormap('gray');axis square; colorbar;

figure(15); clf;
subplot(121); imagesc(diff_mat_pos_freq); title('The associated freq bin on which peak tVal obtained for NW>WA')
colormap('jet');axis square; c = colorbar; ylabel(c,'Hz');
subplot(122); imagesc(diff_mat_neg_freq); title('The associated freq bin on which peak tVal obtained for NW<WA')
colormap('jet');axis square; c = colorbar; ylabel(c,'Hz');
        
  
        


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




