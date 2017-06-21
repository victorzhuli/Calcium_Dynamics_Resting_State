%% label whisker and non-whisker segments
function RS_eff = lz_labeling_nonWhisking(mouse_id, exp_session)

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
        folder_name,'/',mouse_name,'/anglekeeper.mat']); % whisker signals
elseif isunix
    cd('/home/lz206/Dropbox/projects/calcium/Calcium_Dynamics_Resting_State')
    load(['/home/lz206/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/anglekeeper.mat']); % whisker signals
end

threshold = 25;

if (mouse_id == 2 && strcmp(exp_session, 'am')) || (mouse_id == 3 && strcmp(exp_session, 'pm')) 
    nTr = 15;
else
    nTr = 16;
end

sampRate_whisk = 500;
w_ln = 150 ; % whisker window length
w_st = 25 ;  % whisker window step
nLn  = 10000; % length of whisker time-series
nWn  = (nLn - w_ln) / w_st+1; % # total window

% Initiation. 0 will be non-whiskering, 1 will be whisking
whisker_conv_win = nan(nTr, nWn);
for iTr = 1: nTr
    for iWn = 1:nWn
        widx = [(iWn-1)*w_st+1, (iWn-1)*w_st+w_ln];
        data2compute = anglekeeper(iTr,widx(1):widx(end));
        data2compute = inpaint_nans(data2compute, 5); 
            % interpolate NaN points using 8 neighbor average.
        whisker_conv_win(iTr, iWn) = ...
            std(data2compute - mean(data2compute)); % variance
    end
end
whisker_conv_win(whisker_conv_win<threshold) = 0;
whisker_conv_win(whisker_conv_win>threshold) = 1;
diff_whisker_flag = diff(whisker_conv_win');
diff_whisker_flag = diff_whisker_flag';

% find the labels where non-whiskering is start and end
Ind_RSstart = cell(nTr,1);  Ind_RSend = cell(nTr,1);  
for iTr = 1:nTr
    Ind_RSstart{iTr} = find(diff_whisker_flag(iTr,:) == -1);
    Ind_RSend{iTr}   = find(diff_whisker_flag(iTr,:) == 1);
    % detect whisking state at block beginning
    if whisker_conv_win(iTr,1) == 0
        Ind_RSstart{iTr} = [0, Ind_RSstart{iTr}];
    end
    if whisker_conv_win(iTr,end) == 0
        Ind_RSend{iTr} = [Ind_RSend{iTr}, size(whisker_conv_win,2)];
    end
end
    
% set a threshold (th) where only keep non-whisker period that is longer
% than "th"
th_second = 4;
th = th_second * sampRate_whisk/w_st;
% initiate the start point of RS to be kept
RS_eff.ind_start = cell(iTr,1);
RS_eff.ind_end   = cell(iTr,1);
for iTr = 1:nTr
    for iPeriod = 1: length(Ind_RSstart{iTr})
        period = Ind_RSend{iTr}(iPeriod) - Ind_RSstart{iTr}(iPeriod);
        if period >= th
            RS_eff.ind_start{iTr} = [RS_eff.ind_start{iTr}, Ind_RSstart{iTr}(iPeriod)];
            RS_eff.ind_end{iTr}   = [RS_eff.ind_end{iTr},   Ind_RSend{iTr}(iPeriod)];
        end
    end
end

% project whisker index to calcium index
% all indices should add one because they were associated to
% diff_whisker_flag
for iTr = 1:nTr
    RS_eff.ind_start{iTr} = RS_eff.ind_start{iTr} + 1;
    RS_eff.ind_end{iTr}   = RS_eff.ind_end{iTr} + 1;
end

for iTr = 1:nTr
    % one point in whisker_conv_win is associated with w_st point in whisker
    % signal and w_st/5 with Calcium signal ... /5 because sampling rate of
    % whisker signal is 5 times larger than that of Calcium signal.
    RS_eff.ind_start_cal{iTr} = (RS_eff.ind_start{iTr} - 1) * w_st/5 + 1;
    RS_eff.ind_end_cal{iTr}   = (RS_eff.ind_end{iTr} - 1) * w_st/5 + 1;
    % move every end index to (w_ln-1)/5 points backward
    RS_eff.ind_end_cal{iTr} = RS_eff.ind_end_cal{iTr} - (w_ln-1)/5;
    if whisker_conv_win(iTr,end) == 0 && Ind_RSstart{iTr}(end) < (10000 - th_second * 500)/w_st % 500 is the sampling rate of whisker signal
        RS_eff.ind_end_cal{iTr}(end) = size(2000,2); % 2000 is the length of Calcium signal  
    end
end