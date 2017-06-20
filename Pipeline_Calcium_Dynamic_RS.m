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
%% label whisker and non-whisker segments
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
        Ind_RSstart{iTr} = [1, Ind_RSstart{iTr}];
    end
    if whisker_conv_win(iTr,end) == 0
        Ind_RSend{iTr} = [Ind_RSend{iTr}, size(whisker_conv_win,2)];
    end
end
    
% set a threshold (th) where only keep non-whisker period that is longer than the th
th_second = 5;
th = th_second * sampRate_whisk/w_st;
% initiate the start point of RS to be kept
Ind_RSstart_eff = cell(iTr,1);
for iTr = 1:nTr
    for iPeriod = 1: length(Ind_RSstart{iTr})
        period = Ind_RSend{iTr}(iPeriod) - Ind_RSstart{iTr}(iPeriod);
        if period >= th
            Ind_RSstart_eff{iTr} = [Ind_RSstart_eff{iTr}, Ind_RSstart{iTr}(iPeriod)];
        end
    end
end
            


%--------------------------------------------------------------------------
% Check point
figure(1);clf;
subplot(411);plot(Cal(8,:,1),'b');
subplot(412);plot(anglekeeper(2,:)); xlim([0 length(anglekeeper(2,:))])
subplot(413);plot(whisker_conv_win(2,:),'r')
subplot(414);plot(diff_whisker_flag(2,:),'r')
%--------------------------------------------------------------------------





