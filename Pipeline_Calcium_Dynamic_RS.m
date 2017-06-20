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
        folder_name,'/',mouse_name,'/Ca.mat']); % Calcium signals
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/', ...
        folder_name,'/',mouse_name,'/anglekeeper.mat']); % whisker signals
elseif isunix
    
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

w_ln = 150 ; % whisker window length
w_st = 25 ;  % whisker window step
nLn  = 10000; % length of whisker time-series
nWn  = (nLn - w_ln) / w_st+1; % # total window

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
%--------------------------------------------------------------------------
% Check point
figure(1);clf;
subplot(411);plot(Cal(8,:,1),'b');
subplot(412);plot(anglekeeper(1,:)); xlim([0 length(anglekeeper(1,:))])
subplot(413);plot(whisker_conv_win(1,:),'r')
subplot(414);plot(diff_whisker_flag(1,:),'r')
%--------------------------------------------------------------------------





