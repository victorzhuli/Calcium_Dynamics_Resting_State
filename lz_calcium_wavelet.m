%%%% 06/21/2017
%%%% developed from lz_waveletPower.m
%%%% save time-frequency power patterns for each trial.

function [dbconverted, freqIdx] = lz_calcium_wavelet(ts, srate, freq2plot, chan2plot, tr2plot, timeIdx, baselineidx)

%%%% condition: '_stim_L', '_stim_R', '_stim_LN_corrct', '_stim_RN_correct'
%%%% channel: e.g. 'Fp1'
%%%% freq2plot: e.g. [2 50]
%%%% baseline: e.g. [-300 -50] (unit: ms)

%%%% load data
nPt = size(ts,2);
nTr = size(ts,3);

% define baseline period
% baselinetime = baseline; % in ms

% wavelet parameters
min_freq = freq2plot(1);
max_freq = freq2plot(2);
num_frex = 50;

% other wavelet parameters
freqIdx = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet     = length(time);
n_data        = nPt*nTr;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4; 

% FFT of data (note: this doesn't change on frequency iteration)
fft_data = fft(reshape(ts(chan2plot,:,:),1,[]),n_conv_pow2);

% initialize output time-frequency data
tf_data = zeros(length(freqIdx),nPt,nTr);

for fi=1:length(freqIdx)
    
    % create wavelet and get its FFT
    wavelet = (pi*freqIdx(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*freqIdx(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*freqIdx(fi)))^2))/freqIdx(fi);
    fft_wavelet = fft(wavelet,n_conv_pow2);
    
    % run convolution
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
    convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result_fft = reshape(convolution_result_fft,nPt,nTr);
    
    % put power data into time-frequency matrix
    tf_data(fi,:,:) = abs(convolution_result_fft).^2;
end
% 
% % convert baseline window time to indices
% [~,baselineidx(1)]=min(abs(EEG.times-baselinetime(1)));
% [~,baselineidx(2)]=min(abs(EEG.times-baselinetime(2)));
% 
% initial output
dbconverted = NaN(num_frex,nPt,nTr);
% 
% dB-correct
for iTr = 1: nTr
    baseline_power = mean(tf_data(:,baselineidx(1):baselineidx(2),iTr),2);
    dbconverted(:,:,iTr) = 10*log10( bsxfun(@rdivide,tf_data(:,:,iTr),baseline_power));
end
% the following lines of code are equivalent to the previous line:
% dbconverted = 10*( bsxfun(@minus,log10(tf_data),log10(baseline_power)));
% dbconverted = 10*log10( tf_data ./ repmat(baseline_power,1,EEG.pnts) );
% dbconverted = 10*( log10(tf_data) - log10(repmat(baseline_power,1,EEG.pnts)) );

% %%%% alternatively, use percentage convert
% dbconverted = 100 * (tf_data-repmat(baseline_power,1,EEG.pnts))./ repmat(baseline_power,1,EEG.pnts);

% %%%% alternatively, use baseline division
% dbconverted = tf_data ./ repmat(baseline_power,1,EEG.pnts);


% figure(41);clf;
% contourf(timeIdx,freqIdx,dbconverted(:,:,tr2plot),40,'linecolor','none')
% set(gca,'ytick',round(logspace(log10(freqIdx(1)),log10(freqIdx(end)),10)*100)/100,'yscale','log')
% % title([subject, ' ', condition])
% colorbar