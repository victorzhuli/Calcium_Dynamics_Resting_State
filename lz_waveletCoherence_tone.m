function [fv_lk_freq1, fv_lk_freq2, fv_nl_freq1, fv_nl_freq2] = lz_waveletCoherence_tone(caldata, lklabel, SRATE)
% wtc analysis for one subject, 03/29/2018
% @Parameters:
% caldata: calcium data, matrix shape (time x channel x trial)
% lklabel: licking/non-licking labels, cell shape ([l.start,l.end,nl.start,nl.end] x trial)
% @Returns:
% fv_lk_freq1: feature vectors for licking in freq band [0.5,1]
% fv_lk_freq2: feature vectors for licking in freq band (1,2]
% fv_nl_freq1: feature vectors for non-licking in freq band [0.5,1]
% fv_nl_freq2: feature vectors for non-licking in freq band (1,2]

NUM_CH = 30; NUM_CH_PAIR = NUM_CH * (NUM_CH - 1) / 2;
NUM_TR = size(caldata, 3);

% inital feature vectors
fv_lk_freq1 = nan(NUM_TR, NUM_CH_PAIR); % 0.5 - 1 Hz
fv_lk_freq2 = nan(NUM_TR, NUM_CH_PAIR); % 1 - 2 Hz
fv_nl_freq1 = nan(NUM_TR, NUM_CH_PAIR);
fv_nl_freq2 = nan(NUM_TR, NUM_CH_PAIR);

for itr = 1 : NUM_TR
    
    % extract licking/non-licking inital/end labels
    lk.ini = lklabel{1, itr};
    lk.end = lklabel{2, itr};
    nl.ini = lklabel{3, itr};
    nl.end = lklabel{4, itr};
    
    num_lk = length(lk.ini); num_nl = length(nl.ini);
    % compute wc
    iChPair: % index on the feature vector
    for ich = 1 : NUM_CH
        for jch = ich + 1: NUM_CH
            [wcoh, ~, f] = wcoherence(caldata(:, ich, itr), caldata(:, jch, itr), SRATE);
            freq1_idx = find(f >= 0.5 & f <= 1);
            freq2_idx = find(f > 1 & f <= 2);
            
            % compute the feature vector based on wc
            % licking conditions:
            lk_freq1_tmp = nan(num_lk, 1); lk_freq2_tmp = nan(num_lk, 1);
            for ilk = 1 : num_lk
                % freq1:
                lk_freq1_tmp(ilk) = mean(mean(wcoh(freq1_idx, lk.ini(ilk):lk.end(ilk))));
                fv_lk_freq1(itr, iChPair) = mean(lk_freq1_tmp);
                % freq2:
                lk_freq2_tmp(ilk) = mean(mean(wcoh(freq2_idx, lk.ini(ilk):lk.end(ilk))));
                fv_lk_freq2(itr, iChPair) = mean(lk_freq2_tmp);
            end
            
            % non-licking conditions:
            nl_freq1_tmp = nan(num_nl, 1); nl_freq2_tmp = (num_nl, 1);
            for inl = 1 : num_nl
                % freq1:
                nl_freq1_tmp(inl) = mean(mean(wcoh(freq1_idx, nl.ini(inl):nl.end(inl))));
                fv_nl_freq1(itr, iChPair) = mean(nl_freq1_tmp);
                % freq2:
                nl_freq2_tmp(inl) = mean(mean(wcoh(freq2_idx, nl.ini(inl):nl.end(inl))));
                fv_nl_freq2(itr, iChPair) = mean(nl_freq2_tmp);
            end
            
            iChPair = iChPair + 1;
        end
    end
end
           
    
            
            
            
            
            
            