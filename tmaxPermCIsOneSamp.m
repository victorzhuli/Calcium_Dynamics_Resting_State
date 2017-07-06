function [CIs, critTmax]=tmaxPermCIsOneSamp(data,nPerm,span,reports)
% function [CIs, critTmax]=tmaxPermCIsOneSamp(data,nPerm,span,reports)
%  Computes symmetric multiple comparison adjusted confidence intervals (CIs)  
%  of the means of a multivariate set of observations by permuting the residuals 
%  of the data (Groppe, in press). The frequentist probability of one or more CI
%  missing the true mean value is 1-span (if your data are approximately
%  Gaussian).
%
% Usage:
%  >> [CIs, critTmax]=tmaxPermCIsOneSamp(data,nPerm,span,reports);
%
% Required Input:
%  data   - 2D matrix of data (Observation x Variable)
%
% Optional Inputs:
%  nPerm       - Number of permutations used to estimate the complete set
%                of permutations. If the number of observations is less
%                than or equal to 12, all possible permutations are used and
%                this optional input has no effect. If the number of 
%                observations is greater than 12, nPerm specifies the 
%                number of random permutations computed. Manly (1997) 
%                suggests using at least 1000 permutations for an alpha level 
%                of 0.05 and at least 5000 permutations for an alpha level of
%                0.01. {default=5000}
%  span        - Desired CI coverage. Note, because of the finite
%                number of possible permutations, the exact desired coverage
%                may not be possible for small samples. {default=.95}
%  report      - [1 or 0] If non-zero, text updates of the function's progress
%                are output to the MATLAB command line {default: 1}
%
%
% Outputs:
%   CIs        - Matrix of confidence interval bounds. First row is lower
%                bound. Second row is upper bound.
%   critTmax   - The t-score used to derive the confidence interval bounds.
%
%
% References:
%   Groppe, D.M. (in press) Combating the scientific decline effect with 
%   confidence (intervals). Psychophysiology.
%   http://biorxiv.org/content/biorxiv/early/2015/12/10/034074.full.pdf
%
%   Manly, B.F.J. (1997) Randomization, Bootstrap, and Monte Carlo Methods in
%   Biology. 2nd ed. Chapman and Hall, London.
%
%
% Example:
%  data=randn(12,5);
%  CIs=tmaxPermCIsOneSamp(data);
%  disp('95 percent CI lower bound:');
%  disp(CIs(1,:));
%  disp('95 percent CI upper bound:');
%  disp(CIs(2,:));
%
%
% See also mult_comp_perm_t1.m:
%   http://www.mathworks.com/matlabcentral/fileexchange/29782-mult-comp-perm-t1
% For generating p-values analogous to the CIs computed here.
%
%
% Author:
% David M. Groppe
% Feinstein Institute for Medical Research
% Manhasset, NY
% Dec 19, 2015


if nargin<2,
    nPerm=5000;
end

if nargin<3,
    span=.95;
end

if nargin<4,
    reports=1;
end

[nObs, nVar]=size(data);
if reports,
    fprintf('tmaxPermCIsOneSamp: Number of variables: %d\n',nVar);
    fprintf('tmaxPermCIsOneSamp: Number of observations: %d\n',nObs);
    fprintf('t-score degrees of freedom: %d\n',nObs-1);
end

% Remove the mean of the data
mnObs=mean(data);
seObs=std(data)/sqrt(nObs);

rsdls=data-repmat(mnObs,nObs,1); % Residuals

%% Set up permutations
if nObs<=12,
    nPerm=2^nObs; %total number of possible permutations
    exact=1;
    seed_state='exact';
    if reports,
        fprintf('Due to the limited number of observations, all possible permutations of the data will be computed instead of random permutations.\n');
    end
else
    exact=0;
end

if reports,
    fprintf('Executing permutation test with %d permutations...\n',nPerm);
    fprintf('Permutations completed: ');
end


%% Compute permutations

%Constant factor for computing t, speeds up computing t to precalculate
%now
sqrt_nXnM1=sqrt(nObs*(nObs-1));
if exact,
    %Use all possible permutations
    mxt=zeros(1,nPerm);
    for perm=1:nPerm
        if reports,
            if ~rem(perm,100)
                if ~rem(perm-100,1000)
                    fprintf('%d',perm);
                else
                    fprintf(', %d',perm);
                end
                if ~rem(perm,1000)
                    fprintf('\n');
                end
            end
        end
        %set sign of each participant's data
        if perm==1
            bvec=zeros(ceil(log2(nPerm)),1);
        else
            bvec=incBinaryVec(bvec);
        end
        sn=2*bvec-1;
        sn_mtrx=repmat(sn,1,nVar); 
        d_perm=rsdls.*sn_mtrx;
        
        %computes t-score of permuted data across all channels and time points
        sm=sum(d_perm,1);
        mn=sm/nObs;
        sm_sqrs=sum(d_perm.^2,1)-(sm.^2)/nObs;
        stder=sqrt(sm_sqrs)/sqrt_nXnM1;
        t=mn./stder;
        
        %get most extreme t-score
        [dummy mxt_id]=max(abs(t));
        mxt(perm)=t(mxt_id); %get the most extreme t-value with its sign (+ or -)
    end
else
    %Use random permutations
    mxt=zeros(1,nPerm*2);
    for perm=1:nPerm
        if reports,
            if ~rem(perm,100)
                if ~rem(perm-100,1000)
                    fprintf('%d',perm);
                else
                    fprintf(', %d',perm);
                end
                if ~rem(perm,1000)
                    fprintf('\n');
                end
            end
        end
        %randomly set sign of each participant's data
        sn=(rand(nObs,1)>.5)*2-1; 
        sn_mtrx=repmat(sn,1,nVar);
        
        d_perm=rsdls.*sn_mtrx;
        
        %computes t-score of permuted data across all channels and time points
        sm=sum(d_perm,1);
        mn=sm/nObs;
        sm_sqrs=sum(d_perm.^2,1)-(sm.^2)/nObs;
        stder=sqrt(sm_sqrs)/sqrt_nXnM1;
        t=mn./stder;
        
        %get most extreme t-score (sign isn't immportant since we asumme
        %symmetric distribution of null hypothesis for one sample test)
        mxt(perm)=max(abs(t));
    end
    mxt(nPerm+1:2*nPerm)=-mxt(1:nPerm); %add the negative of all values since we assumme
    %null hypothesis distribution is symmetric
end

%End permutations completed line
if reports && rem(perm,1000)
    fprintf('\n');
end


%% Get critical percentile of tmax distribution
alph=1-span;
%critTmax=prctile(mxt,(1-alph/2)*100); % Original line uses in between vals
critTmax=prctile(mxt,(1-(1/nPerm+alph)/2)*100); % can get exact vals

% Compute t-distribution based CIs
CIs=zeros(2,nVar);
CIs(1,:)=mnObs-critTmax*seObs; % Lower bound
CIs(2,:)=mnObs+critTmax*seObs; % Upper bound


%%%% Helper Function %%%%
function bvec=incBinaryVec(bvec)
% Increments a binary vector (representing a base 2 number) by 1 

nEl=length(bvec);
ct=1;
bvec(ct)=bvec(ct)+1;
while bvec(ct)>1
    bvec(ct)=0;
    ct=ct+1;
    if ct>nEl,
        error('Cannot increment binary vector. Max value/length reached.');
    end
    bvec(ct)=bvec(ct)+1;
end
% fprintf('%d: ',a);
% disp(bvec);
