function [tstat, p, cv] = npperm(A, B, iter, correction)

% Written by Tylor J Harlow 8/6/2022.  Analytic assistance from Matthew B Jane
% and credit to Scott C Bressler for benhoch function.
% This function performs non-parametric permutation testing of
% time series data. The function is designed to operate at two
% potential levels: individual or group.
% 
% Inputs:
%     A & B : Cell arrays, of length N subjects, each of which is a 
%     numeric array organized as [trials x time].
% 
%     Iter: number of permutation iterations desired. Performed equally at
%     the level of the individual anf the group.
%     
%     Correction: Scalar isequal to 0 or 1, indiciating the output of the 
%     Benjamini-Hochberg Procedure to address false discovery rate (FDR)



% Setup
if nargin < 3
    iter = 100;
end
if nargin < 4
    correction = 0;
end
N = length(A);

% Shuffle Trial Labels at the individual subject level
for j = 1:iter
    for i = 1:N
        Atrials = size(A{i},1);
        Btrials = size(B{i},1);
        AllTrials = (1:(Atrials + Btrials))';
        subjectAll = cat(1,A{i},B{i});
        [~,idx] = datasample(AllTrials,Atrials,1,"Replace",true);
        [~,idx2] = datasample(AllTrials,Btrials,1,"Replace",true);
        subjectcondA(i,:) = squeeze(mean(subjectAll(idx,:)));
        subjectcondB(i,:) = squeeze(mean(subjectAll(idx2,:)));
        if j ==1
            truA(i,:) = nanmean(A{i},1);
            truB(i,:) = nanmean(B{i},1);
        else
            break
        end
    end

    % Shuffle Trial Labels at the group level
    All = cat(1,subjectcondA,subjectcondB);
    groupcount = (1:size(All,1))';
    [~,idx] = datasample(groupcount,N,1,"Replace",true);
    [~,idx2] = datasample(groupcount,N,1,"Replace",true);
    groupA = All(idx,:);
    groupB = All(idx2,:);

    
    % Establish null t-distribution
    diff = groupA - groupB;
    tstat(j,:) = mean(diff)./(std(diff)/sqrt(N));  
end

% Significance Test
truDiff = truA - truB;
truT = mean(truDiff)./(std(truDiff)/sqrt(N));
z = (truT - mean(tstat))./std(tstat);  p = (1-normcdf(abs(z)))*2;

% Correction
if correction == 1
    [~, cv] = benhoch(p);
    factor = 0.05/cv;
    p = p*factor;
end



%=================== Benhoch ========================================%
function [h,CV] = benhoch(p,FDR)
% [h,CV] = benhoch(p,FDR)
%   Benjamini-Hochberg Procedure to address false discovery rate (FDR)
%
% INPUT VARIABLES
%     p : vector of individual p-values
%   FDR : false discovery rate (default = 0.05)
%
% OUTPUT VARIABLE
%   h : Benjamini-Hochberg adjusted hypothesis test
%       [0=failure to reject null hypothesis, 1=reject null hypothesis]
%  CV : critical (p) value at which to reject null hypothesis
%
% Created: 2018-Sep-07 SCB
% Revised: 2019-May-17 SCB addressed bug in hypothesis assignment

if nargin<2
    FDR = 0.05; % default FDR value
end

[pI,I] = sort(p(:)); % sort p-values
m = numel(pI); % total number of tests
r = (1:m)'; % inidividual p-value ranks
CVs = (r/m)*FDR; % Benjamini-Hochberg critical value
h0 = pI<CVs; % find highest p-value < Benjamini-Hochberg critical value
h0(1:find(h0==1,1,'last')) = 1; % reset all higher ranks to h = 1
CV = CVs(find(h0==1,1,'last'));
h(I) = h0; % re-assign hypothesis test results
h = reshape(h,size(p)); % reshape vector to match original input 'p'





