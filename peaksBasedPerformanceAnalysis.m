function [avgDist,Sensitivity,FDR,p,slope,r2] = peaksBasedPerformanceAnalysis (Positive,testPnt,RR,margin,Flags)
%% Analysis of peaks.
% Benchmark - True Peaks.
% testPnt   - Output of algorithm.

% Positive  - Benchmark peaks are regerded as positives.
% Negative  - 

% Disclaimers:
% ----------------
% (1)   If there are two peak paired to a positive (for example 4 ,7 paired
%       to 6 with margin 3 [3,9]) the one with the closest arithmetic
%       distatnce is chosen (7 will be chosen). This way only one TP can be
%       paired to a Positive.
% (2)   This version is immune to doubled value like [3,3,3, 11,...]. It
%       treats it as one. 
% (3)   
%% Unique Tests
% Take off any duplicate values.
Positive = unique (Positive);
testPnt = unique (testPnt);
refRR = [diff(Positive) nan];
%% True Positive & False Nega
% -----------------
TP  =[];                            % True  Positve
FN  =[];                            % False Negative
truePulse = [];
for k=1:length(Positive)
    %% Positive range:
    TPbound = [Positive(k)-margin,Positive(k)+margin];
    
    % For every Positive range [Positive - margin, Positive + margin] we
    % are searching for a point in its range.
    Pair=testPnt(testPnt>=TPbound(1) & testPnt<=TPbound(2));
    
    % Only one peak can be paried to a Positive, this closest one is
    % chosen.
    if(length(Pair)>1)
        [~,ind]=min(abs(Positive(k)- Pair));
        Pair = Pair(ind);
    end
    
    if(Pair)
        TP= [TP Pair];              % True Positive
        truePulse = [truePulse Positive(k)];
    else
        FN = [FN  Positive(k)];     % False Negative
    end
end

TP = unique(TP);
truePulse = unique(truePulse);
RR = RR(ismember(TP,testPnt));
p = corrcoef(TP,truePulse);
p = p(1,2);
[fitresult, gof] = createFit(TP, truePulse);
slope = fitresult.p1;
r2 = gof.rsquare;
FN=FN(~ismember(FN,Flags));
%% False Positive:
% ------------------
% All of the test points which are not paired to a Positive - are False
% Positive.
FP=  setdiff(testPnt,TP);
FP= FP(~ismember(FP,Flags));
FP=unique(FP);

%% Sensitivity and FDR
Sensitivity = (length(TP)/(length(TP)+length(FN)))*100;
FDR         = (length(FP)/(length(FP)+length(TP)))*100;
HR = 60000./RR;
trueHR = 60000./refRR(ismember(truePulse,Positive));
avgDist = mean(abs((HR(:) - trueHR(:))),'omitnan');
