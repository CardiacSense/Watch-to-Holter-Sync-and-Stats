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

p = corrcoef(TP,truePulse);
p = p(1,2);
[fitresult, gof] = createFit(TP, truePulse);
slope = fitresult.p1;
r2 = gof.rsquare;
%% False Positive:
% ------------------
% All of the test points which are not paired to a Positive - are False
% Positive.
FP=  setdiff(testPnt,TP);
FP=unique(FP);
%%
%gen table
Time = sort([Positive,FP]);
testTime = Time;
testTime(ismember(Time,truePulse)) = TP;
testFlag = -1*ones(size(Time));
testFlag(ismember(testTime,TP)) = 1;
testFlag(ismember(testTime,FP)) = 1;
trueFlag = -1*ones(size(Time));
trueFlag(ismember(Time,truePulse)) = 1;
trueFlag(ismember(Time,FN)) = 1;
trueRR = -1*ones(size(Time));
trueRR(ismember(Time,sort([truePulse FN]))) = refRR;
testRR = -1*ones(size(Time));
testRR(ismember(testTime,sort([TP FP]))) = RR;
testNoise = zeros(size(Time));
testNoise(ismembertol(testTime,Flags,margin,'DataScale',1)) = 1; %how much tol should i use?
%true noise????

%% Sensitivity and FDR
RR = RR(ismember(TP,testPnt));
FN=FN(~ismembertol(FN,Flags,margin,'DataScale',1));
FP = FP(~ismember(FP,Flags));
Sensitivity = (length(TP)/(length(TP)+length(FN)))*100;
FDR         = (length(FP)/(length(FP)+length(TP)))*100;
ppv       = (length(TP)/(length(FP)+length(TP)))*100;
HR = 60000./RR;
trueHR = 60000./refRR(ismember(truePulse,Positive));
avgDist = mean(abs((HR(:) - trueHR(:))),'omitnan');
