function [avgDist,Sensitivity, rrSensitivity,FDR, rrFDR,p,slope,r2,pkData] = peaksBasedPerformanceAnalysis (Positive,testPnt,RR,margin,Flags,lag,Title)
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
Positive = unique(Positive);
testPnt = unique(testPnt);
refRR = [diff(Positive) -1];
%% True Positive & False Nega
% -----------------
TP  =[];                            % True  Positve
FN  =[];                            % False Negative
truePulse = [];
falseRR_FN=[];
trueRR_ref_inx=[];

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
        trueRR_ref_inx=[trueRR_ref_inx,k];                           % Reference index for TP

    else
        FN = [FN  Positive(k)];     % False Negative
        falseRR_FN=[falseRR_FN find(testPnt>Positive(k),1,'first')]; % First peak after FN - false RR
    end
end

[TP, IA] = unique(TP);
%keep unique TPs and assign others as FN...
FN = sort([FN, setdiff(truePulse,truePulse(IA))]);
truePulse = truePulse(IA);
trueRR_ref_inx = trueRR_ref_inx(IA);

p = corrcoef(TP,truePulse);
p = p(1,2);
[fitresult, gof] = createFit(TP, truePulse);
slope = fitresult.p1;
r2 = gof.rsquare;
%p=[];slope=[];r2=[];

% slope = nan;
% r2 = nan;
%% False Positive:
% ------------------
% All of the test points which are not paired to a Positive - are False
% Positive.
[FP,FP_inx] =  setdiff(testPnt,TP);
FP = unique(FP);
falseRR_FP=FP_inx+1;

%%
%gen table
Time = sort([Positive(:);FP(:)]);
testTime = Time; 
testTime(ismember(Time,truePulse)) = TP;
testTimeOrig = testTime - (lag);
testFlag = -1*ones(size(Time));
testFlag(ismember(testTime,TP)) = 1;
testFlag(ismember(testTime,FP)) = 1;
trueFlag = -1*ones(size(Time));
trueFlag(ismember(Time,truePulse)) = 1;
trueFlag(ismember(Time,FN)) = 1;
trueRR = -1*ones(size(Time));
trueRR(ismember(Time,sort([truePulse(:); FN(:)]))) = refRR;
testRR = -1*ones(size(Time));
testRR(ismember(testTime,sort([TP(:); FP(:)]))) = RR;
testNoise = zeros(size(Time));
testNoise(ismembertol(testTime,Flags,margin,'DataScale',1)) = 1; %how much tol should i use?
refNoise = -1*ones(size(Time));
pkData = table(testTime(:),testTimeOrig(:),testRR(:),testFlag(:),...
                refNoise(:),Time(:),trueRR(:),trueFlag(:),testNoise(:));
%true noise????
%% RR TP/FP/FN
rrTP = 0;
rrFN = 0;
rrFP = 0;
for k = 2:length(testRR)
    if trueFlag(k) && trueFlag(k-1)  && testFlag(k) && testFlag(k-1)
        % RR within margin and peaks on both ends were correctly detected
        rrTP = rrTP + 1;
    elseif trueFlag(k) && trueFlag(k-1) && (~testFlag(k) || ~testFlag(k-1))
        % peaks on both ends for holter but not on both ends for watch ->
        % FN
        rrFN = rrFN + 1;
    else
        %either RR not within margin or peaks detected on watch but not
        %holter (either end)
        rrFP = rrFP + 1;
    end
end
rrSensitivity = rrTP/(rrTP+rrFN)*100;
rrFDR         = rrFP/(rrFP+rrTP)*100;
rrppv         = rrTP/(rrFP+rrTP)*100;
%% Sensitivity and FDR
TrueRRinx=find(ismember(testPnt,TP));      % True RR (relate to TP)
TrueRRinx=setdiff(TrueRRinx,falseRR_FP);   % remove false RR's relate to FP
TrueRRinx=setdiff(TrueRRinx,falseRR_FN);   % remove false RR's relate to FN
RR = RR(TrueRRinx);                        % Valid RR

trueRR_ref_inx(falseRR_FP)=[];             % remove false RR's relate to FP
trueRR_ref_inx(falseRR_FN)=[];             % remove false RR's relate to FN

refFlagPnts = ismembertol(FN,Flags,margin,'DataScale',1);
noisePnts = FN(refFlagPnts);
FN=FN(~refFlagPnts);
FP = FP(~ismember(FP,Flags));
Sensitivity = (length(TP)/(length(TP)+length(FN)))*100;
FDR         = (length(FP)/(length(FP)+length(TP)))*100;
ppv         = (length(TP)/(length(FP)+length(TP)))*100;
HR = 60000./RR;
trueHR = 60000./refRR(trueRR_ref_inx);
trueHR(trueHR<0) = nan;
% avgDist = mean(abs((HR(:) - trueHR(:))),'omitnan'); % MAE
avgDist = sqrt((sum((HR(:) - trueHR(:)).^2,'omitnan'))/length(HR)); % RMS
%% Gen new figure 
figure(); hold on;
plot(Positive,ones(size(Positive)),'Color','m','Marker','hexagram','DisplayName','Holter Peaks','MarkerSize',12,'LineStyle','none','LineWidth',1);
plot(TP,ones(size(TP)),'xr','DisplayName',['True Positive ',num2str(numel(TP))],'MarkerSize', 15,'LineWidth',1);
plot(FN,ones(size(FN)),'sk','DisplayName',['False Neg ', num2str(numel(FN))],'MarkerSize',15,'LineWidth',1);
plot(FP,ones(size(FP)),'og','DisplayName',['False Positive ', num2str(numel(FP))],'MarkerSize',15,'LineWidth',1);
plot(noisePnts,ones(size(noisePnts)),'b*','DisplayName','Noise Time','MarkerSize',15,'LineWidth',1);
legend show;
title(Title)


%%
%{
figure; plot(Positive,refRR,'*-','DisplayName','ref'); hold on; 
plot(testPnt,RR,'*-','DisplayName','CS'); 
plot(TP,RR((ismember(testPnt,TP))),'og','DisplayName','TP');
plot(FN,refRR((ismember(Positive,FN))),'om','DisplayName','FN');
plot(FP,RR((ismember(testPnt,FP))),'oy','DisplayName','FP');
WW=setdiff(1:length(RR),TrueRRinx);
plot(testPnt(WW),RR(WW),'sk','DisplayName','false RR');
legend show
%}