function [avgDist,Sensitivity, rrSensitivity,FDR, rrFDR,p,slope,r2,pkData] = peaksBasedPerformanceAnalysis (Positive,testPnt,RR,RRflag,margin,Flags,lag,Title)
%% Analysis of peaks.
% Benchmark - True Peaks.
% 

% Positive  - Benchmark peaks are regerded as positives.
% testPnt   - Output of algorithm (watch)
% RR        - RR output from the watch
% RRflag    - RR flag- 1 for valid RR and 0 otherwise
% margin    - margin for peaks match
% Flags     - noise flag - watch output
% lag       - Lag between watch & ref
% Title     - Title for plots

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
noisePeak=testPnt(RR==(1000*1/128));
testPnt(RR==(1000*1/128))=[];              % peaks for coupling
RRflag(RR==(1000*1/128))=[];
RR(RR==(1000*1/128))=[];

refRR = [-1 diff(Positive)];
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

[TP, IA] = unique(TP);
%keep unique TPs and assign others as FN...
FN = sort([FN, setdiff(truePulse,truePulse(IA))]);
truePulse = truePulse(IA);

p = corrcoef(TP,truePulse);
p = p(1,2);
try
    [fitresult, gof] = createFit(TP, truePulse);
    slope = fitresult.p1;
    r2 = gof.rsquare;
catch
    slope = nan; r2 = nan;
    warning('Unsuccessful fitting');
end
%% False Positive:
% ------------------
% All of the test points which are not paired to a Positive - are False
% Positive.
[FP,FP_inx] =  setdiff(testPnt,TP);
FP = unique(FP);

%% Gen table
whiteFN=FN(ismembertol(FN,noisePeak,margin,'DataScale',1));

Time = unique(sort([Positive(:);FP(:)]));
testTime = Time; 
testTime(ismember(Time,truePulse)) = TP;
testTimeOrig = testTime - (lag);
testFlag = -1*ones(size(Time));                                  % watch peaks
testFlag(ismember(testTime,TP)) = 1;
testFlag(ismember(testTime,FP)) = 1;
trueFlag = -1*ones(size(Time));                                  % reference peaks
trueFlag(ismember(Time,truePulse)) = 1;
trueFlag(ismember(Time,FN)) = 1;
trueRR = -1*ones(size(Time));                                    % reference RR 
trueRR(ismember(Time,sort([truePulse(:); FN(:)]))) = refRR;
testRR = -1*ones(size(Time));                                    % watch RR
testRR(ismember(testTime,sort([TP(:); FP(:)]))) = RR;
testNoise = zeros(size(Time));
testNoise(ismembertol(testTime,Flags,margin,'DataScale',1)) = 1; %how much tol should i use?
testNoise(ismember(testTime,whiteFN)) = 1;                       % invalid peaks printed for debug use
refNoise = -1*ones(size(Time));
RRnoise = ones(size(Time));                                     % unvalid RR flag (watch), 1-unvalid, -1 - false, 0-valid
RRnoise(ismember(testRR,RR(find(~RRflag))))=0;
RRnoise(ismember(testTime,FN(:)))=0;

%true noise????
%% RR stats TP/FP/FN
% visual explanation in pptx in folder
rrTP = 0;
rrFN = 0;
rrFP = 0;
trueRR_inx=[]; 
for k = 2:length(testRR)
    if ~testNoise(k) &&  RRnoise(k)~=1 % valid RR
        if  trueFlag(k)==1 && trueFlag(k-1)==1  && testFlag(k)==1 && testFlag(k-1)==1
            % TP
            % RR within margin and peaks on both ends were correctly detected
            rrTP = rrTP + 1;
            trueRR_inx=[trueRR_inx k];
        else
            if ( (trueFlag(k)==1 && trueFlag(k-1)==1) && (testFlag(k)==-1 || testFlag(k-1)==-1) ) ...
                    || ...
                    ( (trueFlag(k)==1 && trueFlag(k-1)==-1) && (testFlag(k)==1 && testFlag(k-1)==1))
                % FN
                % missing RR- peaks on both ends for holter but not on both ends for watch
                %             OR
                %             both end for watch but not for holter
                rrFN = rrFN + 1;
                RRnoise(k)=-1; % remove wrong RR
            end
            if ( (trueFlag(k)==-1 || trueFlag(k-1)==-1) && (testFlag(k)==1 && testFlag(k-1)==1)) ...
                    || ...
                    ( (trueFlag(k)==1 && trueFlag(k-1)==1) && (testFlag(k)==1 && testFlag(k-1)==-1) )
                % FP
                % wrong RR - both end for watch but not on holter
                %            OR
                %            both end for hlter but left end for watch
                rrFP = rrFP + 1;
                RRnoise(k)=-1; % remove wrong RR
            end
        end
    end
end
rrSensitivity = rrTP/(rrTP+rrFN)*100;
rrFDR         = rrFP/(rrFP+rrTP)*100;
rrppv         = rrTP/(rrFP+rrTP)*100;

trueRR4calc_ref=trueRR(trueRR_inx);
trueRR4calc_test=testRR(trueRR_inx);

% 
pkData = table(testTime(:),testTimeOrig(:),testRR(:),testFlag(:),...
                refNoise(:),Time(:),trueRR(:),trueFlag(:),testNoise(:));%,RRnoise(:));
%pkData=renamevars(pkData,1:10,["testTime","testTimeOrig","testRR","testFlag","refNoise","Time","trueRR","trueFlag","testNoise","RRnoise"]);
%% Sensitivity and FDR
refFlagPnts = ismembertol(FN,[Flags(:);noisePeak(:)],margin,'DataScale',1);
noisePnts = FN(refFlagPnts);
FN=FN(~refFlagPnts);
FP = FP(~ismember(FP,Flags));
Sensitivity = (length(TP)/(length(TP)+length(FN)))*100;
FDR         = (length(FP)/(length(FP)+length(TP)))*100;
ppv         = (length(TP)/(length(FP)+length(TP)))*100;
HR = 60000./trueRR4calc_test;
trueHR = 60000./trueRR4calc_ref;
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