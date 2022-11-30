%% Stats from csv
clear all;
close all;
clc;

Sensitivity=[nan;nan];
FDR=[nan;nan];
avgDist=[nan;nan];

rrSensitivity=[nan;nan];
rrFDR=[nan;nan];
%%
Tag={'PPG','ECG'};
filesPath='C:\Users\adi\OneDrive - CardiacSense\Desktop\General\Recordings\UOM\';
PP='';
for iSig=[1,2]
    [infoF,infoP]=uigetfile([filesPath '\*.csv'],['Choose ' Tag{iSig} ' result file']);
    filesPath=infoP;
    resT=readtable([infoP,infoF]);

    testRR=table2array(resT(:,4));
    testFlag=table2array(resT(:,5));
    testNoise=table2array(resT(:,10));

    trueRR=table2array(resT(:,8));
    trueFlag=table2array(resT(:,9));

    % %% Peaks
    TP=find(testFlag==1 & trueFlag==1);
    FP=find(testFlag==1 & trueFlag==-1 & testNoise==0);
    FN=find(testFlag==-1 & trueFlag==1 & testNoise==0);

    Sensitivity(iSig) = (length(TP)/(length(TP)+length(FN)))*100;
    FDR(iSig)         = (length(FP)/(length(FP)+length(TP)))*100;

    rrTP = 0;
    rrFN = 0;
    rrFP = 0;
    trueRR_inx=[];
    FN_k=[];
    for k = 2:length(testRR)
        if ~testNoise(k) && testRR(k)~=0 % valid RR
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
                    FN_k=[FN_k k];
                end
                if ( (trueFlag(k)==-1 || trueFlag(k-1)==-1) && (testFlag(k)==1 && testFlag(k-1)==1)) ...
                        || ...
                        ( (trueFlag(k)==1 && trueFlag(k-1)==1) && (testFlag(k)==1 && testFlag(k-1)==-1) )
                    % FP
                    % wrong RR - both end for watch but not on holter
                    %            OR
                    %            both end for hlter but left end for watch
                    rrFP = rrFP + 1;
                end
            end
        end
    end
    rrSensitivity(iSig) = rrTP/(rrTP+rrFN)*100;
    rrFDR(iSig)         = rrFP/(rrFP+rrTP)*100;

    trueRR4calc_ref=trueRR(trueRR_inx);
    trueRR4calc_test=testRR(trueRR_inx);

    HR = 60000./trueRR4calc_test;
    trueHR = 60000./trueRR4calc_ref;
    trueHR(trueHR<0) = nan;
    avgDist(iSig) = sqrt((sum((HR(:) - trueHR(:)).^2,'omitnan'))/length(HR)); % RMS
end
%%
TotalResults=table(Sensitivity,FDR,rrSensitivity,rrFDR,avgDist,'RowNames',{'PPG','ECG'});


%% compare results file
[resF,resP]=uigetfile([infoP '\*.csv'],'Choose results summary file');
totResT=readtable([resP,resF]);

PPG_Sensitivity = totResT.PPGSensitivity;
PPG_FDR         = totResT.PPGFDR;
ECG_Sensitivity = totResT.ECGSensitivity;
ECG_FDR         = totResT.ECGFDR;
PPG_RR_Sensitivity = totResT.PPGRRSensitivity;
PPG_RR_FDR         = totResT.PPGRRFDR;
ECG_RR_Sensitivity = totResT.ECGRRSensitivity;
ECG_RR_FDR         = totResT.ECGRRFDR;
PPG_RMS = totResT.PPGRMSBPM;
ECG_RMS = totResT.ECGRMSBPM;

%%
fprintf('\n')
fprintf('PPG Sensitivity distance= %0.4f\n',round(Sensitivity(1),4)-round(PPG_Sensitivity,4));
fprintf('PPG FDR distance= %0.4f\n',round(FDR(1),4)-round(PPG_FDR,4));
fprintf('PPG avgDist distance= %0.4f\n',round(avgDist(1),4)-round(PPG_RMS,4));
fprintf('PPG RR Sensitivity distance= %0.4f\n',round(rrSensitivity(1),4)-round(PPG_RR_Sensitivity,4));
fprintf('PPG RR FDR distance= %0.4f\n',round(rrFDR(1),4)-round(PPG_RR_FDR,4));
fprintf('\n')
fprintf('ECG Sensitivity distance= %0.4f\n',round(Sensitivity(2),4)-round(ECG_Sensitivity,4));
fprintf('ECG FDR distance= %0.4f\n',round(FDR(2),4)-round(ECG_FDR,4));
fprintf('ECG avgDist distance= %0.4f\n',round(avgDist(2),4)-round(ECG_RMS,4));
fprintf('ECG RR Sensitivity distance= %0.4f\n',round(rrSensitivity(2),4)-round(ECG_RR_Sensitivity,4));
fprintf('ECG RR FDR distance= %0.4f\n',round(rrFDR(2),4)-round(ECG_RR_FDR,4));

%% One file for all recordings

recDir=uigetdir();
recFlist=dir(recDir);
recFlist(~[recFlist.isdir],:)=[];
recFlist([1,2],:)=[];

recNo={recFlist.name};
recNo=sort(cellfun(@str2num,recNo,'UniformOutput',true));

splitFecg=@(STR) strsplit(STR,'reportECG.csv');
splitFppg=@(STR) strsplit(STR,'reportPPG.csv');

first=1;
for iRec=recNo
    currDir=dir([recDir '\' num2str(iRec)]);
    currDir(~[currDir.isdir],:)=[];

    if iRec==5 
        currFilePath=[[recDir '\' num2str(iRec)] '\' currDir(end-1).name];
    else
        currFilePath=[[recDir '\' num2str(iRec)] '\' currDir(end).name];
    end
    currDirIn=dir(currFilePath);

    ppgName=currDirIn((cellfun(@length,cellfun(splitFppg,{currDirIn.name},'UniformOutput',false),'UniformOutput',true)==2)).name;
    ecgName=currDirIn((cellfun(@length,cellfun(splitFecg,{currDirIn.name},'UniformOutput',false),'UniformOutput',true)==2)).name;
    
    if first==1
        first=0;
        currPPG_T=readtable(fullfile(currFilePath,ppgName),'ReadRowNames',false);
        currECG_T=readtable(fullfile(currFilePath,ecgName),'ReadRowNames',false);

        writetable(currPPG_T,[recDir '\AllRecReportPPG.csv'],'WriteRowNames',true);
        writetable(currECG_T,[recDir '\AllRecReportECG.csv'],'WriteRowNames',true);
    else
        currPPG_T=readtable(fullfile(currFilePath,ppgName),'ReadRowNames',false);
        currECG_T=readtable(fullfile(currFilePath,ecgName),'ReadRowNames',false);

        writetable(currPPG_T,[recDir '\AllRecReportPPG.csv'],'WriteMode','append','WriteRowNames',true);
        writetable(currECG_T,[recDir '\AllRecReportECG.csv'],'WriteMode','append','WriteRowNames',true);
    end
end

%%
temp3=ismember(temp2(2:end),temp2(1:end-1));

temp3=temp2(2:end)==temp2(1:end-1);

