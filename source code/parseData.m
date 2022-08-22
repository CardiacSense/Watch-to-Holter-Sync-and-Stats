function [FsECG,QRSindParsed,ECG,hdr] = parseData(path,file,fileFormat)

if strcmp(fileFormat,'EDF file') %%%%%%%% Parse EDF data %%%%%%%%
    
    % Get EDF data
    [edfData,record]=edfread([path file]);
    hdr = edfData;
    % Extract signals
    
    ecgRow=strfind(edfData.label,'ECG');
    tf= cellfun('isempty',ecgRow); % true for empty cells
    ecgRow(tf)= {0};
    ecgRow=find(cell2mat(ecgRow));
    
    hrvRow=strfind(edfData.label,'HRV');
    tf= cellfun('isempty',hrvRow); % true for empty cells
    hrvRow(tf)= {0};
    hrvRow=find(cell2mat(hrvRow));
    
    FsECG=edfData.frequency(ecgRow);
    FsRR=edfData.frequency(hrvRow);
    
    ECG=record(ecgRow,:);
    RR=find(record(hrvRow,:));
    
    % Find QRS indices
    QRSind=round(cumsum(record(hrvRow,RR(2:end)))/1000* FsECG)+((RR(1)-1)/ FsRR * FsECG);
    
    % Find shift caused by FsRR
    [~,shiftInd]=max(ECG(QRSind(2)-round(0.5*(QRSind(2)-QRSind(1))):QRSind(2)+round(0.5*(QRSind(3)-QRSind(2)))));
    shift=shiftInd-round(0.5*(QRSind(2)-QRSind(1)))+1;
    QRSindParsed=QRSind+shift;
    
    % figure;
    % plot(ECG);hold on;plot(QRSind,ECG(QRSind),'o');
    eventListForParse = [];
    
else %%%%%%%% Parse MAT data %%%%%%%%
    
    % Get MAT data
    load([path file]);
    
    % Extract signals
    
    FsECG = Fs;
    if ~exist('QRSforTool')
        QRSindParsed = sort(QRS);
    else
        QRSindParsed = sort(QRSforTool);
    end
    if exist('eventList')
        eventListForParse = eventList;
    else
        eventListForParse =[];
    end
    
    
end
end