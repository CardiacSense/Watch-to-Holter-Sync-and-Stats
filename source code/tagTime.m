function [ppgTags,ecgTags] = tagTime(SamplingTime,fsCsEcg)
ecgTags = [];
timeTagsCS = datetime(SamplingTime,'InputFormat','dd/MM/yyyy HH:mm:ss:SSS','Format','HH:mm:ss:SSS');
[h,m,s] = hms(timeTagsCS);
d = [0;diff(day(timeTagsCS))];
ms = 1000*rem(s,1);
s = s - rem(s,1);
msFix = round((ms/4)*(1000/64));
for n=2:length(d)
    d(n)=d(n)+d(n-1);
end

ppgTags = msFix+s*1000+m*60*1000+h*60*60*1000+d*24*60*60*1000;

switch fsCsEcg
    case 256
        ecgCsDiff = diff(ppgTags)/4;
        ecgCsTimeIdxTemp1 = ppgTags +[ecgCsDiff;ecgCsDiff(end)];
        ecgCsTimeIdxTemp2 = ecgCsTimeIdxTemp1 +[ecgCsDiff;ecgCsDiff(end)];
        ecgCsTimeIdxTemp3 = ecgCsTimeIdxTemp2 +[ecgCsDiff;ecgCsDiff(end)];
        
        ecgTags = zeros(4*length(ppgTags),1);
        ecgTags(1:4:end) = ppgTags;
        ecgTags(2:4:end) = ecgCsTimeIdxTemp1;
        ecgTags(3:4:end) = ecgCsTimeIdxTemp2;
        ecgTags(4:4:end) = ecgCsTimeIdxTemp3;
    case 128
        ecgCsDiff = diff(ppgTags)/2;
        ecgCsTimeIdxTemp1 = ppgTags +[ecgCsDiff;ecgCsDiff(end)];
        ecgTags = zeros(2*length(ppgTags),1);
        ecgTags(1:2:end) = ppgTags;
        ecgTags(2:2:end) = ecgCsTimeIdxTemp1;
end
end