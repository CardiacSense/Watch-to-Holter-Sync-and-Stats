function lagFix  = snap2data(Positive,testPnt)
testDiff = [0 ;diff(testPnt(:))];
posDiff = [0;diff(Positive(:))];
%%
[nIdx,nDist] = knnsearch([testPnt(:),testDiff(:)],[Positive(:),posDiff(:)],"Distance","euclidean");
fSigPeaks = testPnt(nIdx);
LM = fitlm(fSigPeaks,Positive,'linear','RobustOpts','on');
if LM.Coefficients.pValue(1)<=0.01
    lagFix = LM.Coefficients.Estimate(1);
else
    lagFix = 0;
end