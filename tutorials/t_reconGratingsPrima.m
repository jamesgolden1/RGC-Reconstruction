% t_reconGratingsPrima
% 
% Generate training data for full-field gratings responses.

clear

%  contrastArr =  .5*[.00002 .0002 .0008 .001 .0016 002 .0035 .005 .008 .02 .04 .06 .08];

% % % % % contrastArr = .001% 5*[0 .00002 .0002 .001 .002 .0025 .003 .004 .005 .008 .02 .04 .06 .08];
%  contrastArr = .5*[.00002 .0002 .0008 .001 .0015 002 .0035 .005 .007 .0085 .009 .0095 .008 .02 .04 .06 .08];
% pool = parpool(length(contrastArr));

freqArr = [.05 .1 .2 .5 1 2 4 5 8 10 16];
% pool = parpool(length(freqArr));

% parfor contrastInd = 1:length(contrastArr)%.5*[0 .02 .04 .06 .08 .11 .15]
    
for freqInd = 1:length(freqArr)
% for contrast =.5* [0  .04  .08 .15]
% %     for contrast =.5* [.02 .06 .11 ]
        reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);
     reconPrimaLandolt.buildPrimaGratings('nTrials',100,'contrast',.01,'gratingSpFreq',freqArr(freqInd),'horizontalFlag',true);

end
delete(gcp);