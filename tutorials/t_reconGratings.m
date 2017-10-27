clear

%%

folderNameTrain = 'aug30';
folderNameTest = 'sep16gratingsFreq';

pixelWidth = 70/8;
currentDecay = 2;
% folderNameTrain = 'july25prima18';
% % folderNameTest  = 'aug1prima18test';
% folderNameTest = 'july30prima18test';

mosaicFile = 'mosaic0';
windowSize = 1;
percentSV = .5;%.5;%.25;%.12;
% shifttime = 2;
shifttime = 15;%15;
dropout = 0;

filterFile  = fullfile(folderNameTrain,...    
    ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);

% pRecon.pixelWidth = pixelWidth;
pRecon.currentDecay = currentDecay;
pRecon.mosaicFile = mosaicFile;
pRecon.filterFile = filterFile;
pRecon.stimFile = folderNameTest;
pRecon.windowSize = windowSize;
pRecon.percentSV = percentSV;
pRecon.dropout = dropout; 
pRecon.folderNameTrain = folderNameTrain;
pRecon.folderNameTest = folderNameTest;
%%

% contrastArr =  .05*.5*[0 .00002 .0002 .0003:.00002:.0004 .001 .002 .005 .008 .02 .04 .06 .08];
% 
% % contrastArr =  .05*.5*[.0003:.00002:.0004];
% % 
% % contrastArr = .5*[ .0025 .003 .004 ];
% % contrastArr =  .5*[0  .008 .04 .08];
% % contrastArr = .5* [.002  .02 .06 .1];
% pool = parpool(length(contrastArr));

% freqArr = 2;%2%[.05 .1 .2 .5 1 2 4 5 8 10 16];
freqArr =[2 4 8.0000 12 16.0000   18.1595   20.6104   23.3921   26.5493   30.1326   34.1995 38.8153  42 44.0541  47 50.0000];

% pool = parpool(length(freqArr));
    
% contrastArr = 1*[0 .00002 .0002 .0003:.00004:.0004 .001 .002 .005 .008 .02 .04 .06 .08];
% contrastArr = [.009:.001:.02];
% contrastArr = [.01125 .0115 .01175 1.0];%.012
contrastArr = .1;%*[0 .00001 .00002 .00004 .00008 .0001 .0002 .0004 .001 .002:.001:.004 .006 .02 .05 .1 ];
% 

% parpool(length(contrastArr));

% parpool(length(freqArr));

 hFlag =0%[0 1]
% par
for freqInd = 1:length(freqArr)

% par
for contrastInd = length(contrastArr)   
%     for contrastInd = length(contrastArr)
    reconGratings = recon();
%     reconGratings.buildGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);

     reconGratings.buildGratings('nTrials',400,'contrast',contrastArr(contrastInd),'gratingSpFreq',freqArr(freqInd),'horizontalFlag',hFlag, pRecon);
end
    end
% delete(gcp);
% contrastArr =  [0 .1]
