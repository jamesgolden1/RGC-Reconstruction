
clear

% t_reconGratingsPrima
% 
% Generate training data for full-field gratings responses.

clear

%  contrastArr =  .5*[.00002 .0002 .0008 .001 .0016 002 .0035 .005 .008 .02 .04 .06 .08];

% % % % % contrastArr = .001% 5*[0 .00002 .0002 .001 .002 .0025 .003 .004 .005 .008 .02 .04 .06 .08];
%  contrastArr = .5*[.00002 .0002 .0008 .001 .0015 002 .0035 .005 .007 .0085 .009 .0095 .008 .02 .04 .06 .08];
% pool = parpool(length(contrastArr));

freqArr = 4;%[.05 .1 .2 .5 1 2 4 5 8 10 16];
% pool = parpool(length(freqArr));

% parfor contrastInd = 1:length(contrastArr)%.5*[0 .02 .04 .06 .08 .11 .15]
  
folderNameTrain = 'aug29prima70';
% folderNameTest = 'sep12prima70gratings_add';
folderNameTest = 'sep20prima70landoltContAdd';
% folderNameTest = 'sep14landolt';

pixelWidth = 70/1;
currentDecay = 2;
% folderNameTrain = 'july25prima18';
% % folderNameTest  = 'aug1prima18test';
% folderNameTest = 'july30prima18test';

mosaicFile = 'mosaic0';
windowSize = 1;
percentSV = .05;%.5;%.5;%.25;%.12;
% shifttime = 2;
shifttime = 4;%15;%3;%15;
dropout = 0;

filterFile  = fullfile(folderNameTrain,...    
    ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);

% pRecon.pixelWidth = pixelWidth;
pRecon.currentDecay = currentDecay;
pRecon.mosaicFile = mosaicFile;
pRecon.filterFile = filterFile;
pRecon.stimFile = folderNameTest;
pRecon.folderNameTest = folderNameTest;
pRecon.folderNameTrain = folderNameTrain;
pRecon.windowSize = windowSize;
pRecon.percentSV = percentSV;
pRecon.dropout = dropout; 

% resizeArr = [.5 1];
% contrastArr = [0 1e-5 .0001 .0003 .001 .003 .01 .015 .02 .025 .03 .035 .04 .05 .06 .1 ];
contrastArr = [.07 .08 .09 .095 .105 .11 .115 .12 .125 .13 .14 .2];
contrastArr = [.0004 .0008 .002 .004 .006 .008 .009 .0095 .00975 .01125 .0125 .01375 .01625 .0175 .01875 .0225];
% contrastArr = [.06 .07 .09:.01:.2];

% resizeArr = sqrt([.0025 .01 .02 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .6 .8 1]);
resizeArr = 1;%sqrt([ .01 .1 .2 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .875 1]);

parpool(length(contrastArr));

hFlag = 0;

% gap = 20-2;
gap = -2;
pRecon.rngInput = 1533420;
for resizeInd = length(resizeArr)
    
parfor contrastInd = 1:length(contrastArr)   
%     parfor hFlag = [0 1]
% for contrast =.5* [0  .04  .08 .15]
% %     for contrast =.5* [.02 .06 .11 ]
        reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);
%      reconPrimaLandolt.buildPrimaGratings('nTrials',400,'contrast',contrastArr(contrastInd),'gratingSpFreq',freqArr(freqInd),'horizontalFlag',hFlag,pRecon);
        reconPrimaLandolt.buildPrimaLandolt('nTrials',200,'contrast',contrastArr(contrastInd),'gap',gap,'resizeFraction',resizeArr(resizeInd),pRecon);   
end
end
delete(gcp);

%%

% contrastArr =  [.002 .008 .02 .04 .06 .08];
% pool = parpool(length(contrastArr));
% for gap = [-2]
%     parfor contrastInd = 1:length(contrastArr)
% 
%         reconLandolt = recon();
%         reconLandolt.buildLandolt('nTrials',200,'contrast',contrastArr(contrastInd),'gap',gap);
%     end
% end

% for contrast = [.02 .04 .06 .08 .11 .15]
%     for gap = [0 ]
%         reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaLandolt('nTrials',200,'contrast',contrast,'gap',gap);
%     end
% end

% for contrast = [.02 .04 .06 .08 .11 .15]
%     for gap = [2]
%         reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaLandolt('nTrials',200,'contrast',contrast,'gap',gap);
%     end
% end
%         