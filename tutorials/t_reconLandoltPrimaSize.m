
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
  
folderNameTrain = 'prosthesis_70_training_aug13';
% folderNameTest = 'sep12prima70gratings_add';

% folderNameTest = 'prosthesis_70_landolt_aug16_gap';
% gap = 20;
folderNameTest = 'prosthesis_70_landolt_aug16_nogap';
gap = -2;

% folderNameTest = 'sep14landolt';

pixelWidth = 70/1;
currentDecay = 2;
% folderNameTrain = 'july25prima18';
% % folderNameTest  = 'aug1prima18test';
% folderNameTest = 'july30prima18test';

mosaicFile = 'mosaic0';
% prosthesis
windowSize = 1;
percentSV = .05;
shifttime = 3;
dropout = .3;

pRecon.pixelWidth = pixelWidth;

filterFile  = fullfile(folderNameTrain,...    
    ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout) '_pitch' sprintf('_%d',pixelWidth) '_decay_2']);


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
contrastVal= .04;%*[0 .0002 .0004 .001 .002:.002:.01 .02:.02:.1 .15 .2];
% 
% contrastArr = [.06 .07 .09:.01:.2];

% resizeArr = sqrt([.0025 .01 .02 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .6 .8 1]);
resizeArr = sqrt([ .01 .1 .2 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .875 1]);
% resizeArr = sqrt([ .01 .1 .2 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .875 1]);

parpool(length(resizeArr));

hFlag = 0;

pRecon.nTrials = 200;
pRecon.contrast = contrastVal;
pRecon.gap = gap;

parfor resizeInd = 1:length(resizeArr)
    
% par
% for contrastInd = length(contrastArr)   
%     parfor hFlag = [0 1]
% for contrast =.5* [0  .04  .08 .15]
% %     for contrast =.5* [.02 .06 .11 ]
        reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);
%      reconPrimaLandolt.buildPrimaGratings('nTrials',400,'contrast',contrastArr(contrastInd),'gratingSpFreq',freqArr(freqInd),'horizontalFlag',hFlag,pRecon);
%         reconPrimaLandolt.buildPrimaLandolt('nTrials',200,'contrast',contrastArr(contrastInd),'gap',gap,'resizeFraction',resizeArr(resizeInd),pRecon);   
        reconPrimaLandolt.buildPrimaLandolt('resizeFraction',resizeArr(resizeInd),pRecon);   

% end
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