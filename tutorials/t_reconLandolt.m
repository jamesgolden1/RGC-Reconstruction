clear

%%

folderNameTrain = 'aug30';
% folderNameTest = 'sep13gratings4_add2';

folderNameTest = 'sep17landoltFreq';

pixelWidth = 70/8;
currentDecay = 2;
pulseFreq = 
pulseDutyCycle
irradianceFraction
currentDecay
% 
% addParameter(p,  'pulseFreq',0,@isscalar);
% addParameter(p,  'pulseDutyCycle',1,@isscalar);
% addParameter(p,  'irradianceFraction',1,@isscalar);
% addParameter(p,  'currentDecay',2,@isscalar);
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

% freqArr = [.05 .1 .2 .5 1 2 4 5 8 10 16];
% pool = parpool(length(freqArr));
    
% contrastArr = 1*[0 .00002 .0002 .0003:.00004:.0004 .001 .002 .005 .008 .02 .04 .06 .08];
% contrastArr = [.009:.001:.02];
% contrastArr = [.01125 .0115 .01175];%.012
% contrastArr = 1*[0 .00002 .0002 .0004 .001 .002:.001:.006 .008 .02:.01:.06];
% 
% 
contrastArr = .1 ;%1*[0 .0002 .0004 .001 .002:.002:.01 .02:.02:.1 .15 .2];

resizeArr = sqrt([.0025 .01 .02 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .6 .8 1]);
% 

% 

% gap = 20;
gap = -2;

parpool(length(resizeArr));

%  hFlag = 1%[0 1]
parfor resizeInd = 1:length(resizeArr)

for contrastInd = length(contrastArr)   
%     for contrastInd = length(contrastArr)
    reconLandolt = recon();
%     reconGratings.buildGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);

%      reconGratings.buildGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',freqArr(freqInd),'horizontalFlag',hFlag, pRecon);
       reconLandolt.buildLandolt('nTrials',100,'contrast',contrastArr(contrastInd),'gap',gap,'resizeFraction',resizeArr(resizeInd),pRecon);   

end
    end
delete(gcp);
% contrastArr =  [0 .1]
