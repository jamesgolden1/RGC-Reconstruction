clear

%%

folderNameTrain = 'aug27';% 'aug30';
% folderNameTest = 'sep13gratings4_add2';

folderNameTest = 'healthy_gratings_aug16_cont';
hFlag = 1;

% folderNameTest = 'healthy_gratings_aug16_nocont';
% hFlag = 0;

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

currentDecay = 2;

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
% freqArr =[2 4 8.0000 12 16.0000   18.1595   20.6104   23.3921   26.5493   30.1326   34.1995 38.8153  42 44.0541  47 50.0000];
freqArr = [.05 .5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 8 12 25];
% pool = parpool(length(freqArr));
    
% contrastArr = 1*[0 .00002 .0002 .0003:.00004:.0004 .001 .002 .005 .008 .02 .04 .06 .08];
% contrastArr = [.009:.001:.02];
% contrastArr = [.01125 .0115 .01175 1.0];%.012
% contrastArr = 1;%*[0 .00001 .00002 .00004 .00008 .0001 .0002 .0004 .001 .002:.001:.004 .006 .02 .05 .1 ];
% 


% contrastArr = 1*[0 .00001 .0001 .0002 .0004 .008 .001 .002 .003 .004 .006 .02 .05 .1 .15 .2];
% contrastArr = [.3 .5 .66 .75 .88 .95 1];%
% contrastArr = [1e-7 2e-6 4e-6 8e-6 9e-6 1e-5 1.25e-5 1.5e-5 1.75e-5 2e-5 4e-5 8e-5 1e-4 4e-4 1e-3 .01];
% contrastArr = [2e-4 3e-4 5e-4 6e-4 7e-4 8e-4 9e-4 2e-3 3e-3 4e-4 6e-4 8e-3 .02 .03];
% contrastArr = [.02 .04 .08 .1 .125 .15 .2 .25 .33 .4 .5 .75];
pRecon.rngInput = 510573;
parpool(length(freqArr));

pRecon.horizontalFlag = hFlag;
  
pRecon.nTrials = 200;
pRecon.contrast = 1;
parfor freqInd = 1:length(freqArr)

%     for contrastInd = length(contrastArr)
    reconGratings = recon();
%     reconGratings.buildGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);

     reconGratings.buildGratings('gratingSpFreq',freqArr(freqInd), pRecon);
end
delete(gcp);
% contrastArr =  [0 .1]
