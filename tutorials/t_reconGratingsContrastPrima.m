% t_reconGratingsPrima
% 
% Generate training data for full-field gratings responses.

clear

%  contrastArr =  .5*[.00002 .0002 .0008 .001 .0016 002 .0035 .005 .008 .02 .04 .06 .08];

% % % % % contrastArr = .001% 5*[0 .00002 .0002 .001 .002 .0025 .003 .004 .005 .008 .02 .04 .06 .08];
%  contrastArr = .5*[.00002 .0002 .0008 .001 .0015 002 .0035 .005 .007 .0085 .009 .0095 .008 .02 .04 .06 .08];
% pool = parpool(length(contrastArr));

freqArr = 2;%[.05 .1 .2 .5 1 2 4 5 8 10 16];
% pool = parpool(length(freqArr));

% parfor contrastInd = 1:length(contrastArr)%.5*[0 .02 .04 .06 .08 .11 .15]
  
folderNameTrain = 'aug29prima70';
folderNameTest = 'sep19prima70gratingsContAdd';

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

% contrastArr = 10*[0 .00002 .0002 .0003:.00004:.0004 .001 .002 .005 .008 .02 .04 .06 .08];
% 
% contrastArr = [.06 .07 .09:.01:.2 .2];

% contrastArr = 1*[0 .0002 .0004 .001 .002:.002:.01 .02:.02:.1 .15 .2];
% contrastArr = .05;%.1;
% freqArr =[2 4 8.0000 12 16.0000   18.1595   20.6104   23.3921   26.5493   30.1326   34.1995 38.8153  42 44.0541  47 50.0000];
% freqArr = [.05 .5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 8 12 25];


% contrastArr = 1*[0 .0001 .001 .008 .01 .0125 .015 .0175 .02 .0225 .025 .0275 .03 .05 .1 .15];

contrastArr = .1;%1*[0 .00004 .0001 .0002 .0003 .0004 .0006 .0008 .001 .002 .004 .006 .008 .01 .02 .05 ];
contrastArr = [ .0008 .0009 .003 .005 .007 .009  ];
% contrastArr = .001;%[.002 .003 .004 .005 .006 .007];
% parpool(18);
parpool(length(contrastArr));

hFlag =0;
for freqInd = length(freqArr)
    
parfor contrastInd = 1:length(contrastArr)   
%     parfor hFlag = [0 1]
% for contrast =.5* [0  .04  .08 .15]
% %     for contrast =.5* [.02 .06 .11 ]
        reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);
     reconPrimaLandolt.buildPrimaGratings('nTrials',200,'contrast',contrastArr(contrastInd),'gratingSpFreq',freqArr(freqInd),'horizontalFlag',hFlag,pRecon);
    end
end
delete(gcp);