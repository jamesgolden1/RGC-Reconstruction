% t_reconTestCV

%%


clear
 
% 
folderNameTrain = 'aug122';
folderNameTest = 'aug122test';
% folderNameTrain = 'aug122prima9';
% folderNameTest = 'aug122prima9test';

pixelWidth = 70/8;
currentDecay = 2;
% folderNameTrain = 'july25prima18';
% % folderNameTest  = 'aug1prima18test';
% folderNameTest = 'july30prima18test';

mosaicFile = 'mosaic0';
windowSize = 1;
percentSV = .5;%.25;%.12;
% shifttime = 2;
shifttime = 17;
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

reconHealthy = recon(pRecon);

mse = reconHealthy.testImagenet(pRecon);
