% t_reconTestCV

%%


clear
 
% 
% folderNameTrain = 'aug27';
% folderNameTest = 'aug27test';

folderNameTrain = 'aug27prima70';
folderNameTest = 'aug27prima70test';

% folderNameTrain = 'aug23sp';
% folderNameTest = 'aug23sptest';
% folderNameTrain = 'august/aug122prima70';
% folderNameTest = 'august/aug122prima70test';

pixelWidth = 70/2;
currentDecay = 2;
% folderNameTrain = 'july25prima18';
% % folderNameTest  = 'aug1prima18test';
% folderNameTest = 'july30prima18test';

mosaicFile = 'mosaic0';
windowSize = 1;
percentSV = .05;%.25;%.12;
% shifttime = 2;
shifttime = 4;%15;%4;%17;
dropout = 0;

filterFile  = fullfile(folderNameTrain,...    
    ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout) '_aug27']);

pRecon.pixelWidth = pixelWidth;
pRecon.currentDecay = currentDecay;
pRecon.mosaicFile = mosaicFile;
pRecon.filterFile = filterFile;
pRecon.stimFile = folderNameTest;
pRecon.windowSize = windowSize;
pRecon.percentSV = percentSV;
pRecon.dropout = dropout; 

reconHealthy = recon(pRecon);

% mseArr = reconHealthy.testImagenet(pRecon);
