% t_reconTestCV

%%


clear


% folderNameTrain = 'aug30';

% folderNameTrain = 'aug29prima70';
% folderNameTest = 'sep8prima70test';

% folderNameTest = 'sep8prima70test';
% folderNameTest = 'aug29prima70test';

% folderNameTrain = 'aug29prima35';
% folderNameTest = 'sep20prima35test';
% folderNameTest = 'sep20prima35hall';

folderNameTrain = 'aug29prima70';
% % folderNameTest = 'sep20prima70test';
folderNameTest = 'sep20prima70hall';

% folderNameTest = 'sep9prima9hall';
% folderNameTest = 'aug29prima35hall';

% folderNameTrain = 'aug27prima9';
% folderNameTest = 'aug29prima9hallrng';

% folderNameTrain = 'aug23sp';
% folderNameTest = 'aug23sptest';

pixelWidth = 70/8;
currentDecay = 2;
% folderNameTrain = 'july25prima18';
% % folderNameTest  = 'aug1prima18test';
% folderNameTest = 'july30prima18test';

mosaicFile = 'mosaic0';
windowSize = 1;
percentSV = .05;%.05;%.02;%.05;%.5;%.5;%.25;%.12;
% shifttime = 2;
shifttime = 4;%15;%4;%15;%4;%15;%3;%15;
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

mse1 = reconHealthy.testImagenet(pRecon);
