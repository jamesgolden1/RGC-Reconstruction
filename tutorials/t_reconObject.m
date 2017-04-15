


pTrain.mosaicFile = 'mosaic_nsPhys';
pTrain.saveFile   = 'nsPhys';
[mosaicFile, saveFile] = trainNaturalScenesPhys(pTrain);

% trainHealthy = recon.build(pTrain);
% recon.mosaicFile
% recon.saveFile / recon.buildFile

mosaicFile = '_mosaicAll_413415'; 
movieFile = 'nsPhys/nsPhys_mov';
spikesFile = 'nsPhys/nsPhys_sp';
filterFile = ['nsPhys/filt1_sv99_' mosaicFile];

pLoad.loadFile = ['nsPhys/' pTrain.saveFile];
pLoad.movieFile = movieFile;
pLoad.spikesFile = spikesFile;
pLoad.mosaicFile = mosaicFile;
% 
[movieFile, spikesFile] = loadSpikesAll(pLoad);
% % % % 
pRecon.movieFile = movieFile;
pRecon.spikesFile = spikesFile;
pRecon.filterFile = filterFile;
% % % % 
pRecon.mosaicFile = mosaicFile;
pRecon.windowSize = 1;
pRecon.percentSV = 0.99;
pRecon.shiftTime = 4;
[filterFile] = runReconstructSVD_fast_all(pRecon);

% trainHealthy = recon.train(pTrain);
% recon.buildFile
% recon.movieFile /recon.stimFile
% recon.spikesFile /recon.respFile
% recon.mosaicFile
% recon.filterFile
% recon.windowSize
% recon.percentSV
% recon.shiftTime

pTest.mosaicFile = mosaicFile;
pTest.filterFile = filterFile;s
testReconNS(pTest);

% testHealthy = recon.test(pTest)