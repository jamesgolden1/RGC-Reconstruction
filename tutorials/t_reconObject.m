% 
% 
% 
% pTrain.mosaicFile = 'mosaic_nsPhys';
% pTrain.saveFile   = 'nsPhys';
% % [mosaicFile, saveFile] = trainNaturalScenesPhys(pTrain);
% 
% % trainHealthy = recon.build(pTrain);
% % recon.mosaicFile
% % recon.saveFile / recon.buildFile
% 
% mosaicFile = '_mosaicAll_413415'; 
% movieFile = 'nsPhys/nsPhys_mov';
% spikesFile = 'nsPhys/nsPhys_sp';
% filterFile = ['nsPhys/filt1_sv99_' mosaicFile];
% 
% pLoad.loadFile = ['nsPhys/' pTrain.saveFile];
% pLoad.movieFile = movieFile;
% pLoad.spikesFile = spikesFile;
% pLoad.mosaicFile = mosaicFile;
% % 
% [movieFile, spikesFile] = loadSpikesAll(pLoad);
% % % % % 
% pRecon.movieFile = movieFile;
% pRecon.spikesFile = spikesFile;
% pRecon.filterFile = filterFile;
% % % % % 
% pRecon.mosaicFile = mosaicFile;
% pRecon.windowSize = 1;
% pRecon.percentSV = 0.99;
% pRecon.shiftTime = 4;
% [filterFile] = runReconstructSVD_fast_all(pRecon);
% 
% % trainHealthy = recon.train(pTrain);
% % recon.buildFile
% % recon.movieFile /recon.stimFile
% % recon.spikesFile /recon.respFile
% % recon.mosaicFile
% % recon.filterFile
% % recon.windowSize
% % recon.percentSV
% % recon.shiftTime
% 
% pTest.mosaicFile = mosaicFile;
% pTest.filterFile = filterFile;s
% testReconNS(pTest);
% 
% % testHealthy = recon.test(pTest)


%%



mosaicFile = '_mosaicAll_1246640' ;%['_' d1(dind).name(11:end-4)];% 
movieFile = 'ns100_r2_10/ns100_jan1_mov3'; 
spikesFile = 'ns100_r2_10/ns100_jan1_sp3';

evArr = [.2 .1 .3 .4];
trainFraction = [.16 .5 .66 .83 1];
evInd = 4; trainFractionInd = 5;
filterFile = ['ns100_r2_10/filters4_ns100_feb6_sh9_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];

clear pRecon
pRecon.buildFile = 'test0';
pRecon.stimFile = movieFile;
pRecon.respFile = spikesFile;
pRecon.filterFile = filterFile;

pRecon.mosaicFile = mosaicFile;
pRecon.windowSize = 1;
pRecon.percentSV = 0.99;

reconHealthy = recon(pRecon);

reconHealthy.build(pRecon);

% reconHealthy.train(pRecon);
% 
% reconHealthy.plot('filters');
% 
% reconHealthy.test(pRecon);