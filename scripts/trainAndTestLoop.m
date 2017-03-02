% trainAndTest
% 
% Run the entire pipeline for the RGC linear decoder. 
%   1. trainNaturalScenes - generate training data - show movies to RGC mosaics
%   2. loadSpikesAll      - rearrange spikes
%   3. runReconstruct     - learn decoding filters
%   4. testReconNS        - decode with learned filters and evaluate
%                               accuracy of reconstructions.
% 

% pTrain.mosaicFile = 'mosaic_wn_Dec3_sp_2';
% pTrain.saveFile   = 'ns100_regmos';
% [mosaicFile, saveFile] = trainNaturalScenes100_r2(pTrain);
% [mosaicFile, saveFile] = trainWhiteNoise(pTrain);

mosaicFile = '_mosaicAll_35336498' ;%['_' d1(dind).name(11:end-4)];% 
movieFile = 'pixium15_100/pix1_ns100_100Hz_mov'; 
spikesFile = 'pixium15_100/pix1_ns100_100Hz_sp';
% filterFile = ['pixium15_100/filters_pix1_nsBig15_1st_sv025_' mosaicFile]
% 
% mosaicFile = '_mosaicAll_1246640' ;%['_' d1(dind).name(11:end-4)];% 
% movieFile = 'ns100_r2_10/ns100_jan1_mov3'; 
% spikesFile = 'ns100_r2_10/ns100_jan1_sp3';
% filterFile = ['ns100_r2_10/filters_ns100_jan1_1st3_sh9_sv30_' mosaicFile];

% mosaicFile = '_mosaicAll_20116' ;%['_' d1(dind).name(11:end-4)];% 
% movieFile = 'ns100_r2_10_regmos/ns100_regmos_mov'; 
% spikesFile = 'ns100_r2_10_regmos/ns100_regmos_sp';
% filterFile = ['ns100_r2_10_regmos/filters_ns100_regmos_1st_sh9_sv30_' mosaicFile];
% % 
% pLoad.loadFile = ['ns100_r2_10_regmos/ns100_regmos'];
pLoad.loadFile = ['pixium15_sm/pix1_ns100_100Hz'];
pLoad.movieFile = movieFile;
pLoad.spikesFile = spikesFile;
pLoad.mosaicFile = mosaicFile;
% [movieFile, spikesFile] = loadSpikesAll(pLoad);
% % % % 
pRecon.movieFile = movieFile;
pRecon.spikesFile = spikesFile;
% pRecon.filterFile = filterFile;
% % % 
pRecon.mosaicFile = mosaicFile;

% evArr = [.2 .1 .3 .4];
evArr = [.01 .02 .03 .04 .1 .2];
trainFraction = [.16 .5 .66 .83 1];
% already did trainFraction =1
for evInd = 1:6
for trainFractionInd = 1:5
    [evInd trainFractionInd]
% if ~(evInd==1 && trainFractionInd==1)

% filterFile = ['ns100_r2_10/filters3_ns100_feb6_sh9_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];
% filterFile = ['pixium15_100/filters3_pix1_feb6_sh0_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];
filterFile = ['pixium15_sm/filters_pix_feb21_sh0_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];

pRecon.filterFile = filterFile;
pRecon.windowSize = 1;
pRecon.percentSV = evArr(evInd);
pRecon.trainFraction = trainFraction(trainFractionInd);
pRecon.shiftTime = 0;
[filterFile] = runReconstructSVD_fast_all(pRecon);
load(filterFile);
close;
figure; for fr = 1:25; subplot(5,5,fr); imagesc(reshape(filterMat(000+fr,:),[100 100]));colormap parula;  end;

% end
end
close all;
end
% % % 
% % % pTest.mosaicFile = mosaicFile;
% % % pTest.filterFile = filterFile;s
% % % testReconNS(pTest);
% % 
% load(filterFile);
% % figure; imagesc(reshape(sum(abs(filterMat),1),[96 96]));
% figure; for fr = 1:64; subplot(8,8,fr); imagesc(reshape(filterMat(1000+fr,:),[100 100]));colormap parula;  end;