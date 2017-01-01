% trainAndTest
% 
% Run the entire pipeline for the RGC linear decoder. 
%   1. trainNaturalScenes - generate training data - show movies to RGC mosaics
%   2. loadSpikesAll      - rearrange spikes
%   3. runReconstruct     - learn decoding filters
%   4. testReconNS        - decode with learned filters and evaluate
%                               accuracy of reconstructions.
% 

pTrain.mosaicFile = 'mosaic_wn_Dec3_sp_2';
pTrain.saveFile   = 'ns100_Dec31';
[mosaicFile, saveFile] = trainNaturalScenes100_r2(pTrain);
% % [mosaicFile, saveFile] = trainWhiteNoise(pTrain);
% 
% mosaicFile = 'mosaicAll_117391';
% % saveFile = 'wnResponses\wn_Dec3_sp';
% saveFile = 'wnResponses\WNstim_response_stx2';
% movieFile = 'wn_Dec3_sp_mov_2';
% spikesFile = 'wn_Dec3_sp_spike_2';

% % pLoad.loadFile = saveFile;
% pLoad.loadFile = 'pixium\pix1_nsBig_100Hz';
% movieFile = 'pixium\pix1_long_mov_nsBig_100hz';
% spikesFile = 'pixium\pix1_long_sp_nsBig_100hz';
% filterFile = 'pixium\pix1_long_filter_nsBig_100hz_4st';

% pLoad.loadFile = 'pixiumBig\pix1_nsBig_100Hz';
% movieFile = 'pixiumBig\pix1_long_mov_nsBig_100hz';
% spikesFile = 'pixiumBig\pix1_long_sp_nsBig_100hz';
% filterFile = 'pixiumBig\pix1_long_filter_nsBig_100hz_4st';


% % mosaicFile = '_mosaicAll_25760187' ;%['_' d1(dind).name(11:end-4)];
% % movieFile = 'ns200/ns_Dec22_mov';
% % spikesFile = 'ns200/ns_Dec22_sp';
% % filterFile = ['ns200/filters_nsDec22_4st_sv005_' mosaicFile];

% mosaicFile = '_mosaicAll_8372855' ;%['_' d1(dind).name(11:end-4)];
% movieFile = 'pixium25/pix1_nsBig_100Hz_mov';
% spikesFile = 'pixium25/pix1_nsBig_100Hz_sp';
% filterFile = ['pixium25/pix1_nsBig_100Hz_1st_sv05_' mosaicFile];
% 
% pLoad.loadFile = ['pixium25/pix1_nsBig_100Hz'];
%     

mosaicFile = '_mosaicAll_19261' ;%['_' d1(dind).name(11:end-4)];% 
movieFile = 'ns100/ns_Dec22_mov';
spikesFile = 'ns100/ns_Dec22_sp';
filterFile = ['ns100/filters_nsDec22_1st_sv05_' mosaicFile];

pLoad.loadFile = ['ns100/ns100_Dec22'];
pLoad.movieFile = movieFile;
pLoad.spikesFile = spikesFile;
pLoad.mosaicFile = mosaicFile;
% [movieFile, spikesFile] = loadSpikesAll(pLoad);
% % 
pRecon.movieFile = movieFile;
pRecon.spikesFile = spikesFile;
pRecon.filterFile = filterFile;
% % 
pRecon.mosaicFile = mosaicFile;
pRecon.windowSize = 1;
pRecon.percentSV = 0.05;
% [filterFile] = runReconstructSVD_fast_all(pRecon);
% 
% pTest.mosaicFile = mosaicFile;
% pTest.filterFile = filterFile;s
% testReconNS(pTest);

% load(filterFile);
% % figure; imagesc(reshape(sum(abs(filterMat),1),[96 96]));
% figure; for fr = 1:64; subplot(8,8,fr); imagesc(reshape(filterMat(00+fr,:),[200 200])); caxis([-.003 .003]); colormap parula;  end;