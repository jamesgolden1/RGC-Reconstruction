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
% pTrain.saveFile   = 'wn_Dec3_sp_2';
% % % [mosaicFile, saveFile] = trainNaturalScenes(pTrain);
% % [mosaicFile, saveFile] = trainWhiteNoise(pTrain);
% 
% mosaicFile = 'mosaicAll_117391';
% % saveFile = 'wnResponses\wn_Dec3_sp';
% saveFile = 'wnResponses\WNstim_response_stx2';
% movieFile = 'wn_Dec3_sp_mov_2';
% spikesFile = 'wn_Dec3_sp_spike_2';

% pLoad.loadFile = saveFile;

d1 = dir([reconstructionRootPath '\dat\pixiumMultmosaicAll_*.mat']);

for dind =  length(d1)%[1:6 8 10 12:length(d1)]
    mosaicFile = ['_' d1(dind).name(11:end-4)];    
    movieFile = 'pixiumMult\pix1_long_mov_nsBig_100hz';
    spikesFile = 'pixiumMult\pix1_long_sp_nsBig_100hz';
    filterFile = ['pixiumMult\pix1_long_filter_nsBig_100hz_4st_' mosaicFile];
    
    pLoad.loadFile = ['pixiumMult\pix1_nsBig_100Hz'];
    pLoad.mosaicFile = mosaicFile;
    pLoad.movieFile = movieFile;
    pLoad.spikesFile = spikesFile;
    
    [movieFile, spikesFile] = loadSpikesAll(pLoad);
    
    
    % filterFile = 'wn_Dec3_sp_filter_new1_drop';%'pixium\pix1filter_100Hz_Parasol';
    pRecon.movieFile = movieFile;
    pRecon.spikesFile = spikesFile;
    pRecon.filterFile = filterFile;
    pRecon.mosaicFile = mosaicFile;
    pRecon.windowSize = 4;
    pRecon.percentSV = 0.125;
    [filterFile] = runReconstructSVD_fast_all(pRecon);
    
    delete([reconstructionRootPath '\dat\' movieFile mosaicFile '.mat']);
    
    %
    % pTest.mosaicFile = mosaicFile;
    % pTest.filterFile = filterFile;
    % testReconNS(pTest);
    
end

load(filterFile);
figure; imagesc(reshape(sum(abs(filterMat),1),[96 96]));
figure; for fr = 1:64; subplot(8,8,fr); imagesc(reshape(filterMat(00+fr,:),[96 96])); caxis([-.003 .003]); colormap parula;  end;