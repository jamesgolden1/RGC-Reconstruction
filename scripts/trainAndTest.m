% trainAndTest
% 
% Run the entire pipeline for the RGC linear decoder. 
%   1. trainNaturalScenes - generate training data - show movies to RGC mosaics
%   2. loadSpikesAll      - rearrange spikes
%   3. runReconstruct     - learn decoding filters
%   4. testReconNS        - decode with learned filters and evaluate
%                               accuracy of reconstructions.
% 

pTrain.mosaicFile = 'mosaic_ns_Nov18';
pTrain.saveFile   = 'ns_Nov18';
[mosaicFile, saveFile] = trainNaturalScenes(pTrain);

pLoad.loadFile = saveFile;
[movieFile, spikesFile] = loadSpikesAll(pLoad);

pRecon.movieFile = movieFile;
pRecon.spikesFile = spikesFile;
[filterFile] = runReconstructSVD_fast_all(pRecon);

pTest.mosaicFile = pTrain.mosaicFile;
pTest.filterFile = filterFile;
testReconNS(pTest);
