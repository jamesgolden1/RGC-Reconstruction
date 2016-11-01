
% Run the reconstruction algorithm to generate the filters

% Must run loadSpikesAll.m to build stimFile and spikesFile

%% run Linear reconstruction for on, off, and joint on/off training 
windowsize = 8;
disp('Loading spike responses...')

stimFileName = 'movie_spikeResp_all0';
spikesFileName = 'spikeResp_all0';
matfON = matfile([reconstructionRootPath '\dat\' spikesFileName]);

movielength = 1*240000;%size(stim,2);
disp(['Total Movie Length in Frames: ' num2str(movielength)]);

fileext = 'mosaic_all';
trainSizeArray = 1;%[.6/8:.6/8:.6]; 
trainInd = 1;
includedComponentsArray = 3952/8;

srON = matfON.spikeResp;
% srOFF = matfOFF.spikeResp;
% clear matfON matfOFF
spikeResp1 = vertcat(srON(:,1:240000));%, srOFF(:,1:240000), matfONP.spikeResp, matfOFFP.spikeResp);

clear matfOFF

for incInd = 1%:length(includedComponentsArray)
    filterMat = linearReconstructSVD_short_midgets_both(stimFileName,spikeResp1,fileext, windowsize,includedComponentsArray(incInd),trainSizeArray(trainInd));
end

figure; imagesc(reshape(filterMat(2+1*8,:),96,96));