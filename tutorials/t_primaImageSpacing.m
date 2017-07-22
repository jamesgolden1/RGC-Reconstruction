% t_primaImageSpacing
% 
% The Prima object allows for the simulation of a subretinal prosthesis.
% The activation of the photovoltaic electrode array in response to a movie
% stimulus is generated. The bipolar current and RGC spikes are generated
% in response to the electrode stimulation. 
% 
% With the RGC spikes, and a retina with a fixed set of parameters, we can
% generate the stimulus reconstruction to examine how it is affected by the
% prosthetic stimulation.
% 
% 
% TOOLBOX DEPENDENCIES - these must be downloaded and added to the matlab
%                               path with subfolders.
%       isetbio:            http://github.com/isetbio/isetbio [bipolar branch]
%       RGC-Reconstruction: https://github.com/Chichilnisky-Lab/RGC-Reconstruction
%       EJLPhosphene:       https://github.com/isetbio/EJLPhosphene
%       RemoteDataToolbox:  https://github.com/isetbio/RemoteDataToolbox
% 

%% Load stimulus movie

% stimFrames = 50;
% movieIn = loadHallStimulus(stimFrames);

% nSteps = 20;
% load([fullfile(phospheneRootPath,'dat','stimuli','silver_small.mat')])
% imMovie = repmat(im1,[1 1 nSteps]);
% imMovie(:,:,1) = zeros(100,100);
% movieIn = imMovie;

%%
blockNum = 1; nStepsIn = 200;
if blockNum <= 288
    movsm = parload(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(mod(blockNum-1,12)+1) '.mat']);
    
    natScenes = movsm(1:100,1:100,nStepsIn*floor((blockNum-1)/12)+randperm(nStepsIn));
else
    movsm = parload(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(12+mod(blockNum-1,12)+1) '.mat']);
    
    natScenes = movsm(1:100,1:100,nStepsIn*(floor((-288+blockNum-1)/12))+randperm(nStepsIn));
end
clear movsm
%%

filterFile = ('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/may22/filters_wmean/filters_mosaic0_sv25_w1_sh4_tr80.mat');
load(filterFile);


% %% Spatial zeroing of filter
lambda = .01;
filterMat2 = zeroFilter(filterMat,lambda);

%%
ii=1;
nSteps = 16;
imMovie = repmat(natScenes(:,:,ii),[1 1 nSteps]);
imMovie(:,:,1) = zeros(100,100);
movieIn = imMovie;
    
%% Simulate bipolar and RGC response to prosthesis stimulation

primaParams.pixelWidth = 1*35e-6; % meters
primaParams.ecc = 1.8;       % deg
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance 

primaRecon = primaArray(movieIn,primaParams);

for ii = 1:200
ii
tic
    imMovie = repmat(natScenes(:,:,ii),[1 1 nSteps]);
imMovie(:,:,1) = zeros(100,100);
movieIn = imMovie;

primaRecon.compute(movieIn)

%% Reconstruct - get spikes and decoding filter
 
spikeResp = mosaicSpikes(primaRecon.innerRetina);
% save('spikeResp_hallway.mat','spikeResp');

% % Remote data toolbox - download decoding filter
% rd = RdtClient('isetbio');
% rd.crp('/resources/data/istim');
% % filterFile = 'filters_mosaic0_sv75_w1_sh2_may26primaSmall';
% filterFile = 'filters_mosaic0_sv20_w1_sh2_dr0';
% % filterFile = 'filters_mosaic0_sv10_w1_sh2_dr0';
% data  = rd.readArtifact(filterFile, 'type', 'mat');
% filterMat = data.filterMat; clear data;
% 
% filterFile = 'filters_mosaic0_sv80_w1_sh4_may22.mat';
% load(filterFile);

% % %% Generate reconstructed movie
spikeAug(1,:) = ones(1,size(spikeResp,2));
spikeAug(1+[1:size(spikeResp,1)],:) = spikeResp;
% movRecon = filterMat'*spikeAug;
% movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
% nFramesPlay = size(spikeResp,2);
% % figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));
% 
% 
% movRecon2 = filterMat2'*spikeAug;
% movReconPlay2 = reshape(movRecon2,[100 100 size(spikeResp,2)]);
% % figure; ieMovie(movReconPlay2(:,:,1:nFramesPlay));


%%


%  27*31+31*35+54*62+63*72
onParasolInd = 1+27*31;
offParasolInd = onParasolInd+31*35;
onMidgetInd = offParasolInd+54*62;
%  load('filters_mosaic0_sv50_w1_sh4_tr80.mat')
spikeNoLearn = spikeAug;
spikeNoLearn(onParasolInd+1:offParasolInd,:) = 0;
spikeNoLearn(onMidgetInd+1:end,:) = 0;
 
% spikeNoLearn(1:onParasolInd,:) = 0;
% spikeNoLearn(offParasolInd+1:onMidgetInd,:) = 0;
movReconNoLearn = filterMat'*spikeNoLearn;
movReconPlayNoLearn = reshape(movReconNoLearn,[100 100 size(spikeResp,2)]);
movReconPlayNoLearn = permute(movReconPlayNoLearn,[2 1 3]);
nFramesPlay = size(spikeResp,2);
% figure; ieMovie(movReconPlayNoLearn(:,:,1:nFramesPlay),'frameRate',10);
% figure; imagesc(movReconPlayNoLearn(:,:,9)); colormap gray; axis image; axis off;

movReconOut(:,:,ii) = movReconPlayNoLearn(:,:,9);
%%
% lambda = .01;
% filterMat2 = zeroFilter(filterMat,lambda);
movReconNoLearn2 = filterMat2'*spikeNoLearn;
movReconPlayNoLearn2 = reshape(movReconNoLearn2,[100 100 size(spikeResp,2)]);

movReconPlayNoLearn2 = permute(movReconPlayNoLearn2,[2 1 3]);
% figure; ieMovie(movReconPlayNoLearn2(:,:,1:nFramesPlay));

movReconOut2(:,:,ii) = movReconPlayNoLearn2(:,:,9);

toc


end

save('primaImageSpacign_pitch1_decay1.mat','movReconOut','movReconOut2','natScenes');