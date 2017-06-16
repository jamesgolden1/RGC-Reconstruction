function obj = movie(obj)
% Generate the reconstruced hallway movie from the recon object.
% 
% 
%%
rsFactor = 1; stimSize = 100;

load([phospheneRootPath '/dat/stimuli/hallMovie.mat'])
szFrames = size(vidFrame,3);
hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
for ii = 1:szFrames
    hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
end

% Set hallway movie stimulus
testmovieshort = (255*ieScale(hallMovieResize)); clear hallMovieResize;


%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.7;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
nSteps = 2000;%12000;%000;
nBlocks = 15;%30;


blockNum =1%1%:nBlocks
tic

blockNum
%     load(filenameRGC, 'innerRetina');

%%
%% Load image
clear coneParams

% One frame of a WN stimulus
% Set parameters for size

% coneParams.nSteps = nSteps;
% coneParams.row = 100; % should be set size to FOV
% coneParams.col = 100;
coneParams.fov = fov;
% % params.vfov = 0.7;

coneParams.startFrames = 0;
    

% iStimNS = ieStimulusMovieCMosaic(rand(100,100,1),coneParams);
iStimNS = ieStimulusMovieCMosaic(testmovieshort(:,:,1:40),coneParams);
cMosaicNS = iStimNS.cMosaic;

%% Bipolar
% clear bpMosaic
% 
% cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
% for cellTypeInd = 1:4
%     clear bpParams
%     bpParams.cellType = cellType{cellTypeInd};
%     
%     % FIX NEGATIVE AND POSITIVE HERE
%     bpParams.ecc = patchEccentricity;
%     bpParams.rectifyType = 1;
%     bpMosaic{cellTypeInd} = bipolar(cMosaicNS, bpParams);
%     bpMosaic{cellTypeInd}.set('sRFcenter',4);
%     bpMosaic{cellTypeInd}.set('sRFsurround',0);
%     bpMosaic{cellTypeInd}.compute(cMosaicNS);
% end

% Make each type of bipolar mosaic
cellType = {'ondiffuse','offdiffuse','onmidget','offmidget'};
clear bpMosaicParams
bpMosaicParams.rectifyType = 1;

bpMosaic  = cell(1,length(cellType));
bpNTrials = cell(1,length(cellType));
for ii = 1:length(cellType)
    
    bpMosaicParams.cellType = cellType{ii};
    
    bpMosaic{ii} = bipolarMosaic(cMosaic, bpMosaicParams);
    bpMosaic{ii}.set('sRFcenter',1);
    bpMosaic{ii}.set('sRFsurround',0);
    
    [~, bpNTrialsCenterTemp, bpNTrialsSurroundTemp] = ...
        bpMosaic{ii}.compute(cMosaic,'coneTrials',alignedC);
    bpNTrials{ii} = bpNTrialsCenterTemp - bpNTrialsSurroundTemp;
    
end
bpL.mosaic = bpMosaic;


% %% RGC
% clear params rgcParams
% params.eyeRadius = patchEccentricity;
% params.eyeAngle = 90;
% innerRetina=ir(bpMosaic,params);
% cellType = {'on parasol','off parasol','on midget','off midget'};
% 
% rgcParams.centerNoise = 0;
% rgcParams.model = 'LNP';
% %     rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% 
% rgcParams.type = cellType{1};
% innerRetina.mosaicCreate(rgcParams);
% rgcParams.type = cellType{2};
% innerRetina.mosaicCreate(rgcParams);
% rgcParams.type = cellType{3};
% innerRetina.mosaicCreate(rgcParams);
% rgcParams.type = cellType{4};
% innerRetina.mosaicCreate(rgcParams);
% 
% innerRetina.compute(bpMosaic);

clear rgcLayer irParams mosaicParams

% Create retina ganglion cell layer object
rgcL = rgcLayer(bpL);

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
mosaicParams.centerNoise = 0;
mosaicParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% mosaicParams.axisVariance = .1;
mosaicParams.type  = cellType;
mosaicParams.model = 'GLM';

% diameters = [5 5 3 3 10];  % In microns.

cellType = {'on parasol','off parasol','on midget','off midget'};
for ii = 1:length(cellType)
%     mosaicParams.rfDiameter = diameters(ii);
    mosaicParams.type = cellType{ii};
    mosaicParams.inMosaic = 1;   % Could switch up and match inputs to outputs
    rgcL.mosaicCreate(mosaicParams);
end

nTrials = 1; rgcL.set('numberTrials',nTrials);


%%
toc
%% Get spikes to make reconstruction movie

spikeResp = mosaicSpikes(innerRetina);

%%

spikesout = double(spikesoutsm);
pointer = 0;%(blockNum-1)*blocklength;

for i = 1:size(spikesoutsm,2)/10
    blocksize = 10;
    endval = i*blocksize;
    if endval > size(spikesout,2)
        endval = size(spikesout,2);
    end
    startval = (i-1)*blocksize + 1;
    spikeResp(:,pointer+i) = sum(spikesout(:,startval:endval),2);
end

% save('spikeResp_hallway.mat','spikeResp');
% filterMat = obj.readFilter()

rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
filterFile = 'filters_mosaic0_sv80_w1_sh4_may22.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterMat = data.filterMat; clear data;

%%
spikeAug(1,:) = ones(1,size(spikeResp,2));
spikeAug(2:9807,:) = spikeResp;
% load('filters__mosaic0.mat')
movRecon = filterMat'*spikeAug;
movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
nFramesPlay = 40;
figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));
% save('hallwayReconMovie.mat','movRecon');

clear movRecon
%%
lambda = .01;
filterMat2 = zeroFilter(filterMat,lambda);
movRecon2 = filterMat2'*spikeAug;

movReconPlay2 = reshape(movRecon2,[100 100 size(spikeResp,2)]);
% nFramesPlay = 200;
figure; ieMovie(movReconPlay2(:,:,1:nFramesPlay));

clear movRecon2 
