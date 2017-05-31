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
clear bpMosaic

cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
for cellTypeInd = 1:4
    clear bpParams
    bpParams.cellType = cellType{cellTypeInd};
    
    % FIX NEGATIVE AND POSITIVE HERE
    bpParams.ecc = patchEccentricity;
    bpParams.rectifyType = 1;
    bpMosaic{cellTypeInd} = bipolar(cMosaicNS, bpParams);
    bpMosaic{cellTypeInd}.set('sRFcenter',4);
    bpMosaic{cellTypeInd}.set('sRFsurround',0);
    bpMosaic{cellTypeInd}.compute(cMosaicNS);
end

%% RGC
clear params rgcParams
params.eyeRadius = patchEccentricity;
params.eyeAngle = 90;
innerRetina=ir(bpMosaic,params);
cellType = {'on parasol','off parasol','on midget','off midget'};

rgcParams.centerNoise = 0;
rgcParams.model = 'LNP';
%     rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta

rgcParams.type = cellType{1};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{2};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{3};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{4};
innerRetina.mosaicCreate(rgcParams);

innerRetina.compute(bpMosaic);


%%
toc
%% Look at covariance matrix
tic
spikesout  = RGB2XWFormat(mosaicGet(innerRetina.mosaic{1},'spikes'));
spikesout2 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{2},'spikes'));
spikesout3 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{3},'spikes'));
spikesout4 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{4},'spikes'));

timeBins = max([size(spikesout,2) size(spikesout2,2) size(spikesout3,2) size(spikesout4,2)]);

spikesoutsm = zeros(size(spikesout,1)+ size(spikesout2,1)+size(spikesout3,1)+size(spikesout4,1), timeBins,'uint8');
spikesoutsm(1:size(spikesout,1) ,1:size(spikesout,2) ) = spikesout;
spikesoutsm(size(spikesout,1)+[1:size(spikesout2,1)],1:size(spikesout2,2) ) = spikesout2;

spikesoutsm(size(spikesout,1)+size(spikesout2,1)+[1:size(spikesout3,1)] ,1:size(spikesout3,2) ) = spikesout3;
spikesoutsm(size(spikesout,1)+size(spikesout2,1)+size(spikesout3,1)+[1:size(spikesout4,1)] ,1:size(spikesout4,2) ) = spikesout4;

%     whiteNoiseSmall = natScenes;
clear  spikesout1 spikesout2 spikesout3 spikesout4 
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

% save('spikeResp_hallway_healthy2.mat','spikeResp');

% load(obj.filterFile);

clear spikesout spikesoutsm

%%

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


[mgr,mgc] = meshgrid(1:100,1:100);

[cmax,cind] = max(abs(filterMat),[],2);
[fmaxc,fmaxr] = ind2sub([100 100],cind);

mgrmat = mgr(:)*ones(1,size(fmaxr,1));
fmaxrmat = ones(size(mgrmat,1),1)*fmaxr';
mgrd = ((mgrmat - fmaxrmat)').^2;

mgcmat = mgc(:)*ones(1,size(fmaxc,1));
fmaxcmat = ones(size(mgcmat,1),1)*fmaxc';
mgcd = ((mgcmat - fmaxcmat)').^2;


% imagesc(signValWN(mind)*staim.*(.1+.9*reshape(exp(-.05*dp(1+paraIndPlus,:).^2),[100 100])));

dp = sqrt(mgrd+mgcd);
clear fmaxrmat fmaxcmat mgrmat mgrd mgcmat mgcd 
% filterMat2 = filterMat;
% filterMat2(dp>5) = 0;
expFilter = 1.2*(exp(-.03*dp.^2));
expFilter(expFilter>1) = 1;
filterMat2 = filterMat.*expFilter;

movRecon2 = filterMat2'*spikeAug;

movReconPlay2 = reshape(movRecon2,[100 100 size(spikeResp,2)]);
% nFramesPlay = 200;
figure; ieMovie(movReconPlay2(:,:,1:nFramesPlay));

clear movRecon2 dp expfilter