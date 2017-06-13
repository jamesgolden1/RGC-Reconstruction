% function obj = movie(obj)

% 26*29+30*35+51*59+60*69
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
fov = 3.2;

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


% iStimNS = ieStimulusMovieCMosaic(rand(100,100,1),coneParams);
iStimNS = ieStimulusMovieCMosaic(testmovieshort(:,:,1:100),coneParams);
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

% save('spikeResp_hallway_high2.mat','spikeResp');

%%
% load('spikeResp_hallway_low100fr.mat')
load('spikeResp_hallway_high.mat')
spikeAug(1,:) = zeros(1,size(spikeResp,2));
% spikeAug(2:5177,:) = spikeResp;
spikeAug(2:8954,:) = spikeResp;
% load('filters__mosaic0.mat')
% load('filters_mosaic0_sv25_w1.mat')
% load('filters_mosaic0_sv15_w1_sh16.mat')
% load('filters_mosaic0_sv40_w1_sh16.mat');

mosaicFile = '_mosaic0';
windowSize = 1;
percentSV = .2;
shifttime = 16;
filterFile = ['filters576'  mosaicFile sprintf('_sv%2.2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime)];
load(filterFile);

movRecon = filterMat'*spikeAug;

movReconRS = reshape(movRecon,[100 100 616]);
% figure; ieMovie(reshape(movRecon2,[100 100 616]));
figure; ieMovie(255*ieScale(movReconRS(:,:,25:end-300)));
% save('hallwayReconMovie_sv40_sh16.mat','movRecon');

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

dp = sqrt(mgrd+mgcd);
dpCeil = 1.15*exp(-.0005*dp.^2); 
% dpCeil = 1.5*exp(-.0005*dp.^2); 
dpCeil(dpCeil>1)=1;
filterMat2 = filterMat.*dpCeil;
% filterMat2 = filterMat.*exp(-.003*dp.^2);
% filterMat2(dp>16) = 0;
filterMat2(1,:) = 1*filterMat(1,:);

movRecon2 = filterMat2'*spikeAug;
movReconRS2 = reshape(movRecon2,[100 100 616]);
% figure; ieMovie(reshape(movRecon2,[100 100 616]));
figure; ieMovie(255*ieScale(movReconRS2(:,:,25:end-300)));
% save('hallwayReconMovie2_sv40_sh16.mat','movRecon2');

%%

% p.vname = [reconstructionRootPath 'hallwayReconMovie_sv40_sh16.avi'];
% p.fps = 30;
% p.save = true;
% figure; ieMovie(movReconRS,p);