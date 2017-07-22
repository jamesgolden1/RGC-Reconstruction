function obj = test(obj, varargin)
%TEST - test accuracy of RGC reconstruction on hallway movie



%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.7;% 3.2;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
nSteps = 500;%000;
nBlocks = 15;%30;

tic    %
rsFactor =1; stimSize = 100;
load([phospheneRootPath '/dat/stimuli/hallMovie.mat'])
szFrames = size(vidFrame,3);
hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
for ii = 1:szFrames
    hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
end

% Set hallway movie stimulus
testmovieshort = (255*ieScale(hallMovieResize)); clear hallMovieResize;

%%
%% Load image

if ~exist('spikeResp_hallway.mat','file')

clear coneParams

% One frame of a WN stimulus
% Set parameters for size

% coneParams.nSteps = nSteps;
% coneParams.row = 100; % should be set size to FOV
% coneParams.col = 100;
coneParams.fov = fov;
% % params.vfov = 0.7;
coneParams.startFrames = 0;

iStimNS = ieStimulusMovieCMosaic(testmovieshort(:,:,1:100),coneParams);
% iStimNS = ieStimulusMovieCMosaic(natScenes,coneParams);
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
rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta

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

% whiteNoiseSmall = natScenes;

end

%%
% % % % % % % % % Feb 6
% evArr = [.2 .1 .3 .4];
%
% % evArr = [.01 .02 .03 .04];
% trainFraction = [.16 .5 .66 .83 1];
%
% for evInd = 1:4
% for trainFractionInd = 1:5
% [evInd trainFractionInd]
% filterFile = ['ns100_r2_10/filters3_ns100_feb6_sh9_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];

% % % % % % % % % % % % Feb 22
% evArr = [.01 .02 .03 .04 .1 .2];
% trainFraction = [.16 .5 .66 .83 1];


mosaicFile = '_mosaic0';
% shiftArr = [4 ]; 
shiftArr = 2;
shiftTime = shiftArr; shiftind = 1;
windowSize = 1;
percentSVarr =.05%[ .75 .50 .25 ]%[.05 .075 .1 .125]% [.2 .4 .6 .8 1];
% percentSVarr = [.62 .38 .12]
trainSizeArr = .8%[.2 .4 .6 .8];

for percentSVind = 1%z:length(percentSVarr)
    for trainSizeInd = 1:length(trainSizeArr)
        
        percentSV = percentSVarr(percentSVind);
        shifttime = shiftArr(shiftind);
        trainSize = trainSizeArr(trainSizeInd);
%         filterFile = ['may22/filters/filters'  mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_tr%d',100*trainSize)];
        filterFile = ['june16prima/filters_wmean/filters'  mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_tr%d',100*trainSize)];
     
        load(filterFile);
%     shiftval = 9;
%     shiftval = shiftTime+9;
    
    spikeAug(1,:) = ones(1,size(spikeResp,2));
    spikeAug(2:9807,:) = spikeResp;
    % load('filters__mosaic0.mat')
    movRecon = filterMat'*spikeAug;
    movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     nFramesPlay = 40;
%     figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));
    % save('hallwayReconMovie.mat','movRecon');
    
    lambda = .01;
    filterMat2 = zeroFilter(filterMat,lambda);
    movRecon2 = filterMat2'*spikeAug;

%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     imSingle = testmovieshort(:,:,40); imMean = mean(imSingle(:));
%     mtest = max(abs(imSingle(:)-imMean));
%     imSingleRecon = movReconPlay(:,:,40); mtestrecon = max(abs(imSingleRecon(:)));
%     figure; subplot(131); imagesc((testmovieshort(:,:,40)-imMean)/mtest); colorbar; colormap gray;
%     
%     subplot(132); imagesc((movReconPlay(:,:,40)/mtestrecon)); colorbar; colormap gray
%     subplot(133); imagesc((movReconPlay(:,:,40)/mtestrecon)-(testmovieshort(:,:,40)-imMean)/mtest); colorbar; colormap gray

%     clear movRecon
    %
%     lambda = .01;
%     filterMat2 = zeroFilter(filterMat,lambda);
%     movRecon2 = filterMat2'*spikeAug;

    szLen = 550;
    shiftval = shiftArr+1;
%     shiftval = 5;
    % shiftval = 3;
    mc1 = (movReconPlay(:,:,shiftval+1:szLen+1));
    mc2 = (testmovieshort(:,:,1:szLen-shiftval+1));

%     mc1 = ieScale(movReconPlay(:,:,shiftval+1:szLen+1));
%     mc2 = ieScale(testmovieshort(:,:,1:szLen-shiftval+1));
%     clear bigMovie
%     bigMovie(1:100,101:200,:) = mc1;  bigMovie(1:100,1:100,:) = mc2;
%     errmov =(mc1-mean(mc1(:)))-(mc2-mean(mc2(:)));
%     errtot = ((errmov.^2));

    mc1rs = RGB2XWFormat(mc1);
    mc2rs = RGB2XWFormat(mc2);
    
    mc1rz = mc1rs;% - ones(size(mc1rs,1),1)*mean(mc1rs,1);    
    mc2rz = mc2rs;%%%% - ones(size(mc2rs,1),1)*mean(mc2rs,1);
    
%     std_recon = std(mc1rz(:));
%     std_test = std(mc2rz(:));
%     
%     mc1rz = mc1rz*std_test/std_recon;
    
    
    errmov =(mc1rz)-(mc2rz);
    errtot = ((errmov.^2));

    mse(percentSVind,trainSizeInd) = sqrt(mean(errtot(:)));
    mss(percentSVind,trainSizeInd) = sqrt(var(errtot(:)));
    
    trainSizeMat(percentSVind,trainSizeInd) = trainSizeArr(trainSizeInd);
    evMat(percentSVind,trainSizeInd) = percentSVarr(percentSVind);
    
    clear mc1 mc2 errmov errtot movieRecon
end
mse(percentSVind,:)
end


figure; 
plot(1e6*trainSizeMat',mse'/255,'-x','linewidth',4)
hold on;
% plot(1e6*trainSizeMat',mseh','-x','linewidth',4)
grid on; % axis([0 1e6 .13 .19]);
xlabel('Training Set Size','fontsize',15); 
ylabel('Mean Square Error - Reconstruction','fontsize',15);
set(gca,'fontsize',15)

legend('25% SVs','50% SVs','75% SVs','location','sw')
% legend('10% SVs','20% SVs','30% SVs','40% SVs');