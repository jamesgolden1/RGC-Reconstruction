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

if 1%~exist('spikeResp_hallway.mat','file')

clear coneParams

% One frame of a WN stimulus
% Set parameters for size
    coneParams.cmNoiseFlag = 'random';
    coneParams.osNoiseFlag = 'random';
    % % params.vfov = 0.7;
    
    
    iStimNS = ieStimulusMovieCMosaic(testmovieshort(1:55),coneParams);
    cMosaicNS = iStimNS.cMosaic;
    
    %% Bipolar
    %% Create a set of bipolar cell types in the bipolar mosaic
    
    clear bpL
    
    bpL = bipolarLayer(cMosaicNS);
    
    % Make each type of bipolar mosaic
    cellType = {'on diffuse','off diffuse','on midget','off midget','on sbc'};
    
    % Stride isn't influencing yet.s
    clear bpMosaicParams
    bpMosaicParams.rectifyType = 1;  % Experiment with this
    bpMosaicParams.spread  = 1;  % RF diameter w.r.t. input samples
    bpMosaicParams.stride  = 1;  % RF diameter w.r.t. input samples
    bpMosaicParams.spreadRatio  = 10;  % RF diameter w.r.t. input samples
    bpMosaicParams.ampCenter = 1.3;%1.5 _2
    bpMosaicParams.ampSurround = 1;%.5
    % Maybe we need a bipolarLayer.compute that performs this loop
    for ii = 1:length(cellType)
        bpL.mosaic{ii} = bipolarMosaic(cMosaicNS, cellType{ii}, bpMosaicParams);
        bpL.mosaic{ii}.compute();
    end
    
%     bpL.window;

    
    %% RGC
    
    clear rgcL rgcParams
    
    % Create retina ganglion cell layer object
    rgcL = rgcLayer(bpL);
    
    % There are various parameters you could set.  We will write a script
    % illustrating these later.  We need a description.
    rgcParams.centerNoise = 0;
    rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
    % mosaicParams.axisVariance = .1;
    
    % 27*31+31*35+54*62+63*72
    onPdiameter = 9.4;
    diameters = [onPdiameter onPdiameter*.9 onPdiameter*.5 onPdiameter*.45];  % In microns.
    
    cellType = {'on parasol','off parasol','on midget','off midget'};
    for ii = 1:length(cellType)
        rgcParams.rfDiameter = diameters(ii);
        rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{ii},cellType{ii},rgcParams);
    end
    
    nTrials = 1; rgcL.set('numberTrials',nTrials);
    
    %% Compute the inner retina response and visualize
    
    % Every mosaic has its input and properties assigned so we should be able
    % to just run through all of them.
    rgcL = rgcL.compute('bipolarScale',50,'bipolarContrast',1);


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