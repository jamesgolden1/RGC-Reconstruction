function obj = buildLeaf(obj, varargin)
%BUILD - builds training set for the recon object
% Run natural scenes through the RGC array for the big four RGC cell types.
%
% inputs:
%   mosaicFile - a string that is used to save the mosaic file
%   buildFile - a string that is used to store the spikes and the movie stim
%
% See also: trainAndTest.m
%
p = inputParser;
p.addParameter('mosaicFile',[],@ischar);
p.addParameter('buildFile',[],@ischar);
p.addParameter('stimFile',[],@ischar);
p.addParameter('respFile',[],@ischar);
p.addParameter('blockIn',1,@isnumeric);
p.addParameter('stimTypeBuild','ns',@ischar);
%
% p.addParameter('pixelWidth',70,@isnumeric);
% p.addParameter('currentDecay',2,@isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;
stimTypeBuild = p.Results.stimTypeBuild;
% pixelWidth = p.Results.pixelWidth;
% currentDecay = p.Results.currentDecay;

if isempty(buildFile)
    buildFile = ['NS_training_' num2str(round(cputime*100))];
end

tic
%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.7;% 3.2;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
% nSteps = 10;%000;
nBlocks = 15;%30;


tic


% natScenes = 255*ieScale(loadHallStimulus(nSteps));

nSteps = 85;
load([fullfile(phospheneRootPath,'dat','stimuli','silver_small.mat')])
imMovie = repmat(im1,[1 1 nSteps]);
imMovie(:,:,1) = zeros(100,100);
% imMovie(:,:,nSteps+1:100) = zeros(100,100,length(nSteps+1:100));
natScenes = (imMovie);
%%
%% Load image
clear coneParams

% One frame of a WN stimulus
% Set parameters for size

% coneParams.nSteps = nSteps;
% coneParams.row = 100; % should be set size to FOV
% coneParams.col = 100;
coneParams.fov = fov;
coneParams.cmNoiseFlag = 'none';
coneParams.osNoiseFlag = 'none';
% % params.vfov = 0.7;


iStimNS = ieStimulusMovieCMosaic(natScenes,coneParams);
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
rgcL.compute('bipolarScale',50,'bipolarContrast',1);

spikeResp = mosaicSpikes(rgcL);

%% Get reconstruction filters

% % rd = RdtClient('isetbio');
% % rd.crp('/resources/data/istim');
% % filterFile = 'filtersmosaic0_sv50_w1_sh17_dr0.mat';
% % data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% % filterMat = data.filterMat; clear data;

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug27/filtersmosaic0_sv50_w1_sh15_dr0.mat')

rd = RdtClient('isetbio');
rd.crp('/resources/data/reconstruction');
filterFile = 'filtersmosaic0_sv50_w1_sh15_dr0_aug27.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterMat = data.filterMat; clear data;

% % % % % % % % % %
%     if isempty(pixelWidth)
%     filterFileFull  = fullfile(reconstructionRootPath,'/dat/', ...
%     [filterFile '.mat']);
%
%     else
%     filterFileFull  = fullfile(reconstructionRootPath,'/dat/', ...
%     [filterFile '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '.mat']);
%     end
%         load(filterFileFull);
% % % % % % %


%% Build reconstruction

%     shiftval = 9;
%     shiftval = shiftTime+9;

spikeAug(1,:) = ones(1,size(spikeResp,2));
spikeAug(2:9716+1,:) = spikeResp;
% load('filters__mosaic0.mat')

movRecon = filterMat'*spikeAug;
mse = 0;
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     nFramesPlay = 40;
%     figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));
% save('hallwayReconMovie.mat','movRecon');

lambda = .0025;%.01;%.0075;
filterMat2 = zeroFilter(filterMat,lambda);
movRecon2 = filterMat2'*spikeAug;

%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon2.avi']
%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon_pros18_filt0025.avi']
%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon_pros35_filt0005.avi']
%     p2.save = true;
%     p2.FrameRate = 25;
%
%     figure; ieMovie(reshape(movRecon2,[100 100 size(movRecon2,2)]),p2);
%     figure; ieMovie(reshape(movRecon2,[100 100 size(movRecon2,2)]));

figure; imagesc(reshape(movRecon2(:,26),[100 100])); axis image; axis off; colormap gray
clim1 = caxis;
set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_healthy.pdf']);


%% Build on reconstruction

%    28*32+31*35+55*63+61*70
spikeOn = zeros(size(spikeAug));
spikeOn(1:1+28*32,:) = spikeAug(1:1+28*32,:);
spikeOn(28*32+31*35+1:28*32+31*35+55*63+1,:) = spikeAug(28*32+31*35+1:28*32+31*35+55*63+1,:);

movReconOn = filterMat2'*spikeOn;

%     figure; ieMovie(reshape(movReconOn,[100 100 size(movRecon2,2)]));

figure; imagesc(reshape(movReconOn(:,26),[100 100])); axis image; axis off; colormap gray
caxis(clim1);
set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_healthy_on.pdf']);

%% Separate on reconstruction into positive and negative components

% filterMatNegInd = (filterMat2<0);
% filterMatNeg = zeros(size(filterMat2));
% filterMatNeg(filterMatNegInd) = filterMat2(filterMatNegInd);
% movReconOffNeg = filterMatNeg'*spikeOn;
% figure; imagesc(reshape(movReconOffNeg(:,26),[100 100])); axis image; axis off; colormap gray
% 
% filterMatPosInd = (filterMat2>0);
% filterMatPos = zeros(size(filterMat2));
% filterMatPos(filterMatPosInd) = filterMat2(filterMatPosInd);
% movReconOffPos = filterMatPos'*spikeOn;
% figure; imagesc(reshape(movReconOffPos(:,26),[100 100])); axis image; axis off; colormap gray
% 
% figure; imagesc(reshape(movReconOffNeg(:,26)+movReconOffPos(:,26),[100 100])); axis image; axis off; colormap gray
% 
% figure;
% subplot(131);
% imagesc(reshape(movReconOffPos(:,26),[100 100])); caxis([-120 120]); axis image; axis off; colormap gray
% 
% title('On Response, Center');
% subplot(132);
% 
% imagesc(reshape(movReconOffNeg(:,26),[100 100])); caxis([-120 120]); axis image; axis off; colormap gray
% title('On Response, Surround');
% subplot(133);
% imagesc(reshape(movReconOffNeg(:,26)+movReconOffPos(:,26),[100 100])); caxis([-45 45]); axis image; axis off; colormap gray

%% Build off reconstruction

spikeOff = zeros(size(spikeAug));

spikeOff(1+28*32:1+28*32+31*35,:) = spikeAug(1+28*32:1+28*32+31*35,:);
spikeOff(1+28*32+31*35+55*63:end,:) = spikeAug(1+28*32+31*35+55*63:end,:);

movReconOff = filterMat2'*spikeOff;

%     figure; ieMovie(reshape(movReconOff,[100 100 size(movRecon2,2)]));

figure; imagesc(reshape(movReconOff(:,26),[100 100])); axis image; axis off; colormap gray
caxis(clim1);
set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_healthy_off.pdf']);

%% Separate off reconstruction into positive and negative components
% filterMatNegInd = (filterMat2<0);
% filterMatNeg = zeros(size(filterMat2));
% filterMatNeg(filterMatNegInd) = filterMat2(filterMatNegInd);
% movReconOffNeg = filterMatNeg'*spikeOff;
% figure; imagesc(reshape(movReconOffNeg(:,26),[100 100])); axis image; axis off; colormap gray
% 
% filterMatPosInd = (filterMat2>0);
% filterMatPos = zeros(size(filterMat2));
% filterMatPos(filterMatPosInd) = filterMat2(filterMatPosInd);
% movReconOffPos = filterMatPos'*spikeOff;
% figure; imagesc(reshape(movReconOffPos(:,26),[100 100])); axis image; axis off; colormap gray
% 
% figure; imagesc(reshape(movReconOffNeg(:,26)+movReconOffPos(:,26),[100 100])); axis image; axis off; colormap gray
% 
% figure;
% subplot(131);
% imagesc(reshape(movReconOffNeg(:,26),[100 100])); axis image; axis off; colormap gray
% title('Off Response, Center');
% subplot(132);
% imagesc(reshape(movReconOffPos(:,26),[100 100])); axis image; axis off; colormap gray
% title('Off Response, Surround');
% subplot(133);
% imagesc(reshape(movReconOffNeg(:,26)+movReconOffPos(:,26),[100 100])); axis image; axis off; colormap gray

%% Build activation maps for each type of RGC

for typeInd = 1:4
    rgcL.mosaic{typeInd}.plot('mosaicFill'); title(''); colorbar off; axis image; axis off;
    set(gcf,'PaperPositionMode','auto')
    % %            print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures/spatial_mos' num2str(typeInd) '.pdf']);
    %             print(gcf, '-dsvg',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_activation_mos' num2str(typeInd) '.svg']);
    
end

