function obj = buildPrimaLeaf(obj, varargin)
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
p.addParameter('pixelWidth',70,@isnumeric);
p.addParameter('currentDecay',2,@isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;
stimTypeBuild = p.Results.stimTypeBuild;
pixelWidth = p.Results.pixelWidth;
currentDecay = p.Results.currentDecay;

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

%% Load image stimulus

nSteps = 85;
load([fullfile(phospheneRootPath,'dat','stimuli','silver_small.mat')]);
% im1 = 255*ieScale( imresize(IMAGESr(:,:,1),[100 100]));
% load('/Users/james/Documents/MATLAB/attractor/IMAGES_RAW.mat')
% im1 = 255*ieScale( imresize(IMAGESr(:,:,1),[100 100]));
% load('wheel.mat'); im1 = reshape(im1,[100 100]);
imMovie = repmat(im1,[1 1 nSteps]);
imMovie(:,:,1) = zeros(100,100);
% imMovie(:,:,nSteps+1:100) = zeros(100,100,length(nSteps+1:100));
movieIn = (imMovie);

%% Compute response

primaParams.pixelWidth = pixelWidth*1e-6; % meters
primaParams.ecc = 1.8;       % mm
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance

primaParams.currentDecay = currentDecay;%2.6;%

primaRecon = primaArray(movieIn,primaParams);

primaRecon.compute(movieIn)
%     innerRetina = primaRecon.innerRetina;

spikeResp = mosaicSpikes(primaRecon.innerRetina);

%% Load reconstruction filters

%     load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug27/filtersmosaic0_sv50_w1_sh15_dr0.mat')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug27prima70/filtersmosaic0_sv 5_w1_sh4_dr0_pitch_70_decay_2_aug27.mat')

rd = RdtClient('isetbio');
rd.crp('/resources/data/reconstruction');
filterFile = 'filtersmosaic0_sv5_w1_sh4_dr0_pitch_70_decay_2_aug27.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterMat = data.filterMat; clear data;

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0.mat');
% filterNS = filterMat; clear filterMat
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv 5_w1_sh4_dr0_pitch_70_decay_2.mat')
%     load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug27prima35/filtersmosaic0_sv10_w1_sh5_dr0_pitch_35_decay_2.mat');

% rd = RdtClient('isetbio');
% rd.crp('/resources/data/istim');
% filterFile = 'filtersmosaic0_sv50_w1_sh17_dr0.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% filterMat = data.filterMat; clear data;

% rd = RdtClient('isetbio');
% rd.crp('/resources/data/reconstruction');
% filterFile = 'filtersmosaic0_sv50_w1_sh15_dr0_aug27.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% filterMat = data.filterMat; clear data;

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

%% Compute learning reconstruction

%     load('ws_prima_leaf_recon_sep20.mat');
%     load('filtersmosaic0_sv 2_w1_sh3_dr0.mat')
%        load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug27prima70/filtersmosaic0_sv 5_w1_sh4_dr0_pitch_70_decay_2_aug27.mat')

%     lambda = .0005;%.0025;%.01;%.0075;
% lambda = .005;
%     lambda = .0025;
% %     lambda = .001;
%     filterMat2 = zeroFilter(filterMat,lambda);
movRecon2 = filterMat'*spikeAug;

%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon2.avi']
%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon_pros18_filt0025.avi']
%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon_pros35_filt0005.avi']
%     p2.save = true;
%     p2.FrameRate = 25;
%
%     figure; ieMovie(reshape(movRecon2,[100 100 size(movRecon2,2)]),p2);
%     figure; ieMovie(reshape(movRecon2,[100 100 size(movRecon2,2)]));

figure; imagesc(reshape(movRecon2(:,3),[100 100])); axis image; axis off; colormap gray
clim1 = caxis;
set(gcf,'PaperPositionMode','auto')
% %     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_healthy.pdf']);
% %     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_prima70.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_sep20_prima35_decay2_bigSpread_learning.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_sep20_prima70_decay2_bigSpread_learning.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/wheel_recon_sep20_prima35_decay2_bigSpread_learning.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/wheel_recon_sep20_prima70_decay2_bigSpread_learning.pdf']);
%% 3, 8, 14
% mi = mi+1;
%      imagesc(reshape(movRecon2(:,mi),[100 100])'); axis image; axis off; colormap gray
% title(sprintf('%d',mi'));

%% Compute no-learning reconstruction

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug27/filtersmosaic0_sv50_w1_sh15_dr0.mat')

rd = RdtClient('isetbio');
rd.crp('/resources/data/reconstruction');
filterFile = 'filtersmosaic0_sv50_w1_sh15_dr0_aug27.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterMat = data.filterMat; clear data;

lambda = .0025;
%     lambda = .001;
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

figure; imagesc(reshape(movRecon2(:,9),[100 100])); axis image; axis off; colormap gray
clim1 = caxis;
set(gcf,'PaperPositionMode','auto')
% %     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_healthy.pdf']);
% %     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_prima70.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_sep20_prima70_decay2_bigSpread_noLearning.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_sep20_prima35_decay2_bigSpread_noLearning.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/wheel_recon_sep20_prima70_decay2_bigSpread_noLearning.pdf']);
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/wheel_recon_sep20_prima35_decay2_bigSpread_noLearning.pdf']);

%% Compute no-learning only on reconstruction

%    28*32+31*35+55*63+61*70
spikeOn = zeros(size(spikeAug));
spikeOn(1:1+28*32,:) = spikeAug(1:1+28*32,:);
spikeOn(28*32+31*35+1:28*32+31*35+55*63+1,:) = spikeAug(28*32+31*35+1:28*32+31*35+55*63+1,:);

movReconOn = filterMat2'*spikeOn;

%     figure; ieMovie(reshape(movReconOn,[100 100 size(movRecon2,2)]));
figure; imagesc(reshape(movReconOn(:,9),[100 100])); axis image; axis off; colormap gray
%     figure; imagesc(reshape(movReconOn(:,26)+movReconOff(:,26),[100 100])); axis image; axis off; colormap gray
%     caxis(clim1);
set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_healthy_on.pdf']);
%          print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_sep20_prima35_decay2_bigSpread_NoLearning_on.pdf']);
%          print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/wheel_recon_sep20_prima70_decay2_bigSpread_NoLearning_on.pdf']);
%          print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/wheel_recon_sep20_prima35_decay2_bigSpread_NoLearning_on.pdf']);
% figure; imagesc(imresize(im1,[15 15])); colormap gray; axis image axis off;
% figure; imagesc(imresize(im1,[8 8])); colormap gray; axis image; axis off;
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_prima_resize70.pdf']);
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_prima_on_learning.pdf']);
%      print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_prima_on.pdf']);


%% Separate no-learning, only on reconstruction into positive and negative components

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

%% Compute no-learning only off reconstruction

spikeOff = zeros(size(spikeAug));

spikeOff(1+28*32:1+28*32+31*35,:) = spikeAug(1+28*32:1+28*32+31*35,:);
spikeOff(1+28*32+31*35+55*63:end,:) = spikeAug(1+28*32+31*35+55*63:end,:);

movReconOff = filterMat2'*spikeOff;

%     figure; ieMovie(reshape(movReconOff,[100 100 size(movRecon2,2)]));

figure; imagesc(reshape(movReconOff(:,9),[100 100])); axis image; axis off; colormap gray
%     caxis(clim1);
set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_healthy_off.pdf']);
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_prima_off.pdf']);
%     print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_prima_off_learning.pdf']);
%         print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_recon_sep20_prima70_decay2_bigSpread_NoLearning_off.pdf']);
ph=1;

%% Separate no-learning, only off reconstruction into positive and negative components
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
%%

for typeInd =1:4
    primaRecon.innerRetina.mosaic{typeInd}.plot('mosaicFill'); title(''); colorbar off; axis image; axis off;
    %            rgcL.mosaic{typeInd}.plot('mosaicFill'); title(''); colorbar off; axis image; axis off;
    set(gcf,'PaperPositionMode','auto')
    % %            print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures/spatial_mos' num2str(typeInd) '.pdf']);
    %             print(gcf, '-dsvg',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/leaf_activation_sep20_prima_decay2_bigSpread_mos' num2str(typeInd) '.svg']);
    
end

