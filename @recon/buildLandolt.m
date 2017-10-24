function obj = buildLandolt(obj, varargin)
%BUILD - builds training set for the recon object
% Run natural scenes through the RGC array for the big four RGC cell types.
%
% inputs:
%   mosaicFile - a string that is used to save the mosaic file
%   buildFile - a string that is used to store the spikes and the movie stim
%
% See also: trainAndTest.m

%%

p = inputParser;
p.addParameter('mosaicFile',[],@ischar);
p.addParameter('buildFile',[],@ischar);
p.addParameter('stimFile',[],@ischar);
p.addParameter('respFile',[],@ischar);
p.addParameter('blockIn',1,@isnumeric);
p.addParameter('stimTypeBuild','ns',@ischar);
p.addParameter('nTrials',100,@isnumeric);
p.addParameter('contrast',0.1,@isnumeric);

p.addParameter('pixelWidth',70,@isnumeric);
p.addParameter('currentDecay',2,@isnumeric);

p.addParameter('gratingSpFreq',25,@isnumeric);
p.addParameter('horizontalFlag',false,@isnumeric);
p.addParameter('folderNameTrain','',@ischar);
p.addParameter('folderNameTest','',@ischar);
p.addParameter('filterFile','',@ischar);
p.addParameter('windowSize','',@isnumeric);
p.addParameter('percentSV','',@isnumeric);
p.addParameter('shifttime','',@isnumeric);

p.addParameter('gap',0,@isnumeric);

p.addParameter('resizeFraction',1,@isnumeric);

p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;
stimTypeBuild = p.Results.stimTypeBuild;
nTrials = p.Results.nTrials;
contrast = p.Results.contrast;
horizontalFlag = p.Results.horizontalFlag;

gratingSpFreq = p.Results.gratingSpFreq;

pixelWidth = p.Results.pixelWidth;
currentDecay = p.Results.currentDecay;
folderNameTrain = p.Results.folderNameTrain;
folderNameTest = p.Results.folderNameTest;
filterFile = p.Results.filterFile;
windowSize = p.Results.windowSize;
percentSV = p.Results.percentSV;
shifttime  = p.Results.shifttime;
gap = p.Results.gap;
resizeFraction = p.Results.resizeFraction;
%%


% p = inputParser;
% p.addParameter('mosaicFile',[],@ischar);
% p.addParameter('buildFile',[],@ischar);
% p.addParameter('stimFile',[],@ischar);
% p.addParameter('respFile',[],@ischar);
% p.addParameter('blockIn',1,@isnumeric);
% p.addParameter('stimTypeBuild','ns',@ischar);
% p.addParameter('nTrials',100,@isnumeric);
% p.addParameter('contrast',0.1,@isnumeric);
% 
% p.addParameter('gratingSpFreq',25,@isnumeric);
% 
% p.addParameter('horizontalFlag',false,@islogical);
% 
% p.KeepUnmatched = true;
% p.parse(varargin{:});
% mosaicFile = p.Results.mosaicFile;
% buildFile = p.Results.buildFile;
% stimFile = p.Results.stimFile;
% respFile = p.Results.respFile;
% blockIn = p.Results.blockIn;
% stimTypeBuild = p.Results.stimTypeBuild;
% nTrials = p.Results.nTrials;
% contrast = p.Results.contrast;
% gratingSpFreq = p.Results.gratingSpFreq;
% 
% horizontalFlag = p.Results.horizontalFlag;
% if isempty(buildFile)
%     buildFile = ['primaLandolt_training_' num2str(round(cputime*100))];
% end

tic

%%

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.7;% 3.2/2;

nSteps = 16;

% orient0 = 'down';
% gapsize0 = gap;
% landoltStim = imgLandoltC('orientation',orient0,'gapsize',gapsize0);

%%


nSteps = 16;
orient0 = 'down'; 
gapsize0 = gap;
landoltStim = imgLandoltC('orientation',orient0,'gapsize',gapsize0);
landoltRed = zeros(100,100);

if resizeFraction == 1
    
landoltResize = (landoltStim(1:100,1:100));
else
landoltResize = imresize(landoltStim(1:100,1:100),resizeFraction*[100 100]);
end

sizeL = length([ceil(51-resizeFraction*50):ceil(50+resizeFraction*50)]);%,ceil(51-resizeFraction*50):ceil(50+resizeFraction*50)]);

landoltRed(ceil(51-resizeFraction*50):ceil(50+resizeFraction*50),...
    ceil(51-resizeFraction*50):ceil(50+resizeFraction*50)) = ...
    landoltResize(1:sizeL,1:sizeL);


landoltMovie = (repmat(landoltRed,[1 1 nSteps]));
landoltMovie = (landoltMovie - mean(landoltMovie(:)));
landoltMovie = landoltMovie./max(landoltMovie(:));

landoltMovieContrast = 127+127*(contrast)*(landoltMovie(1:100,1:100,:));%+0.5;

% landoltMovieContrast = (contrast)*(landoltMovie(1:100,1:100,:) - 0.5);%+0.5;
% landoltMovieContrast(end,end,:) = .5; landoltMovieContrast(1,end,:) = -.5;
landoltMovieContrast(:,:,1) = zeros(100,100);
% landoltMovieContrast = 128+254*landoltMovieContrast;

%% Reconstruction filter parameters


load(filterFile);


% lambda = .01;
% filterMat2 = zeroFilter(filterMat,lambda);
%% Load image
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
    coneParams.startFrames =0;
    coneParams.endFrames = 0;
    % % params.vfov = 0.7;
    
    
    load('rngSeedTraining.mat'); rng(s1);
    iStimNS = ieStimulusMovieCMosaic(landoltMovieContrast,coneParams);
    cMosaicNS = iStimNS.cMosaic;
    
    % Bipolar
    % Create a set of bipolar cell types in the bipolar mosaic
    
    clear bpL
    
    bpL = bipolarLayer(cMosaicNS);
    
    % Make each type of bipolar mosaic
    cellType = {'on diffuse','off diffuse','on midget','off midget','on sbc'};
    
    % Stride isn't influencing yet.s
    clear bpMosaicParams
    bpMosaicParams.rectifyType = 1;  % Experiment with this
    bpMosaicParams.spread  = 1;  % RF diameter w.r.t. input samples
    bpMosaicParams.stride  = 1;  % RF diameter w.r.t. input samples
    bpMosaicParams.spreadRatio  = 9;  % RF diameter w.r.t. input samples
    bpMosaicParams.ampCenter = 1;%1.3;%1.5 _2
    bpMosaicParams.ampSurround = .5;%1;%.5
    % Maybe we need a bipolarLayer.compute that performs this loop
    for ii = 1:length(cellType)
        bpL.mosaic{ii} = bipolarMosaic(cMosaicNS, cellType{ii}, bpMosaicParams);
        bpL.mosaic{ii}.compute();
    end
    
%     bpL.window;

    
    % RGC
    
    clear rgcL rgcParams
    
    % Create retina ganglion cell layer object
    rgcL = rgcLayer(bpL);
    
    % There are various parameters you could set.  We will write a script
    % illustrating these later.  We need a description.
    rgcParams.centerNoise = 0;
    rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
    % mosaicParams.axisVariance = .1;
    
    % 28*32+31*35+55*63+61*70
    onPdiameter = 9.4;
    diameters = [onPdiameter onPdiameter*.9 onPdiameter*.5 onPdiameter*.45];  % In microns.
    
    cellType = {'on parasol','off parasol','on midget','off midget'};
    for ii = 1:length(cellType)
        rgcParams.rfDiameter = diameters(ii);
        rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{ii},cellType{ii},rgcParams);
    end
    
    rgcL.set('numberTrials',1);
    
    % Compute the inner retina response and visualize
    
    % Every mosaic has its input and properties assigned so we should be able
%     % to just run through all of them.
%     rgcL.compute('bipolarScale',50,'bipolarContrast',1);
    
%     rng(90845);

for iTrial = 1:nTrials
    tic
    iTrial

    cMosaicNS.computeCurrent();
    for ii = 1:length(cellType)
%         bpL.mosaic{ii} = bipolarMosaic(cMosaicNS, cellType{ii}, bpMosaicParams);
        bpL.mosaic{ii}.compute();
    end
    rgcL.compute('bipolarScale',50,'bipolarContrast',1);
    
    %% Get spikes and make reconstruction
    
    spikeResp = mosaicSpikes(rgcL);
    % save('spikeResp_primaLandolt.mat','spikeResp');
    
    % load(filterFile);
    spikeAug(1,:) = ones(1,size(spikeResp,2));
    spikeAug(2:9717,:) = spikeResp;
    movRecon = filterMat'*spikeAug;
    
    movTrials(:,iTrial) = movRecon(:,11); % tested 11th frame w code below
    
    spikeTrials(:,iTrial) = spikeAug(:,11);
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
    % figure; ieMovie(movReconPlay);
    % save('hallwayReconMovie.mat','movRecon');
    
    % lambda = .05;
    % filterMat2 = zeroFilter(filterMat,lambda);
    % movRecon2 = filterMat2'*spikeAug;
    toc
    ph=1;
% end

if mod(iTrial,50)==0
    trialReconPlay = reshape(movTrials,[100 100 size(movTrials,2)]);
    
    % save(fullfile(reconstructionRootPath,'dat','july17gratings','healthyH',[sprintf('recon_%4d.mat',10000*contrast)]),'trialReconPlay');
    
    
    if gap < 0
        save(fullfile(reconstructionRootPath,'dat',folderNameTest,'cont',[sprintf('recon_freqH%4d_cont%6d_rs%2d.mat',100*gratingSpFreq,1e5*contrast,100*resizeFraction)]),'trialReconPlay','contrast','spikeTrials');
        %   save(fullfile(reconstructionRootPath,'dat',folderNameTest,[sprintf('primaRecon_freqH%4d_cont%4d.mat',100*gratingSpFreq,1000*contrast)]),'trialReconPlay');
        %
    else
        save(fullfile(reconstructionRootPath,'dat',folderNameTest,'cont',[sprintf('recon_freq%4d_cont%6d_rs%2d.mat',100*gratingSpFreq,1e5*contrast,100*resizeFraction)]),'trialReconPlay','contrast','spikeTrials');
        
    end
    
end

end
%% Find best frame
% gratingsMovieImage = landoltMovieContrast(:,:,2);
% movRecon0 = movRecon - ones(10000,1)*mean(movRecon);
% stdRecon0 = std(movRecon0);
% std0 = std(single(landoltMovieContrast(:)));
% movReconNorm = std0*movRecon0./(ones(10000,1)*stdRecon0);
% % movRecon0 = movRecon - ones(10000,1)*me an(movRecon);
% 
% errMov = (movReconNorm - (single(gratingsMovieImage(:))*ones(1,size(spikeResp,2))));
% 
% % errMov = (movReconNorm - RGB2XWFormat(gratingsMovieContrast));
% size(errMov)
% mse1 = mean(errMov.^2);
% figure; plot(mse1)
% [mv,mi] = min(mse1)