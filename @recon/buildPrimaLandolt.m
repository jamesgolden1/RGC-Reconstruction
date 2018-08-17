function obj = buildPrimaLandolt(obj, varargin)
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

p.addParameter('rngInput',1,@isnumeric);
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
rngInput = p.Results.rngInput;
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

if isempty(buildFile)
    buildFile = ['primaLandolt_training_' num2str(round(cputime*100))];
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
nSteps = 16;%000;
% nBlocks = 15;%30;


%%


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

%% Reconstruction filter parameters

% mosaicFile = 'mosaic0';
% % windowSize = 1;
% % percentSV = .05;%.5;%.5;%.25;%.12;
% % % shifttime = 2;
% % shifttime = 4;%15;%3;%15;
% dropout = 0;
% 
% filterFile  = fullfile(folderNameTrain,...    
%     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);

filterMat = load(filterFile);

% lambda = .005;
% filterMat2 = zeroFilter(filterMat,lambda);
% clear filterMat;

%% Load image

primaParams.pixelWidth = pixelWidth*1e-6; % meters
primaParams.ecc = 1.8;       % mm
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance


primaRecon = primaArray(landoltMovieContrast,primaParams);

primaRecon.compute(landoltMovieContrast);
    rng(rngInput);
for iTrial = 1:nTrials
    iTrial
    tic
    % primaRecon = primaArray(landoltMovieContrast,primaParams);
    
    primaRecon.computeRGC();
    innerRetina = primaRecon.innerRetina;
    
    %% Get spikes and make reconstruction
    
    spikeResp = mosaicSpikes(primaRecon.innerRetina);
    % save('spikeResp_primaLandolt.mat','spikeResp');
    
    % load(filterFile);
    spikeAug(1,:) = ones(1,size(spikeResp,2));
    spikeAug(2:9717,:) = spikeResp;
%     movRecon = filterMat'*spikeAug;
%     toc
%     movTrials(:,iTrial) = -movRecon(:,6);
    
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     figure; ieMovie(movReconPlay);
    % save('hallwayReconMovie.mat','movRecon');
    
    if isstruct(filterMat)
        droputIndices = filterMat.dropoutIndices;
        filterMatSm = filterMat.filterMat;
        filterMat=zeros(9717,10000);
        filterMat(droputIndices,:)=filterMatSm;
        filterMatSm=[];        
        
%         spikeAugFull = spikeAug;
%         spikeAug = zeros(size(spikeAug));
%         spikeAug(droputIndices,:) = spikeAugFull(droputIndices,:);
    end
    
    movRecon = filterMat'*spikeAug;
    toc
    movTrials(:,iTrial) = movRecon(:,8);
    spikeTrials(:,iTrial) = spikeAug(:,8);
    clear spikeResp spikeAug


    trialReconPlay = reshape(movTrials,[100 100 size(movTrials,2)]);
    
    if mod(iTrial,50)==0
        % save(fullfile(reconstructionRootPath,'dat','july17gratings','primaH',[sprintf('primaRecon_%4d.mat',10000*contrast)]),'trialReconPlay');
        
        if gap < 0
            save(fullfile(reconstructionRootPath,'dat',folderNameTest,'cont',[sprintf('primaRecon_freqH%4d_cont%6d_rs%2d.mat',100*gratingSpFreq,1e5*contrast,100*resizeFraction)]),'trialReconPlay','contrast','spikeTrials');
            
        else
            save(fullfile(reconstructionRootPath,'dat',folderNameTest,'cont',[sprintf('primaRecon_freq%4d_cont%6d_rs%2d.mat',100*gratingSpFreq,1e5*contrast,100*resizeFraction)]),'trialReconPlay','contrast','spikeTrials');
            
        end
        
    end

end
%% Find best frame
% gratingsMovieImage = landoltMovieContrast(:,:,2);
% movRecon0 = movRecon - ones(10000,1)*mean(movRecon);
% stdRecon0 = std(movRecon0);
% std0 = std(single(gratingsMovieImage(:)));
% movReconNorm = std0*movRecon0./(ones(10000,1)*stdRecon0);
% % movRecon0 = movRecon - ones(10000,1)*me an(movRecon);
% 
% errMov = (movReconNorm - (single(gratingsMovieImage(:))*ones(1,16)));
% 
% % errMov = (movReconNorm - RGB2XWFormat(gratingsMovieContrast));
% size(errMov)
% mse1 = mean(errMov.^2);
% figure; plot(mse1)
% [mv,mi] = min(mse1)