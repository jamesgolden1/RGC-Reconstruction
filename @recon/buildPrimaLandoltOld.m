function obj = buildPrimaLandoltOld(obj, varargin)
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

p.addParameter('gap',0,@isnumeric);
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
gap = p.Results.gap;
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

orient0 = 'down';
gapsize0 = gap;
landoltStim = imgLandoltC('orientation',orient0,'gapsize',gapsize0);
landoltMovie = repmat(landoltStim,[1 1 nSteps]);

% contarr = [.05 .01 .025 .05 .005]; contind = 1;
landoltMovieContrast = (contrast)*(landoltMovie(1:100,1:100,:) - 0.5);%+0.5;
landoltMovieContrast(end,end,:) = .5; landoltMovieContrast(1,end,:) = -.5;
landoltMovieContrast(:,:,1) = zeros(100,100);
            
%% Reconstruction filter parameters


mosaicFile = '_mosaic0';
windowSize = 1;
shiftArr = 2;
shiftTime = shiftArr; shiftind = 1;
percentSVarr =[.25]; percentSVind = 1;
trainSizeArr = [.8]; trainSizeInd = 1;

percentSV = percentSVarr(percentSVind);
shifttime = shiftArr(shiftind);
trainSize = trainSizeArr(trainSizeInd);
% filterFile = ['may22/filters_wmean/filters'  mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_tr%d',100*trainSize)];
filterFile = ['june16prima/filters/filters'  mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_tr%d',100*trainSize)];

load(filterFile);

%% Load image

primaParams.pixelWidth = 1*35e-6; % meters
primaParams.ecc = 1.8;       % mm
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance


primaRecon = primaArray(landoltMovieContrast,primaParams);

primaRecon.compute(landoltMovieContrast)
    
for iTrial = 1:nTrials
    iTrial
    tic
    % primaRecon = primaArray(landoltMovieContrast,primaParams);
    
    primaRecon.computeRGC()
    innerRetina = primaRecon.innerRetina;
    
    %% Get spikes and make reconstruction
    
    spikeResp = mosaicSpikes(primaRecon.innerRetina);
    % save('spikeResp_primaLandolt.mat','spikeResp');
    
    % load(filterFile);
    spikeAug(1,:) = ones(1,size(spikeResp,2));
    spikeAug(2:9807,:) = spikeResp;
    movRecon = filterMat'*spikeAug;
    toc
    movTrials(:,iTrial) = movRecon(:,8);
    
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     figure; ieMovie(movReconPlay);
    % save('hallwayReconMovie.mat','movRecon');
    
    % lambda = .05;
    % filterMat2 = zeroFilter(filterMat,lambda);
    % movRecon2 = filterMat2'*spikeAug;
    clear spikeResp spikeAug
end

trialReconPlay = reshape(movTrials,[100 100 size(movTrials,2)]);

save(fullfile(reconstructionRootPath,'dat','landoltC',[sprintf('primaLandoltRecon_gap%1d_%2d.mat',2+(gapsize0),100*contrast)]),'trialReconPlay');