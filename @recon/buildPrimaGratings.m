function obj = buildPrimaGratings(obj, varargin)
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

p.addParameter('gratingSpFreq',25,@isnumeric);
p.addParameter('horizontalFlag',false,@islogical);
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

numCols = 100; numRows = 100; % gratingSpFreq = 50;
% gratingsMovieRow = sin((2*pi/(numCols/gratingSpFreq))*[1:numCols]);
gratingsMovieRow = sin((2*pi/(numCols/gratingSpFreq))*[1:numCols] + pi/2);
if horizontalFlag
    gratingsMovieIm= (ones(numRows,1)*gratingsMovieRow)';%*tempFreq;
else
    
    gratingsMovieIm= (ones(numRows,1)*gratingsMovieRow);%*tempFreq;
end
    
gratingsMovie = repmat(gratingsMovieIm,[1 1 nSteps]);

% Need factor of two on contrast for accurate match
gratingsMovieContrast = (2*contrast)*(.5*gratingsMovie(1:100,1:100,:));
gratingsMovieContrast(end,end,:) = .5; gratingsMovieContrast(1,end,:) = -.5;
gratingsMovieContrast(:,:,1) = zeros(100,100);
gratingsMovieContrast = 128+127*gratingsMovieContrast;
%% Reconstruction filter parameters

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

primaParams.pixelWidth = 2*35e-6; % meters
primaParams.ecc = 1.8;       % mm
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance


primaRecon = primaArray(gratingsMovieContrast,primaParams);

primaRecon.compute(gratingsMovieContrast)
    
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
    movTrials(:,iTrial) = -movRecon(:,6);
    
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     figure; ieMovie(movReconPlay);
    % save('hallwayReconMovie.mat','movRecon');
    
    % lambda = .05;
    % filterMat2 = zeroFilter(filterMat,lambda);
    % movRecon2 = filterMat2'*spikeAug;
    clear spikeResp spikeAug
end

trialReconPlay = reshape(movTrials,[100 100 size(movTrials,2)]);


% save(fullfile(reconstructionRootPath,'dat','july17gratings','primaH',[sprintf('primaRecon_%4d.mat',10000*contrast)]),'trialReconPlay');

if horizontalFlag
    save(fullfile(reconstructionRootPath,'dat','july17gratings','primaFreqH',[sprintf('primaRecon_freqH%4d.mat',100*gratingSpFreq)]),'trialReconPlay');

else
        save(fullfile(reconstructionRootPath,'dat','july17gratings','primaFreq',[sprintf('primaRecon_freq%4d.mat',100*gratingSpFreq)]),'trialReconPlay');

end
%% Find best frame
% gratingsMovieImage = gratingsMovieContrast(:,:,2);
% movRecon0 = movRecon - ones(10000,1)*mean(movRecon);
% stdRecon0 = std(movRecon0);
% std0 = std(single(gratingsMovieContrast(:)));
% movReconNorm = std0*movRecon0./(ones(10000,1)*stdRecon0);
% % movRecon0 = movRecon - ones(10000,1)*me an(movRecon);
% 
% errMov = (-movReconNorm - (single(gratingsMovieImage(:))*ones(1,16)));
% 
% % errMov = (movReconNorm - RGB2XWFormat(gratingsMovieContrast));
% size(errMov)
% mse = mean(errMov.^2);
% figure; plot(mse)