function obj = buildGratings(obj, varargin)
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
gratingSpFreq = p.Results.gratingSpFreq;

horizontalFlag = p.Results.horizontalFlag;
if isempty(buildFile)
    buildFile = ['primaLandolt_training_' num2str(round(cputime*100))];
end

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

numCols = 100; numRows = 100; % gratingSpFreq = 50;
% gratingsMovieRow = sin((2*pi/(numCols/gratingSpFreq))*[1:numCols]);
gratingsMovieRow = sin((2*pi/(numCols/gratingSpFreq))*[1:numCols] + pi/2);
% gratingsMovieIm= (ones(numRows,1)*gratingsMovieRow)';%*tempFreq;
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


mosaicFile = '_mosaic0';
windowSize = 1;
shiftArr = 4;
shiftTime = shiftArr; shiftind = 1;
percentSVarr =[.25]; percentSVind = 1;
trainSizeArr = [.8]; trainSizeInd = 1;

percentSV = percentSVarr(percentSVind);
shifttime = shiftArr(shiftind);
trainSize = trainSizeArr(trainSizeInd);
filterFile = ['may22/filters_wmean/filters'  mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_tr%d',100*trainSize)];
% filterFile = ['june16prima/filters/filters'  mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_tr%d',100*trainSize)];

load(filterFile);

%% Load image
clear coneParams
coneParams.fov = fov;

coneParams.cmNoiseFlag = 'random';
coneParams.osNoiseFlag = 'random';
iStimNS = ieStimulusMovieCMosaic((gratingsMovieContrast),coneParams);
cMosaicNS = iStimNS.cMosaic;

% Bipolar
clear bpMosaic

cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
for cellTypeInd = 1:4
    clear bpParams
    bpParams.cellType = cellType{cellTypeInd};
    
    bpParams.ecc = patchEccentricity;
    bpParams.rectifyType = 1;
    bpMosaic{cellTypeInd} = bipolar(cMosaicNS, bpParams);
    bpMosaic{cellTypeInd}.set('sRFcenter',1);
    bpMosaic{cellTypeInd}.set('sRFsurround',0);
    bpMosaic{cellTypeInd}.compute(cMosaicNS);
end

% RGC
clear params rgcParams
params.eyeRadius = patchEccentricity;
params.eyeAngle = 90;
innerRetina=ir(bpMosaic,params);
cellType = {'on parasol','off parasol','on midget','off midget'};

rgcParams.centerNoise = 0;
rgcParams.model = 'LNP';
%     rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
rgcParams.axisVariance = 0;
%     rgcParams.rfDiameter = 2;
rgcParams.type = cellType{1};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{2};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{3};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{4};
innerRetina.mosaicCreate(rgcParams);

for iTrial = 1:nTrials
    tic
    iTrial
    innerRetina.compute(bpMosaic);
%     innerRetina.mosaic{4}.window
    
    %% Get spikes and make reconstruction
    
    spikeResp = mosaicSpikes(innerRetina);
    % save('spikeResp_primaLandolt.mat','spikeResp');
    
    % load(filterFile);
    spikeAug(1,:) = ones(1,size(spikeResp,2));
    spikeAug(2:9807,:) = spikeResp;
    movRecon = filterMat'*spikeAug;
    
    movTrials(:,iTrial) = -movRecon(:,8); % tested 8th frame w code below
    
    movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
    figure; ieMovie(movReconPlay);
    % save('hallwayReconMovie.mat','movRecon');
    
    % lambda = .05;
    % filterMat2 = zeroFilter(filterMat,lambda);
    % movRecon2 = filterMat2'*spikeAug;
    toc
    ph=1;
end

trialReconPlay = reshape(movTrials,[100 100 size(movTrials,2)]);

% save(fullfile(reconstructionRootPath,'dat','july17gratings','healthyH',[sprintf('recon_%4d.mat',10000*contrast)]),'trialReconPlay');


if horizontalFlag
    save(fullfile(reconstructionRootPath,'dat','july17gratings','healthyFreqH',[sprintf('recon_freqH%4d.mat',100*gratingSpFreq)]),'trialReconPlay');

else
        save(fullfile(reconstructionRootPath,'dat','july17gratings','healthyFreq',[sprintf('recon_freq%4d.mat',100*gratingSpFreq)]),'trialReconPlay');

end
%% Find best frame
% gratingsMovieImage = gratingsMovieContrast(:,:,2);
% movRecon0 = movRecon - ones(10000,1)*mean(movRecon);
% stdRecon0 = std(movRecon0);
% std0 = std(single(gratingsMovieContrast(:)));
% movReconNorm = std0*movRecon0./(ones(10000,1)*stdRecon0);
% % movRecon0 = movRecon - ones(10000,1)*me an(movRecon);
% 
% errMov = (movReconNorm - (single(gratingsMovieImage(:))*ones(1,32)));
% 
% % errMov = (movReconNorm - RGB2XWFormat(gratingsMovieContrast));
% size(errMov)
% mse = mean(errMov.^2);
% figure; plot(mse)