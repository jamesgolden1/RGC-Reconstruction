function obj = build(obj, varargin)
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
p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;

if isempty(buildFile)
    buildFile = ['NS_training_' num2str(round(cputime*100))];
end

tic
%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 3.2;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
nSteps = 500;%12000;%000;
nBlocks = 15;%30;


for blockNum =blockIn%1%:nBlocks
    tic
    
    blockNum
%     load(filenameRGC, 'innerRetina');
    
    %      load([ reconstructionRootPath  '/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
%     load(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
    movsm = parload(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(mod(blockNum-1,12)+1) '.mat']);
    %           load([ reconstructionRootPath  '\dat\imagenetBlocks\movsm_' num2str(blockNum) '.mat'],'movsm');
    
    natScenes = movsm(1:100,1:100,nSteps*floor((blockNum-1)/12)+randperm(nSteps));
    
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
    iStimNS = ieStimulusMovieCMosaic(natScenes,coneParams);
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
    
    whiteNoiseSmall = natScenes;

    if ismac || isunix
        filename1 = [reconstructionRootPath '/dat/' buildFile '_block_' num2str(blockNum) '' mosaicFile '.mat'];
    else
        filename1 = [reconstructionRootPath '\dat\ns100/' buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    end
    % save(filename1, 'spikesoutsm','whiteNoiseSmall');
    parsave(filename1, spikesoutsm, whiteNoiseSmall);
    toc
    close all
end


% Working STA!
% cd /Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/test4
% load('sp_may5_test4_mosaic_may5.mat')
% load('mov_may5_test4_mosaic_may5.mat')
% stimzm = (single(stim)-(ones(size(stim,2),1)*mean(stim,2)')');
% for j = 2:30; sta2(:,j-1) = ((stimzm(:,1:end-(j-1))))*single(spikeResp(700,j:end)'); end;
% figure; imagesc(sta2);
% figure; ieMovie(reshape(sta2,[100 100 29]));