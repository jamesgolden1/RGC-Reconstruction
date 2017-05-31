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
p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;

if isempty(buildFile)
    buildFile = ['NS_training_' num2str(round(cputime*100))];
end

tic
%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 3.2/4;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
nSteps = 20;
nBlocks = 15;%30;
rng(1504);

for blockNum =1%:nBlocks
    tic
    
    blockNum
%     load(filenameRGC, 'innerRetina');
    
    %      load([ reconstructionRootPath  '/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
%     load('/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_1.mat','movsm');
    %           load([ reconstructionRootPath  '\dat\imagenetBlocks\movsm_' num2str(blockNum) '.mat'],'movsm');
    
%     natScenes = movsm(1:100,1:100,randperm(nSteps));
    
    rsFactor = 1; stimSize = 100;
%     load('C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\stimuli\hallMovie.mat')
%     load([phospheneRootPath '/dat/stimuli/hallMovie.mat'])
    load('/Users/james/Documents/matlab/EJLPhosphene/dat/stimuli/hallMovie.mat')
    szFrames = size(vidFrame,3);
    hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
    for ii = 1:szFrames
        hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
    end
    
    % Set hallway movie stimulus
    testmovieshort = (255*ieScale(hallMovieResize)); clear hallMovieResize;
    
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
    
    
    iStimNS = ieStimulusMovieCMosaic(testmovieshort(:,:,1:50),coneParams);
    cMosaicNS = iStimNS.cMosaic;
    
    %% Bipolar
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
    
    %% RGC
    clear params rgcParams
    params.eyeRadius = .5;%patchEccentricity;
    params.eyeAngle = 90;
    innerRetina=ir(bpMosaic,params);
    cellType = {'on parasol','off parasol','on midget','off midget'};
    
    rgcParams.centerNoise = 0;
    rgcParams.model = 'LNP';
%     rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
    rgcParams.axisVariance = -.1;
%     rgcParams.rfDiameter = 2;
    rgcParams.type = cellType{1};
    innerRetina.mosaicCreate(rgcParams);
%     rgcParams.type = cellType{2};
%     innerRetina.mosaicCreate(rgcParams);
%     rgcParams.type = cellType{3};
%     innerRetina.mosaicCreate(rgcParams);
%     rgcParams.type = cellType{4};
%     innerRetina.mosaicCreate(rgcParams);
    
    innerRetina.compute(bpMosaic);
%     innerRetina.mosaic{4}.window
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
        filename1 = [reconstructionRootPath '/dat/test/' buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    else
        filename1 = [reconstructionRootPath '\dat\ns100/' buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    end
    save(filename1, 'spikesoutsm','whiteNoiseSmall');
    toc
    close all
end


%%%%%%%%%%%%%%%%%%%%%

%%%%% loadSpikesAll %%%%5
% load('test0_block_1_mosaicAll_23282.mat')
%%
if isempty(stimFile)
    stimFile = ['movie_' num2str(round(cputime*100))];
end

if isempty(respFile)
    respFile = ['spikes_' num2str(round(cputime*100))];
end

loadFile = buildFile;
%% 
if ismac || isunix
%     dNames = (dir([phospheneRootPath '/dat/' loadFile '*block_*' mosaicFile '.mat']));
 dNames = (dir([reconstructionRootPath '/dat/' loadFile '*block_*' mosaicFile '.mat']));
else
%     dNames = (dir([phospheneRootPath '\dat\' loadFile '*block_*' mosaicFile '.mat']));
 dNames = (dir([reconstructionRootPath '\dat\' loadFile '*block_*' mosaicFile '.mat']));
end
% blocklength = 12000;
numReps = length(dNames);
% numCells= 36+64+169+225;
% spikeResp = zeros(numCells, blocklength*numReps);
if isunix || ismac
%     filename1 = [phospheneRootPath '/dat/' loadFile '_block_' num2str(1) mosaicFile '.mat'];
    filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(1) mosaicFile '.mat'];
else
%     filename1 = [phospheneRootPath '\dat\' loadFile '_block_' num2str(1) mosaicFile '.mat'];
filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(1) mosaicFile '.mat'];
end

matf = matfile(filename1);
szMov = size(matf.whiteNoiseSmall);
blocklength = szMov(3);
stim = zeros(szMov(1)*szMov(2),blocklength*numReps,'uint8');
clear matf
blockNum = 0;
for blockNumInd =[1:length(dNames) ]
% for blockNumInd =[1:12 21:50]
    blockNum = blockNum+1
    % filename1 = [reconstructionRootPath '\dat\WNstim_response_stx2_block_' num2str(blockNum) '.mat'];    
    % filename1 = [reconstructionRootPath '\dat\WNstim_response_block_' num2str(blockNumInd) '.mat'];    
    % filename1 = [reconstructionRootPath '/dat/nsResponses/NSstim_response_betha_ns0_block_' num2str(blockNum) '.mat'];    
    % filename1 = [reconstructionRootPath '\dat\NSstim_response_overlap0_block_' num2str(blockNum) '.mat'];

    if isunix || ismac
%         filename1 = [phospheneRootPath '/dat/' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
 filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
    else
%         filename1 = [phospheneRootPath '\dat\' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
 filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
    end
    matf = matfile(filename1);
    spikesoutsm = matf.spikesoutsm;
    % Spikes in this variable for each block
    spikesout = double(spikesoutsm);
    pointer = (blockNum-1)*blocklength;
    
    for i = 1:blocklength
        blocksize = 10;
        endval = i*blocksize;
        if endval > size(spikesout,2)
            endval = size(spikesout,2);
        end
        startval = (i-1)*blocksize + 1;
        spikeResp(:,pointer+i) = sum(spikesout(:,startval:endval),2);
    end
    clear spikesout
    clear spikesoutsm
    szMov = size(matf.whiteNoiseSmall);
    stimtmp = reshape(matf.whiteNoiseSmall,szMov(1)*szMov(2),blocklength);
%     stimtmp = uint8(128+double(stimtmp) - ones(size(stimtmp,1),1)*mean(stimtmp,1));
    stim(:,(blockNum-1)*blocklength +1 : blockNum*blocklength) = stimtmp;

    
    % Stimulus here
    % whiteNoiseSmall;
end


if isunix || ismac
    
    save([reconstructionRootPath '/dat/' respFile mosaicFile ],'spikeResp','-v7.3');
    save([reconstructionRootPath '/dat/' stimFile mosaicFile],'stim','-v7.3')
else
    
    save([reconstructionRootPath '\dat\' respFile mosaicFile],'spikeResp','-v7.3');
    save([reconstructionRootPath '\dat\' stimFile mosaicFile],'stim','-v7.3')
end