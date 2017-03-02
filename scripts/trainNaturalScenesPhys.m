function [mosaicFile, saveFile] = trainNaturalScenesPhys(varargin)
%
% Run binary white noise through the RGC array for the big four RGC cell types.
%       based on trainNaturalScenes100_r2
% inputs:
%   mosaicFile - a string that is used to save the mosaic file
%   saveFile - a string that is used to store the spikes and the movie stim
% 
% See also: trainAndTest.m
%
p = inputParser;
p.addParameter('mosaicFile',[],@ischar);
p.addParameter('saveFile',[],@ischar);
p.parse(varargin{:});
% mosaicFile = p.Results.mosaicFile;
saveFile = p.Results.saveFile;

if isempty(saveFile)
    saveFile = ['NS_training_' num2str(round(cputime*100))];
end

%% Initialize
% clear;
% ieInit;

addpath(genpath('/Volumes/Lab/Users/james/isetbio'));
addpath(genpath('/Volumes/Lab/Users/james/RGC-Reconstruction'));

% % Set path
cd(reconstructionRootPath);

tic
%% Parameters to alter

% Retinal patch eccentricity
% patchEccentricity = 12; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 2*1.6;

% Stimulus length = nSteps*nBlocks;
nSteps = 12000;
nBlocks = 15;%30;

%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = nSteps;
params.row = 100;
params.col = 100;
params.fov = fov;
% % params.vfov = 0.7;

%%% Grating subunit stimulus

iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;

%% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

experimentI   = 1;       % Choose dataset to load parameters and spikes
cellTypeI     = 3;%:2    % Choose On Parasol (1) or Off Parasol (2)
stimulusTestI = 2;%:2     % Choose WN test stimulus (1) or NSEM test stimulus (2)
    
% Switch on the conditions indices
% Experimental dataset
switch experimentI
    case 1; experimentID = '2013-08-19-6';
    otherwise; error('Data not yet available');
end
% The other experimental data will be added to the RDT in the future.

% Stimulus: white noise or natural scene movie with eye movements
switch stimulusTestI
    case 1; stimulusTest = 'WN';
    case 2; stimulusTest = 'NSEM';
end

% Cell type: ON or OFF Parasol
% switch cellTypeI
%     case 1; cellType = 'prosthesis selective';
%     case 2; cellType = 'prosthesis off parasol';
%     case 3; cellType = 'prosthesis on parasol';
% end


%% Show raw stimulus for osIdentity
% % figure;
% % for frame1 = 1:size(whiteNoise.sceneRGB,3)
% %     imagesc(squeeze(whiteNoise.sceneRGB(:,:,frame1,:)));
% %     colormap gray; drawnow;
% % end
% % % close;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('displayrgb');

sceneSize = sceneGet(whiteNoise.scene,'size');
retinalPatchWidth = sensorGet(whiteNoise.sensor,'width','m');
% retinalPatchHeight = sensorGet(whiteNoise.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

% timeStep = sensorGet(whiteNoise.sensor,'time interval','sec');
timeStep = (1/125)/1;
os = osSet(os, 'timeStep', timeStep);

os = osSet(os, 'rgbData', whiteNoise.sceneRGB);
% os = osCompute(absorptions);

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

retinalPatchSize = osGet(os,'size');

%% Build RGC array

% if isempty(mosaicFile)
    % clear paramsIR innerRetina
    
%     paramsIR.name    = 'Macaque inner retina 1'; % This instance
%     paramsIR.eyeSide   = 'left';   % Which eye
%     paramsIR.eyeRadius = 1.8;        % Radius in mm
%     paramsIR.eyeAngle  = 90;       % Polar angle in degrees
%     
%     model   = 'LNP';    % Computational model
%     innerRetina = irCreate(os,paramsIR);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
    
    
    % Generate RGC object for simulated GLM prediction of response
    % Set the parameters for the inner retina RGC mosaic. For the inner retina
    % type irPhys, the values for eyeSide, eyeRadius and eyeAngle have no
    % effect, because those are dependent on the properties of the retinal
    % piece used in the Chichilnisky Lab experiment.
    
    % Set parameters
    paramsIR.name = 'macaque phys';
    paramsIR.eyeSide = 'left';
    paramsIR.eyeRadius = 12;
    paramsIR.eyeAngle = 0; ntrials = 0;
    
    % Determined at beginning to allow looping
    paramsIR.experimentID = experimentID; % Experimental dataset
    paramsIR.stimulusTest = stimulusTest; % WN or NSEM
    
    cellType = 'prosthesis on parasol';
    paramsIR.cellType = cellType;         % ON or OFF Parasol    
    innerRetina = irPhys(os, paramsIR);
    nTrials = 1;
    innerRetina.mosaic{1} = innerRetina.mosaic{1}.set('numberTrials',nTrials);
    for cellnum = 1:length(innerRetina.mosaic{1}.cellLocation); newtCenter{cellnum} = -innerRetina.mosaic{1}.tCenter{cellnum}; end;
    innerRetina.mosaic{1} = mosaicSet(innerRetina.mosaic{1},'tCenter',newtCenter);
    
    cellType = 'prosthesis off parasol';
    paramsIR.cellType = cellType;         % ON or OFF Parasol    
    innerRetina2 = irPhys(os, paramsIR);%     
    innerRetina2.mosaic{1} = innerRetina2.mosaic{1}.set('numberTrials',nTrials);
%     innerRetina.mosaic{2} = innerRetina2.mosaic{1}; clear innerRetina2;
    
    
%     cellType = 'prosthesis on midget';
%     paramsIR.cellType = cellType;         % ON or OFF Parasol    
%     innerRetina3 = irPhys(os, paramsIR);
%     innerRetina3.mosaic{1} = innerRetina3.mosaic{1}.set('numberTrials',nTrials);
%     
%     innerRetina.mosaic{3} = innerRetina3.mosaic{1}; clear innerRetina3;
    
    % innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);
    
    % innerRetina.mosaic{1}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{2}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{3}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{4}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{5}.mosaicSet('numberTrials',1);
    
    % irPlot(innerRetina,'mosaic');
    
    mosaicFile = ['mosaicAll_' num2str(round(cputime*100))];
    filenameRGC = [reconstructionRootPath '/dat/nsPhys/' mosaicFile '.mat'];
%     innerRetina.mosaic{1}
%     innerRetina.mosaic{2}
%     innerRetina.mosaic{3}
%     innerRetina.mosaic{4}
    save(filenameRGC, 'innerRetina','innerRetina2');
%     mosaicFile = 'mosaicAll_1246640';
%     filenameRGC = [reconstructionRootPath '/dat/ns100_r2_10_regmos/' mosaicFile '.mat'];
% else
%     
%     % filenameRGC = [reconstructionRootPath '\dat\mosaic_all.mat'];
%     % filenameRGC = [reconstructionRootPath '/dat/mosaic_all_bertha_ns0.mat'];
%     filenameRGC = [reconstructionRootPath '/dat/' mosaicFile '.mat'];
% end

for blockNum =1:nBlocks
    tic
    % clear psthNorm spikesout spikesoutM spikesoutsm whiteNoiseSmall whiteNoise iStim absorptions innerRetina
    
    
    blockNum
    %%% Grating subunit stimulus
    load(filenameRGC, 'innerRetina');
%     iStim = ieStimulusBinaryWhiteNoise(params);
%     absorptions = iStim.sensor;
%     whiteNoise = iStim;
    
    % % If training mosaics separately and need to load white noise movie
    % filename1 = [reconstructionRootPath '\dat\WNstim_response_OnParasol_block_' num2str(blockNum) '.mat'];
    % load(filename1); clear spikesoutsm;
    % whiteNoise.sceneRGB = double(whiteNoiseSmall);
    
    
%      load([ reconstructionRootPath  '\dat\movsm_' num2str(blockNum) '.mat'],'movsm');

     load([ reconstructionRootPath  '/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
%      load([ reconstructionRootPath  '/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
%           load([ reconstructionRootPath  '\dat\imagenetBlocks\movsm_' num2str(blockNum) '.mat'],'movsm');
    natScenes = movsm(1:40,1:80,randperm(nSteps));
    natScenes(:,:,nSteps+1:2*nSteps) = movsm(50+[1:40],1:80,randperm(nSteps));
    os = osSet(os, 'rgbData', double(natScenes));
    
    innerRetina = irCompute(innerRetina,os);

%     spikesout  = RGB2XWFormat(mosaicGet(innerRetina.mosaic{1},'spikes'));
%     clear innerRetina
    
    innerRetina2 = irCompute(innerRetina2,os);
    
%     spikesout2 = RGB2XWFormat(mosaicGet(innerRetina2.mosaic{1},'spikes'));
%     clear innerRetina2
%     innerRetina3 = irCompute(innerRetina3,os);
    
    % irPlot(innerRetina, 'linear');
    % irPlot(innerRetina, 'psth');
    toc
    %% Look at covariance matrix
    tic
    spikesout  = RGB2XWFormat(mosaicGet(innerRetina.mosaic{1},'spikes'));
    spikesout2 = RGB2XWFormat(mosaicGet(innerRetina2.mosaic{1},'spikes'));
%     spikesout3 = RGB2XWFormat(mosaicGet(innerRetina3.mosaic{3},'spikes'));
%     spikesout4 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{4},'spikes'));
    
    timeBins = max([size(spikesout,2) size(spikesout2,2)]);% size(spikesout3,2) ]);%size(spikesout4,2)]);
    
%     spikesoutsm = zeros(size(spikesout,1)+ size(spikesout2,1)+size(spikesout3,1)+size(spikesout4,1), timeBins,'uint8');
    
    spikesoutsm = zeros(size(spikesout,1)+ size(spikesout2,1), timeBins,'uint8');
    spikesoutsm(1:size(spikesout,1) ,1:size(spikesout,2) ) = spikesout;
    spikesoutsm(size(spikesout,1)+[1:size(spikesout2,1)],1:size(spikesout2,2) ) = spikesout2;
    
%     spikesoutsm(size(spikesout,1)+size(spikesout2,1)+[1:size(spikesout3,1)] ,1:size(spikesout3,2) ) = spikesout3;
%     spikesoutsm(size(spikesout,1)+size(spikesout2,1)+size(spikesout3,1)+[1:size(spikesout4,1)] ,1:size(spikesout4,2) ) = spikesout4;
    
%     whiteNoiseSmall = uint8(squeeze(whiteNoise.sceneRGB(:,:,:,1)));
    whiteNoiseSmall = natScenes;
    
%     spikesoutsm = uint8(spikesoutmat);
    
    % filename1 = [reconstructionRootPath '\dat\NSstim_response_overlap0_block_' num2str(blockNum) '.mat'];
    % filename1 = [reconstructionRootPath '/dat/nsResponses/NSstim_response_betha_ns0_block_' num2str(blockNum) '.mat'];
%     filename1 = [reconstructionRootPath '/dat/nsResponses/' saveFile '_block_' num2str(blockNum) '.mat'];
    if ismac || isunix
        filename1 = [reconstructionRootPath '/dat/nsPhys/' saveFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    else
        filename1 = [reconstructionRootPath '\dat\ns100/' saveFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    end
    save(filename1, 'spikesoutsm','whiteNoiseSmall');
    toc
    close all
    figure; imagesc(double(spikesoutsm)*double(spikesoutsm')); colormap parula; drawnow;
end