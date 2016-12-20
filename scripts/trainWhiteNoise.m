function [mosaicFile, saveFile] = trainWhiteNoise(varargin)
%
% Run binary white noise through the RGC array for the big four RGC cell types.
% 
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
mosaicFile = p.Results.mosaicFile;
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
patchEccentricity = 12; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6;

% Stimulus length = nSteps*nBlocks;
nSteps = 12000;
nBlocks = 40;

%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = nSteps;
params.row = 96;
params.col = 96;
params.fov = fov;
% % params.vfov = 0.7;

%%% Grating subunit stimulus

iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;
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
    paramsIR.name    = 'Macaque inner retina 1'; % This instance
    paramsIR.eyeSide   = 'left';   % Which eye
    paramsIR.eyeRadius = 4;        % Radius in mm
    paramsIR.eyeAngle  = 90;       % Polar angle in degrees
    
    model   = 'LNP';    % Computational model
    innerRetina = irCreate(os,paramsIR);
    innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
    innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
    innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
    innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
    
    % innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);
    
    % innerRetina.mosaic{1}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{2}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{3}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{4}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{5}.mosaicSet('numberTrials',1);
    
    % irPlot(innerRetina,'mosaic');
    
    mosaicFile = ['mosaicAll_' num2str(round(cputime*100))];
    filenameRGC = [reconstructionRootPath '/dat/' mosaicFile '.mat'];
    save(filenameRGC, 'innerRetina');

% else
%     
%     % filenameRGC = [reconstructionRootPath '\dat\mosaic_all.mat'];
%     % filenameRGC = [reconstructionRootPath '/dat/mosaic_all_bertha_ns0.mat'];
%     filenameRGC = [reconstructionRootPath '/dat/' mosaicFile '.mat'];
% end

for blockNum = 0+[1:nBlocks]
    
    % clear psthNorm spikesout spikesoutM spikesoutsm whiteNoiseSmall whiteNoise iStim absorptions innerRetina
    
    
    blockNum
    %%% Grating subunit stimulus
    load(filenameRGC, 'innerRetina');
    iStim = ieStimulusBinaryWhiteNoise(params);
    absorptions = iStim.sensor;
    whiteNoise = iStim;
    
    % % If training mosaics separately and need to load white noise movie
    % filename1 = [reconstructionRootPath '\dat\WNstim_response_OnParasol_block_' num2str(blockNum) '.mat'];
    % load(filename1); clear spikesoutsm;
    % whiteNoise.sceneRGB = double(whiteNoiseSmall);
    
    
    
    
    
    
    os = osSet(os, 'rgbData', whiteNoise.sceneRGB);
    
    innerRetina = irCompute(innerRetina,os);
    
    % irPlot(innerRetina, 'linear');
    % irPlot(innerRetina, 'psth');
    
    %% Look at covariance matrix
    
    spikesout  = RGB2XWFormat(mosaicGet(innerRetina.mosaic{1},'spikes'));
    spikesout2 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{2},'spikes'));
    spikesout3 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{3},'spikes'));
    spikesout4 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{4},'spikes'));
    
    timeBins = max([size(spikesout,2) size(spikesout2,2) size(spikesout3,2) size(spikesout4,2)]);
    
    spikesoutmat = zeros(size(spikesout,1)+ size(spikesout2,1)+size(spikesout3,1)+size(spikesout4,1), timeBins);
    spikesoutmat(1:size(spikesout,1) ,1:size(spikesout,2) ) = spikesout;
    spikesoutmat(size(spikesout,1)+[1:size(spikesout2,1)],1:size(spikesout2,2) ) = spikesout2;
    
    spikesoutmat(size(spikesout,1)+size(spikesout2,1)+[1:size(spikesout3,1)] ,1:size(spikesout3,2) ) = spikesout3;
    spikesoutmat(size(spikesout,1)+size(spikesout2,1)+size(spikesout3,1)+[1:size(spikesout4,1)] ,1:size(spikesout4,2) ) = spikesout4;
    
    whiteNoiseSmall = uint8(squeeze(whiteNoise.sceneRGB(:,:,:,1)));
    
    spikesoutsm = uint8(spikesoutmat);
    
%     filename1 = [reconstructionRootPath '\dat\WNstim_response_stx2_block_' num2str(blockNum) '.mat'];
    filename1 = [reconstructionRootPath '\dat\nsResponses/' saveFile '_block_' num2str(blockNum) '.mat'];
    
    save(filename1, 'spikesoutsm','whiteNoiseSmall');
    toc
    close
end