function [rgcL, bpL, cMosaicNS, iStimNS] = mosaicResponse(stimScenes, fov)
% MOSAICRESPONSE


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
% % params.vfov = 0.7;


iStimNS = ieStimulusMovieCMosaic(stimScenes,coneParams);
cMosaicNS = iStimNS.cMosaic;

%% Bipolar
%% Create a set of bipolar cell types in the bipolar mosaic

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


%% RGC

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

nTrials = 1; rgcL.set('numberTrials',nTrials);

%% Compute the inner retina response and visualize

% Every mosaic has its input and properties assigned so we should be able
% to just run through all of them.
rgcL.compute('bipolarScale',50,'bipolarContrast',1);
