% t_reconDemo
% 
% The recon object is used to reconstruct a movie stimulus from RGC spikes.
% This demo computes the RGC response to a set of image stimuli with
% isetbio, and then generates the reconstructed images using pre-trained
% reconstruction filters. The mean RMS error and median reconstruction
% correlation coefficient (over every image) are reported and compared with
% a reference value.
% 
% TOOLBOX DEPENDENCIES - these must be downloaded and added to the matlab
%                               path with subfolders.
%       isetbio:            http://github.com/isetbio/isetbio [branch:phosphenePaper]
%       RGC-Reconstruction: https://github.com/Chichilnisky-Lab/RGC-Reconstruction
%       EJLPhosphene:       https://github.com/isetbio/EJLPhosphene
%       RemoteDataToolbox:  https://github.com/isetbio/RemoteDataToolbox
% 
% (c) isetbio team 2017 JRG 

%% 
clear;

% Choose demo  type
demoType = 'healthy';
% demoType = 'prosLearning';
% demoType = 'prosNoLearning'

% Set compute flag
% 0 means RGC spikes have already been computed once with this script,
% which takes 10-12 minutes on a laptop.
% 1 bypasses spike computation, script takes only 1 minute
spikeComputeFlag = 0;

%% Build objects

% Choose a folder name, will be created in reconRootPath/dat/ directory
folderName = ['demo_' demoType]; 

% Stimulus is 25 natural scene images presented for 20 frames each at 0.001 sec
pRecon.stimTypeBuild = 'ns500'; 

% Set destination files for movie and spike data
pRecon.stimFile  = fullfile(folderName, 'mov');
pRecon.respFile  = fullfile(folderName, 'sp');
pRecon.buildFile = fullfile(folderName, 'raw','build');

%% Compute spikes

% Create recon object
reconHealthy = recon(pRecon);

if spikeComputeFlag
switch demoType
    
    case 'healthy'
        
        % Use local filter file - change folder name here
        % filterFolderName = 'aug27';
        % mosaicFile = 'mosaic0';        
        % windowSize = 1;
        % percentSV = .5;
        % shifttime = 15;
        % dropout = 0;
        % filterFile  = fullfile(filterFolderName,...    
        %     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV)...
        %         sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);
        % pRecon.filterFile = filterFile;

        % Build RGC mosaic, find spikes in repsonse to stimulus
        reconHealthy.build(pRecon);
        
    case 'prosLearning'
        % Build RGC mosaic with subretinal prosthesis, find spikes in repsonse to stimulus
        pRecon.pixelWidth = 70; % microns 
        reconHealthy.buildPrima(pRecon);
        
    case 'prosNoLearning'
        % Build RGC mosaic with subretinal prosthesis, find spikes in repsonse to stimulus
        pRecon.onlyOnFlag = 1;
        reconHealthy.buildPrima(pRecon);
end

% Format spikes
reconHealthy.loadSpikes(pRecon);
end

% reconHealthy.plot('filters');

%% Compute test measures
% Compute normalized RMS error and corr coeff for each image.
% 
% Note that these measurements differ slightly from values in the paper
% because they are over only 25 images. The test set responses take a long
% time to compute, so only a sample is presented here.

[mse1, cc] = reconHealthy.testImagenet(pRecon);

switch demoType
    case 'healthy'        
        disp('Healthy: 0.1025    0.8305'); % mean       
%         disp('Healthy: 0.0952    0.8305'); % median       
    case 'prosLearning'
        disp('Prosthesis, learning: 0.1551    0.6629'); % mean
%         disp('Prosthesis, learning: 0.1420    0.6629'); % median
    case 'prosNoLearning'        
        disp('Prosthesis, no learning: 0.1530    0.6796'); % mean
%         disp('Prosthesis, no learning: 0.1462    0.6796'); % median
end
