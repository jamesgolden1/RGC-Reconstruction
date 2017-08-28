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
p.addParameter('stimTypeBuild','ns',@ischar);
p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;
stimTypeBuild = p.Results.stimTypeBuild;

if isempty(buildFile)
    buildFile = ['NS_training_' num2str(round(cputime*100))];
end

tic
%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.7;% 3.2/2;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
nSteps = 500;
nBlocks = 15;%30;
rng(1504);

for blockNum =blockIn%:nBlocks
    tic
    
    blockNum
    %     natScenes = 255*ieScale(loadHallStimulus(20));
    if stimTypeBuild == 'ns'
        if blockNum <= 288
            movsm = parload(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(mod(blockNum-1,12)+1) '.mat']);
            
            natScenes = movsm(1:100,1:100,nSteps*floor((blockNum-1)/12)+randperm(nSteps));
        else
            movsm = parload(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(12+mod(blockNum-1,12)+1) '.mat']);
            
            natScenes = movsm(1:100,1:100,nSteps*(floor((-288+blockNum-1)/12))+randperm(nSteps));
        end
    elseif stimTypeBuild == 'wn'
        natScenesRaw = (rand(100,100,nSteps));
        natScenes = 255*round(natScenesRaw); clear natScenesRaw;
            
    end

%     testInds = [1:10:500-10];
%     natScenesAll = natScenes;
%     natScenes = zeros(size(natScenesAll));
%     for ti = 0:9
%     natScenes(:,:,testInds+ti) = natScenesAll(:,:,testInds);
%     end
    %%
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
    
    
    iStimNS = ieStimulusMovieCMosaic(natScenes,coneParams);
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
    rgcL.compute('bipolarScale',250,'bipolarContrast',1);
    
    %%
    toc
    %% Look at covariance matrix
    tic
    spikesout  = RGB2XWFormat(rgcL.mosaic{1}.get('spikes'));
    spikesout2 = RGB2XWFormat(rgcL.mosaic{2}.get('spikes'));
    spikesout3 = RGB2XWFormat(rgcL.mosaic{3}.get('spikes'));
    spikesout4 = RGB2XWFormat(rgcL.mosaic{4}.get('spikes'));
    
    timeBins = max([size(spikesout,2) size(spikesout2,2) size(spikesout3,2) size(spikesout4,2)]);
    
    spikesoutsm = zeros(size(spikesout,1)+ size(spikesout2,1)+size(spikesout3,1)+size(spikesout4,1), timeBins,'uint8');
    spikesoutsm(1:size(spikesout,1) ,1:size(spikesout,2) ) = spikesout;
    spikesoutsm(size(spikesout,1)+[1:size(spikesout2,1)],1:size(spikesout2,2) ) = spikesout2;
    
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+[1:size(spikesout3,1)] ,1:size(spikesout3,2) ) = spikesout3;
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+size(spikesout3,1)+[1:size(spikesout4,1)] ,1:size(spikesout4,2) ) = spikesout4;
    
    whiteNoiseSmall = natScenes;

    if ismac || isunix
        filename1 = [reconstructionRootPath '/dat/' buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    else
        filename1 = [reconstructionRootPath '\dat\ns100/' buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    end
    save(filename1, 'spikesoutsm','whiteNoiseSmall');
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

% sta = stimzm(:,1:end-14)*single(spikeResp(:,15:end)');
% figure; imagesc(reshape(sta(:,4700),[100 100]))

% % for computing mean spike rate
% s1= (spikeResp(1:837,:)); mean(s1(:))
% sum(length(vertcat(innerRetina.mosaic{2}.cellLocation{:})))
% s1= (spikeResp(838:838+1085,:)); mean(s1(:))
% sum(length(vertcat(innerRetina.mosaic{3}.cellLocation{:})))
%  s1= (spikeResp(838+1085+1:838+1085+3348,:)); mean(s1(:))
% sum(length(vertcat(innerRetina.mosaic{4}.cellLocation{:})))
% s1= (spikeResp(838+1085+3348+1:838+1085+3348+4526,:)); mean(s1(:))