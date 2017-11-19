function obj = build(obj, varargin)
%BUILD - builds training set for the recon object
% Run natural scenes through the RGC array for the big four RGC cell types.
% 
% Loads appropriate stimulus based on 'stimTypeBuild' setting with the
% 'generateStimulus' function. 'ns500' runs a set of images for demo
% purposes, where the reconstruction is computed and the error is
% caluculated and compared with benchmarks.'ns' loads natural scenes from
% the imagenet database, but only when the script is being run from the
% appropriate Stanford server.
% 
% 'mosaicResponse' builds the large cone, bipolar and RGC mosaics in
% isetbio used for reconstruction.
% 
% Once the responses are computed, the spikes and stimuli are saved to a
% mat file.
% 
%
p = inputParser;
p.addParameter('mosaicFile','mosaic0',@ischar);
p.addParameter('buildFile',[],@ischar);
p.addParameter('stimFile',[],@ischar);
p.addParameter('respFile',[],@ischar);
p.addParameter('blockIn',1,@isnumeric);
p.addParameter('startInd',1,@isnumeric);
p.addParameter('testFlag',0,@isnumeric);
p.addParameter('stimTypeBuild','ns',@ischar);

p.KeepUnmatched = true;
p.parse(varargin{:});

mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;
startInd = p.Results.startInd;
testFlag = p.Results.testFlag;
stimTypeBuild = p.Results.stimTypeBuild;

if isempty(buildFile)
    buildFile = ['NS_training_' num2str(round(cputime*100))];
end

tic
%% Parameters to alter

% Retinal patch eccentricity
% patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.7;% 3.2/2;

% Stimulus length = nSteps*nBlocks;
nSteps = 500;

rng(1504);

for blockNum =blockIn%:nBlocks
    %% Get stimulus
    
    tic
    
    blockNum
    %     natScenes = 255*ieScale(loadHallStimulus(20));
    
    stimScenes = generateStimulus(stimTypeBuild, blockNum, nSteps, testFlag);
     
    %% Build mosaic and get response
    
    [rgcL, bpL, cMosaicNS, iStimNS] = mosaicResponse(stimScenes, fov);

    toc
    
    %% Get spikes
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
    
    %% Save spikes and movie
    whiteNoiseSmall = stimScenes;

    if ~exist(fullfile(reconstructionRootPath,'dat',buildFile(1:end-5)),'dir')
        mkdir(fullfile(reconstructionRootPath,'dat',buildFile(1:end-5)));
    end

    if testFlag        
            filename1 = fullfile(reconstructionRootPath,'dat',[ buildFile '_block_' num2str(blockNum) '_start_' num2str(startInd) '_' mosaicFile '.mat']);
    else
            filename1 = fullfile(reconstructionRootPath,'dat',[ buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat']);
    end
    
    save(filename1, 'spikesoutsm','whiteNoiseSmall');
    toc
    close all
end
end

%%
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