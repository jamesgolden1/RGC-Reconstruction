function obj = buildPrima(obj, varargin)
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
fov = 1.7;% 3.2;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
nSteps = 500;%000;
nBlocks = 15;%30;


for blockNum =blockIn%1%:nBlocks
    tic
    
    blockNum
    %     load(filenameRGC, 'innerRetina');
    
    %      load([ reconstructionRootPath  '/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
    %     load(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
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
        natScenes = 192*round(natScenesRaw); clear natScenesRaw;
%         natScenes = round(192*natScenesRaw); clear natScenesRaw;
            
    end
%               load([ reconstructionRootPath  '\dat\imagenetBlocks\movsm_' num2str(blockNum) '.mat'],'movsm');
%     
%     rsFactor =1; stimSize = 100;
%         load([phospheneRootPath '/dat/stimuli/hallMovie.mat'])
%     szFrames = size(vidFrame,3);
%     hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
%     for ii = 1:szFrames
%         hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
%     end
%     
%     % Set hallway movie stimulus
%     testmovieshort = (255*ieScale(hallMovieResize)); clear hallMovieResize;
    
    %%
    %% Load image       
    
    primaParams.pixelWidth = 1*35e-6; % meters
    primaParams.ecc = 1.8;       % mm
    primaParams.fov = 1.7/1;     % deg
    
    primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
    primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
    primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance
    
    primaRecon = primaArray(natScenes,primaParams);
    
    primaRecon.compute(natScenes)
    innerRetina = primaRecon.innerRetina;
    
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