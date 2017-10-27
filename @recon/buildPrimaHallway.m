function obj = buildPrimaHallway(obj, varargin)
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

p.addParameter('pixelWidth',70,@isnumeric);
p.addParameter('currentDecay',2,@isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;
stimTypeBuild = p.Results.stimTypeBuild;
pixelWidth = p.Results.pixelWidth;
currentDecay = p.Results.currentDecay;

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
nSteps =550;%000;
nBlocks = 15;%30;


    tic
    
    
movieIn = 255*ieScale(loadHallStimulus(nSteps));

    %%
    %% Load image       
    
    primaParams.pixelWidth = pixelWidth*1e-6; % meters
    primaParams.ecc = 1.8;       % mm
    primaParams.fov = 1.7/1;     % deg
    
    primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
    primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
    primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance
    
    primaParams.currentDecay = currentDecay;
    
    primaRecon = primaArray(movieIn,primaParams);
    
    primaRecon.compute(movieIn)
    innerRetina = primaRecon.innerRetina;
    
    %%
    toc
    %% Look at covariance matrix
    tic
%     spikesoutsm = mosaicSpikes(innerRetina);
    
    spikesout   = RGB2XWFormat(innerRetina.mosaic{1}.get('spikes'));    
    spikesout2  = RGB2XWFormat(innerRetina.mosaic{2}.get('spikes'));    
    spikesout3  = RGB2XWFormat(innerRetina.mosaic{3}.get('spikes'));
    spikesout4  = RGB2XWFormat(innerRetina.mosaic{4}.get('spikes'));
    
%     spikesout  = RGB2XWFormat(mosaicGet(innerRetina.mosaic{1},'spikes'));
%     spikesout2 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{2},'spikes'));
%     spikesout3 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{3},'spikes'));
%     spikesout4 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{4},'spikes'));
    
    timeBins = max([size(spikesout,2) size(spikesout2,2) size(spikesout3,2) size(spikesout4,2)]);
    
    spikesoutsm = zeros(size(spikesout,1)+ size(spikesout2,1)+size(spikesout3,1)+size(spikesout4,1), timeBins,'uint8');
    spikesoutsm(1:size(spikesout,1) ,1:size(spikesout,2) ) = spikesout;
    spikesoutsm(size(spikesout,1)+[1:size(spikesout2,1)],1:size(spikesout2,2) ) = spikesout2;
    
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+[1:size(spikesout3,1)] ,1:size(spikesout3,2) ) = spikesout3;
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+size(spikesout3,1)+[1:size(spikesout4,1)] ,1:size(spikesout4,2) ) = spikesout4;
    
    whiteNoiseSmall = movieIn;

    %     filename1 = [reconstructionRootPath '/dat/' buildFile '_block_' num2str(blockNum) '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat'];

    if ismac || isunix
        filename1 = [reconstructionRootPath '/dat/' buildFile '_hallway' '_' mosaicFile '.mat'];
    else
        filename1 = [reconstructionRootPath '\dat\ns100/' buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    end
    %     save(filename1, 'spikesoutsm','whiteNoiseSmall');
    %     toc
    %     close all
    
    %     filename1 = [reconstructionRootPath '/dat/' buildFile '_block_' num2str(blockNum) '_' mosaicFile '.mat']
    save(filename1, 'spikesoutsm','whiteNoiseSmall');
%     save(filename1, spikesoutsm, whiteNoiseSmall);
    
    loadSpikesPitch(buildFile,stimFile,respFile,mosaicFile,pixelWidth,currentDecay,filename1)
    

function loadSpikesPitch(buildFile,stimFile,respFile,mosaicFile,pixelWidth,currentDecay,filename1)
loadFile = buildFile;

% dNames = (dir([reconstructionRootPath '/dat/' loadFile '*block_*' '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) mosaicFile '.mat']));

numReps = 1;%length(dNames);
% filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(451) '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) mosaicFile '.mat'];

matf = matfile(filename1);
szMov = size(matf.whiteNoiseSmall);
blocklength = szMov(3);
stim = zeros(szMov(1)*szMov(2),blocklength*numReps,'uint8');
clear matf
blockNum = 0;

for blockNumInd =1%[1:length(dNames) ]
    % for blockNumInd =[1:12 21:50]
    blockNum = blockNum+1
    
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



save([reconstructionRootPath '/dat/' respFile '_hallway_' mosaicFile ],'spikeResp','-v7.3');
save([reconstructionRootPath '/dat/' stimFile '_hallway_' mosaicFile],'stim','-v7.3')