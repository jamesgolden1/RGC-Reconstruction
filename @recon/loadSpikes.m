function loadSpikes(obj, varargin) %#ok<INUSL>
% loadspike(reconObject, parameters)
% 
% This function takes the saved spikes and downsamples them in order to
% compute the image reconstruction. The spiking code generates spikes in 10
% bins for each time step; if the cone mosaic interval is set to 0.001
% seconds, the spikes are generated in bins of 0.0001 seconds, and are
% downsampled for reconstruction.
% 
% There is some manipulation of the file name in order to save them
% appropriately by blocks when a larger number of files are saved in
% training.
% 
p = inputParser;
p.addParameter('stimFile',[],@ischar);
p.addParameter('respFile',[],@ischar);
p.addParameter('filterFile',[],@ischar);
p.addParameter('windowSize',8,@isnumeric);
p.addParameter('percentSV',0.5,@isnumeric);
p.addParameter('mosaicFile','mosaic0',@ischar);
p.addParameter('trainFraction',1,@isnumeric);
p.addParameter('shiftTime',0,@isnumeric);
p.addParameter('dropout',0,@isnumeric);
p.addParameter('stimType','ns',@ischar);
p.addParameter('buildFile',[],@ischar);

p.addParameter('pixelWidth',70,@isnumeric);
p.addParameter('currentDecay',2,@isnumeric);

p.addParameter('numTest',1,@isnumeric);

p.addParameter('testFlag',0,@isnumeric);

p.KeepUnmatched = true;
p.parse(varargin{:});
filterFile = p.Results.filterFile;
stimFileName = p.Results.stimFile;
respFileName = p.Results.respFile;
windowSize = p.Results.windowSize;
percentSV = p.Results.percentSV;
mosaicFile = p.Results.mosaicFile;
trainSizeArray = p.Results.trainFraction;
shiftTime = p.Results.shiftTime;
stimType = p.Results.stimType;
buildFile = p.Results.buildFile;
dropout = p.Results.dropout;

pixelWidth = p.Results.pixelWidth;
currentDecay = p.Results.currentDecay;
numTest = p.Results.numTest;
testFlag = p.Results.testFlag;

if isempty(filterFile)
    filterFile = ['filters_' num2str(round(cputime*100))];
end
% %

stimFileTest = [ stimFileName '_' mosaicFile];

% stimFileTest = [ stimFileName '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile];

% if exist([stimFileTest '.mat'],'file') ~= 2
%     loadSpikes(buildFile,stimFileName,respFileName,mosaicFile,pixelWidth,currentDecay,numTest)
% end


loadFile = buildFile;

%% Assign a filename if none has been

if isempty(stimFileName)
    stimFile = ['movie_' num2str(round(cputime*100))];
end

if isempty(respFileName)
    respFile = ['spikes_' num2str(round(cputime*100))];
end

loadFile = buildFile;

%% Loop through numbered blocks and starting indices
for fi = 1:numTest
    
    % Find all files in directory
    if ~testFlag
        dNames = dir(fullfile(reconstructionRootPath,'dat', [loadFile '*block_*' mosaicFile '.mat']));
    else
        dNames = dir(fullfile(reconstructionRootPath,'dat', [loadFile '*block_*' '_start_' num2str(fi) '_' mosaicFile '.mat']));
    end
    %  dNames = (dir([reconstructionRootPath '/dat/' loadFile '*block_*' '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat']));
    
    numReps = length(dNames);
    
    % Blocks 1-468 of 500 images are used in training
    % Blocks 469-576 are used in testing
    if isempty(strfind(dNames(1).name,'469'))
        blockStart = 0;
    else
        blockStart = 468;
    end
    
    % Different file names for train vs. test, test includes start index
    if ~testFlag
        filename1 = fullfile(reconstructionRootPath, 'dat', [loadFile '_block_' num2str(blockStart+1) '_' mosaicFile '.mat']);
    else
        filename1 = fullfile(reconstructionRootPath, 'dat', [loadFile '_block_' num2str(blockStart+1) '_start_' num2str(fi) '_' mosaicFile '.mat']);
    end
        
    % We've loaded our files, now downsample spikes and save as new files
    matf = matfile(filename1);
    szMov = size(matf.whiteNoiseSmall);
    blocklength = szMov(3);
    stim = zeros(szMov(1)*szMov(2),blocklength*numReps,'uint8');
    clear matf
    blockNum = 0+blockStart;
    for blockNumInd =[1:length(dNames) ]
        blockNum = blockNum+1        
        
        % Different file names for train vs. test, test includes start index
        if ~testFlag            
            filename1 = fullfile(reconstructionRootPath, 'dat', [loadFile '_block_' num2str(blockNum) '_' mosaicFile '.mat']);
        else
            filename1 = fullfile(reconstructionRootPath, 'dat', [loadFile '_block_' num2str(blockNum) '_start_' num2str(fi) '_' mosaicFile '.mat']);
        end
        
        matf = matfile(filename1);
        spikesoutsm = matf.spikesoutsm;
        % Spikes in this variable for each block
        spikesout = double(spikesoutsm);
        pointer = (blockNum-1-blockStart)*blocklength;
        
        % Downsample 
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
        stim(:,(blockNum-1-blockStart)*blocklength +1 : (blockNum-blockStart)*blocklength) = stimtmp;
        
    end
    
    % Different file names for train vs. test, test includes start index
    if ~testFlag
        save(fullfile(reconstructionRootPath, 'dat', [respFileName '_' mosaicFile '.mat']),'spikeResp','-v7.3');
        save(fullfile(reconstructionRootPath, 'dat', [stimFileName '_' mosaicFile '.mat']),'stim','-v7.3')
    else
        save(fullfile(reconstructionRootPath, 'dat', [respFile '_' num2str(fi) '_' mosaicFile '.mat']),'spikeResp','-v7.3');
        save(fullfile(reconstructionRootPath, 'dat', [stimFile '_' num2str(fi) '_' mosaicFile '.mat']),'stim','-v7.3')
    end
    
    
    
end
