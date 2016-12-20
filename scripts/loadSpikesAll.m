function [movieFile, spikesFile] = loadSpikesAll(varargin)
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
p.addParameter('loadFile',[],@ischar);
p.addParameter('movieFile',[],@ischar);
p.addParameter('spikesFile',[],@ischar);
p.addParameter('mosaicFile',[],@ischar);
p.parse(varargin{:});
loadFile = p.Results.loadFile;
movieFile = p.Results.movieFile;
spikesFile = p.Results.spikesFile;
mosaicFile = p.Results.mosaicFile;

if isempty(movieFile)
    movieFile = ['movie_' num2str(round(cputime*100))];
end

if isempty(spikesFile)
    spikesFile = ['spikes_' num2str(round(cputime*100))];
end

%% 
if ismac || isunix
    dNames = (dir([phospheneRootPath '/dat/' loadFile '*block_*' mosaicFile '.mat']));
else
    dNames = (dir([phospheneRootPath '\dat\' loadFile '*block_*' mosaicFile '.mat']));
end
blocklength = 12000;
numReps = length(dNames);
% numCells= 36+64+169+225;
% spikeResp = zeros(numCells, blocklength*numReps);
if isunix || ismac
    filename1 = [phospheneRootPath '/dat/' loadFile '_block_' num2str(1) mosaicFile '.mat'];
else
    filename1 = [phospheneRootPath '\dat\' loadFile '_block_' num2str(1) mosaicFile '.mat'];
end

matf = matfile(filename1);
szMov = size(matf.whiteNoiseSmall);
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
        filename1 = [phospheneRootPath '/dat/' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
    else
        filename1 = [phospheneRootPath '\dat\' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
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
    
    save([reconstructionRootPath '/dat/' spikesFile mosaicFile ],'spikeResp','-v7.3');
    save([reconstructionRootPath '/dat/' movieFile mosaicFile],'stim','-v7.3')
else
    
    save([reconstructionRootPath '\dat\' spikesFile mosaicFile],'spikeResp','-v7.3');
    save([reconstructionRootPath '\dat\' movieFile mosaicFile],'stim','-v7.3')
end

% save([reconstructionRootPath '/dat/NSspikeResp_40reps_ns0'],'spikeResp','-v7.3');
% save([reconstructionRootPath '/dat/NSmovie_40reps_ns0'],'stim','-v7.3')
% save([reconstructionRootPath '\dat\WNspikeResp_70reps_overlap0'],'spikeResp','-v7.3');
% save([reconstructionRootPath '\dat\WNmovie_spikeResp_70reps_overlap0'],'stim','-v7.3');

% % save('../dat/spikeResp_offParasol','spikeResp')
% save('../dat/spikeResp_onParasol_fast','spikeResp')
% % save('../dat/spikeResp_onMidget','spikeResp')
% % save('spikeResp_onMidget_long300','spikeResp')
% save('../dat/movie_spikeResp_onParasol_fast','stim','-v7.3')
% 
% %figure; imagesc(spikesout*spikesout');

