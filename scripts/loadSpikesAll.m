% clear all;
blocklength = 12000;
numReps = 70;
numCells= 36+64+169+225;
% spikeResp = zeros(numCells, blocklength*numReps);
stim = zeros(96*96,blocklength*numReps,'uint8');
blockNum = 0;
for blockNumInd =[1:12 21:50]
    blockNum = blockNum+1
%     filename1 = [reconstructionRootPath '\dat\WNstim_response_stx2_block_' num2str(blockNum) '.mat'];

filename1 = [reconstructionRootPath '\dat\WNstim_response_block_' num2str(blockNumInd) '.mat'];

% filename1 = [reconstructionRootPath '\dat\NSstim_response_overlap0_block_' num2str(blockNum) '.mat'];
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
    stimtmp = reshape(matf.whiteNoiseSmall,96*96,blocklength);
%     stimtmp = uint8(128+double(stimtmp) - ones(size(stimtmp,1),1)*mean(stimtmp,1));
    stim(:,(blockNum-1)*blocklength +1 : blockNum*blocklength) = stimtmp;

    
    % Stimulus here
    % whiteNoiseSmall;
end


save([reconstructionRootPath '\dat\WNspikeResp_70reps_overlap0'],'spikeResp','-v7.3');

save([reconstructionRootPath '\dat\WNmovie_spikeResp_70reps_overlap0'],'stim','-v7.3')

% % save('../dat/spikeResp_offParasol','spikeResp')
% save('../dat/spikeResp_onParasol_fast','spikeResp')
% % save('../dat/spikeResp_onMidget','spikeResp')
% % save('spikeResp_onMidget_long300','spikeResp')
% save('../dat/movie_spikeResp_onParasol_fast','stim','-v7.3')
% 
% %figure; imagesc(spikesout*spikesout');

