% clear all;
% folderLocation  = '../dat/may26_off/';
blocklength = 2400;
numReps = 100
numCells= 36;
% spikeResp = zeros(numCells, blocklength*numReps);
stim = zeros(96*96,blocklength*numReps,'uint8');
for blockNum = 1:numReps
    blockNum
%     filename1 = [folderLocation 'WNstim_response_OffParasol_block_' num2str(blockNum) '.mat'];

    filename1 = ['C:\Users\James\Documents\GitHub\onParasol1_all\WNstim_response_OnParasol_block_' num2str(blockNum) '.mat'];
%     filename1 = ['C:\Users\James\Documents\GitHub\offParasol1_all\WNstim_response_OffParasol_block_' num2str(blockNum) '.mat'];
%     filename1 = ['C:\Users\James\Documents\GitHub\onMidget1\WNstim_response_OnMidget_block_' num2str(blockNum) '.mat'];
%     filename1 = ['C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_block_' num2str(blockNum) '.mat'];
    matf = matfile(filename1);
    spikesoutsm = matf.spikesoutsm;
    % Spikes in this variable for each block
    spikesout = double(spikesoutsm);
    pointer = (blockNum-1)*blocklength;
    
    for i = 1:blocklength
        blocksize = 100;
        endval = i*blocksize;
        if endval > size(spikesout,2)
            endval = size(spikesout,2);
        end
        startval = (i-1)*blocksize + 1;
        spikeResp(:,pointer+i) = sum(spikesout(:,startval:endval),2);
    end
    clear spikesout
    clear spikesoutsm
%     stimtmp = reshape(matf.whiteNoiseSmall,96*96,blocklength);
%     stim(:,(blockNum-1)*blocklength +1 : blockNum*blocklength) = stimtmp;

    
    % Stimulus here
    % whiteNoiseSmall;
end

% save('../dat/spikeResp_offParasol','spikeResp')
save('../dat/spikeResp_onParasol','spikeResp')
% save('../dat/spikeResp_onMidget','spikeResp')
% save('spikeResp_onMidget_long300','spikeResp')
% save('../dat/movie_onMidget_long300','stim','-v7.3')

%figure; imagesc(spikesout*spikesout');

