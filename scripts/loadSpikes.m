clear all;
folderLocation  = '../dat/may26_off/';
blocklength = 2400;
spikeResp = zeros(64, blocklength*150);
stim = zeros(96*96,blocklength*150,'uint8');
for blockNum = 1:150
    blockNum
    filename1 = [folderLocation 'WNstim_response_OffParasol_block_' num2str(blockNum) '.mat'];
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
    stimtmp = reshape(matf.whiteNoiseSmall,96*96,blocklength);
    stim(:,(blockNum-1)*blocklength +1 : blockNum*blocklength) = stimtmp;

    
    % Stimulus here
    % whiteNoiseSmall;
end
save('../dat/spikeResp_may26_off','spikeResp')
%save('../dat/movie_may26','stim','-v7.3')
%figure; imagesc(spikesout*spikesout');

