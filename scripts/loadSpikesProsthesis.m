
% folderLocation  = '../dat/may26_off/';
% folderLocation = '';
blocklength = 120000;
numCells = 48;
numberBlocks = 2;
spikeResp = zeros(numCells, blocklength*numberBlocks);
stim = zeros(40*80,blocklength*numberBlocks,'uint8');

for blockNum = 1:numberBlocks
    blockNum
%     filename1 = [folderLocation 'WNstim_response_OffParasol_block_' num2str(blockNum) '.mat'];
%     filename1 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/WNresponse/prosthesis_off_parasol/WNstim_response_ProsOFFParasol_block_' num2str(blockNum) '.mat'];
filename1 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/WNresponse/prosthesis_on_parasol/WNstim_response_ProsONParasol_stixds4_short6_' num2str(blockNum) '.mat'];

%     filename1 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/WNresponse/prosthesis_on_parasol/WNstim_response_ProsONParasol_block_' num2str(blockNum) '.mat'];
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
    stimtmp = reshape(matf.whiteNoiseSmall, 40*80,blocklength);
    stim(:,(blockNum-1)*blocklength +1 : blockNum*blocklength) = stimtmp;

    
    % Stimulus here
    % whiteNoiseSmall;
end

% save('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/spikeResp_prosthesis_on_parasol','spikeResp');
% save('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/movie_prosthesis_on_parasol','stim','-v7.3')

save('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/spikeResp_prosthesis_stixds4_short6','spikeResp');
save('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/movie_prosthesis_stixds4_short6','stim','-v7.3');

%figure; imagesc(spikesout*spikesout');

