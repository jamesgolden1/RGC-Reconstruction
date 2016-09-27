function [ spikeResp] = downSampRespSetLength(spikesout, numcells)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    blocksize = 100;
    blocklength = round(size(spikesout,2)/blocksize);
spikeResp = zeros(numcells, blocklength);
for i = 1:blocklength

    endval  = i*blocksize;
    if endval > size(spikesout,2)
        endval = size(spikesout,2);
    end
    startval = (i-1)*blocksize+1;
    spikeResp(:,i) = sum(spikesout(:,startval:endval),2);
end

end

