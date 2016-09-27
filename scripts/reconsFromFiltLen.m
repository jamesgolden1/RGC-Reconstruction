function [ recons_stim ] = reconsFromFiltLen(filters, spikeResp, numbins)
%reconstruct a stimulus given the filters and spike responses
% numbins=12;
trainTimes = 1:size(spikeResp,2)-numbins;
respTest = zeros(length(trainTimes),size(spikeResp,1)*numbins+1);
respTest(:,1) = ones(length(trainTimes),1);
for t = 1:length(trainTimes)
    starttime = trainTimes(t)+1;
    endtime = trainTimes(t)+numbins;
    respTest(t,2:end) = reshape(spikeResp(:,starttime:endtime)',1,size(spikeResp,1)*numbins);
end

recons_stim = (respTest*filters)';

end

