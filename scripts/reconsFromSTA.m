function [ recons_stim_out ] = reconsFromSTA(sta, spikeResp)
%reconstruct a stimulus given the filters and spike responses
% numbins=30;
% trainTimes = 1:size(spikeResp,2)-numbins;
% 
% respTest = zeros(length(trainTimes),size(spikeResp,1)*numbins+1);
% respTest(:,1) = ones(length(trainTimes),1);
% 
% for t = 1:length(trainTimes)
%     starttime = trainTimes(t)+1;
%     endtime = trainTimes(t)+numbins;
%     respTest(t,2:end) = reshape(spikeResp(:,starttime:endtime)',1,size(spikeResp,1)*numbins);
% end
% 
% recons_stim = (respTest*filters)';


recons_stim = zeros(96,96,size(spikeResp,2)+30);
recons_ctr = zeros(96,96,size(spikeResp,2)+30);

for cellind = 1:36%size(spikeResp,1)
    clear spTimes
    spTimes = find(spikeResp(cellind,:)==1);
    staRS = reshape(sta{cellind},[96 96 30]);
    for spInd = 1:length(spTimes)
        recons_stim(:,:,spTimes(spInd):spTimes(spInd)+29) = recons_stim(:,:,spTimes(spInd):spTimes(spInd)+29) + staRS(:,:,30:-1:1).*(abs(staRS(:,:,30:-1:1))>.01);
        recons_ctr(:,:,spTimes(spInd):spTimes(spInd)+29)= recons_ctr(:,:,spTimes(spInd):spTimes(spInd)+29) + ones(96,96,30);
    end
    
    
end

recons_stim_out = recons_stim;%./recons_ctr;
% recons_stim_out(isnan(recons_stim_out)) = 0;

end

