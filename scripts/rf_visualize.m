

spTimes = find(spikeResp(500,:)==1);

% figure; imagesc(reshape(stim(:,1000),[100 100]));

% % % Zero mean for NS
% for blockNum = 1:floor(size(stim,2)/12000)
%     stim2(:,(blockNum-1)*12000+1:blockNum*12000) = ...
%         uint8(128+127*(double(stim(:,(blockNum-1)*12000+1:blockNum*12000)) - ones(size(stim,1),1)*mean(stim(:,(blockNum-1)*12000+1:blockNum*12000),1)));
% end



rf = (zeros(10000,16));
for spind = 10:10000%length(spTimes)-10;
    imtemp =  double(stim(:,spTimes(spind):spTimes(spind)+15));
    rf = rf + -mean(imtemp(:)) + imtemp;
end

figure;
for fr = 1:31
    subplot(4,4,fr);
    imagesc(reshape(rf(:,fr),[100 100]));
%     caxis([-1e5 1e5]);
end
figure; plot(rf');