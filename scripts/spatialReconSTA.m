
numcells = 48;
for i=1:numcells;
    c1(i,:) = innerRetina.mosaic{1}.cellLocation{i};
end;

figure; scatter(c1(:,1),c1(:,2));

blocklength = 352;

m1 = max(c1(:,1)); m2 = max(c1(:,2));
spikesMov = zeros(m1,m2, blocklength);

spikeRespOff = downSampRespPhys(spikesoutsm, numcells, blocklength);
% spikeRespOff = downSampResp(spikesoutsm, numcells, blocklength);
l1 = innerRetina.mosaic{1}.get('responselinear');
for i=1:numcells;
%     spikesMov(c1(i,1),c1(i,2),1:blocklength) = spikeRespOff(i,:);
    spikesMov(c1(i,1),c1(i,2),1:blocklength) = l1(i,1,:);
end

spikesMov = permute(spikesMov,[2 1 3]);

figure; ieMovie(spikesMov);

%%


numDownSampleBlocks = ceil(size(spikesoutsm,2)*0.1);
val = zeros(1, numDownSampleBlocks);
lenBlock = round(1./.1);
whiteNoiseSmall = uint8(squeeze(testmovieRS(:,:,:,1)));
for ii=1:48%:nCells(1)
    for jj = 1%:nCells(2)
        for kk = 1:numDownSampleBlocks-1
            downSampleIndStart = (kk-1)*lenBlock + 1;
            downSampleIndEnd   = (kk)*lenBlock;
            val(ii,kk) = sum(spikesoutsm(ii,downSampleIndStart:downSampleIndEnd));
        end
        kk = numDownSampleBlocks;downSampleIndStart = (kk-1)*lenBlock + 1;
        val(ii,kk) = sum(spikesoutsm(ii,jj,downSampleIndStart:end));
    end
end
            
figure;
for ii = 1%:48
    clear evIdx
    evIdx = find(val(ii,:)>0);
    rf = zeros(80,40,30);
    
    for evind = 1:length(evIdx)-4
        rf = rf + double((whiteNoiseSmall(:,:,evIdx(evind):evIdx(evind)+29)))-.5;
    end
    
    rf2 = RGB2XWFormat(rf);
    hold on;
%     figure; 
    plot(mean(rf2,1))
end
figure; ieMovie(rf);

rf2 = RGB2XWFormat(rf);
figure; plot(mean(rf2,1))
  

%%

% spikesoutsm = reshape(spikesoutsm,36,11999912);

numDownSampleBlocks = ceil(size(spikesoutsm,2)*0.01);
val = zeros(1, numDownSampleBlocks);
lenBlock = round(1./.01);
% whiteNoiseSmall = uint8(squeeze(testmovieRS(:,:,:,1)));
for ii=1:12%:nCells(1)
    for jj = 1%:nCells(2)
        for kk = 1:numDownSampleBlocks-1
            downSampleIndStart = (kk-1)*lenBlock + 1;
            downSampleIndEnd   = (kk)*lenBlock;
            val(ii,kk) = sum(spikesoutsm(ii,downSampleIndStart:downSampleIndEnd));
        end
        kk = numDownSampleBlocks;downSampleIndStart = (kk-1)*lenBlock + 1;
        val(ii,kk) = sum(spikesoutsm(ii,jj,downSampleIndStart:end));
    end
end
            
figure;
for ii = 8%:48
    clear evIdx
    evIdx = find(val(ii,:)>0);
    rf = zeros(40,80,30);
    
    for evind = 1:length(evIdx)-4
        rf = rf + double((whiteNoiseSmall(:,:,evIdx(evind):evIdx(evind)+29)))-.5;
    end
    
    rf2 = RGB2XWFormat(rf);
    hold on;
%     figure; 
    plot(mean(rf2,1))
end
figure; ieMovie(rf);

rf2 = RGB2XWFormat(rf);
figure; plot(mean(rf2,1))

%%

spikesoutsm = matfON.spikeResp;
% spikesoutsm = spikesout;
figure;
for ii = 3%:14%:48
    clear evIdx
evIdx = find(spikesoutsm(ii,:)>0);
rf0 = zeros(80*40,30);

for evind = 1:length(evIdx)-5
    rf0 = rf0 + spikesoutsm(ii,evIdx(evind))*double((stim(:,evIdx(evind):evIdx(evind)+29)))-.5;
end

%     rf2 = RGB2XWFormat(rf);
    hold on;
%     figure; 
    plot(mean(rf0,1))
end

rf = reshape(rf0,40,80,30);
  figure; ieMovie(rf);
  
%   rf2 = RGB2XWFormat(rf);
%   figure; plot(mean(rf2,1))

%%
figure;
for ii = 1:48
    clear evIdx
evIdx = find(spikesoutsm(ii,:)>0);
rf0 = zeros(80*40,30);

for evind = 1:length(evIdx)-10
    rf0 = rf0 + spikesoutsm(ii,evIdx(evind))*double((stim(:,evIdx(evind):evIdx(evind)+29)))-.5;
end

subplot(7,7,ii);
imagesc(reshape(rf0(:,18),40,80)); drawnow

end

%%

figure;
for ii = 1%:48
    clear evIdx
% evIdx = find(respTrain(:,2+30)>0);
evIdx = find(respTrain(2+30,:)>0);
rf0 = zeros(80*40,30);

for evind = 1:length(evIdx)-5
    rf0 = rf0 + double((stim(:,evIdx(evind):evIdx(evind)+29)))-.5;
end

%     rf2 = RGB2XWFormat(rf);
    hold on;
%     figure; 
    plot(mean(rf0,1))
end

rf = reshape(rf0,40,80,30);
  figure; ieMovie(rf);
%%

sRF= innerRetina.mosaic{1}.get('sRFcenter');
tC= innerRetina.mosaic{1}.get('tCenter');

% sRF= ir.mosaic{1}.get('sRFcenter');
% tC= ir.mosaic{1}.get('tCenter');
% sRF= obj.get('sRFcenter');
% tC= obj.get('tCenter');
figure;
for ci = 1:length(sRF)
    clear outerprodtemp maxv maxi
    sRFmat(ci,:) = sRF{ci}(:);
    tmat(ci,:) = tC{ci};
    outerprodtemp = sRFmat(ci,:)'*tmat(ci,:);
     outerprod(ci,:,:) =  outerprodtemp;
     [maxv,maxi] = max(abs(outerprodtemp(:)));
     sv(ci) = sign(outerprodtemp(maxi))
     hold on;
     plot(mean(outerprodtemp,1));
end

%%

% sRF= innerRetina.mosaic{1}.get('sRFcenter');
% tC= innerRetina.mosaic{1}.get('tCenter');
figure;
for ci = 1:length(mosaicGLM)
    clear outerprodtemp maxv maxi
    sRFmat(ci,:) = mosaicGLM{ci}.linearfilters.Stimulus.space_rk1(:);
    tmat(ci,:) = mosaicGLM{ci}.linearfilters.Stimulus.time_rk1(:);
    outerprodtemp = sRFmat(ci,:)'*tmat(ci,:);
%      outerprod(ci,:,:) =  outerprodtemp; 
     [maxv,maxi] = max(abs(outerprodtemp(:)));
     sv(ci) = sign(outerprodtemp(maxi))
     outerprodtemp(maxi)
     hold on;
     plot(mean(outerprodtemp,1));
end
