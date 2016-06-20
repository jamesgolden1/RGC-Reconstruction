load filtersRF.mat
load respTest.mat
load spikeResp.mat
load movie.mat
batchsize = 96;
sizeStim = size(stim,2);
trainSize = 0.6;
validSize = 0.2;

%indices to block out train and validation data from the full stimulus
i1 = floor(trainSize*sizeStim);
i2 = floor((trainSize+validSize)*sizeStim);

%number of time bins used in reconstruction from single-timepoint
numbins = 30;

inds = find(spikeResp(25,:) > 0);
testTimes = i1+1:i2-numbins;
indstest = find(ismember(testTimes, inds)==1);
recons_test = zeros(size(stim,1), length(indstest(1:100)), 'uint8');

recons_test = (respTest(indstest(1:100),1:30)*filterMat(1:30,:))';

load staFull.mat
staFrames = cell(size(spikeResp,1),1);
ttf = zeros(30,1);
ttfMat = zeros(30,size(spikeResp,1));
for icell = 1:size(spikeResp,1)
    [m,a] = max(sta{icell}(:,26));
    ttfMat(:,icell) = sta{icell}(a,:)';
    staFrames{icell} = reshape(sta{icell}(:,26),96,96);
end

ttf = mean(ttfMat,2);
clear ttfMat
movoutrecons = temporal_filter(recons_test, ttf,0.7);
movoutrecons = reshape(movoutrecons,96,96,100);
movoutrecons = spatial_filter(movoutrecons, staFrames, 0.1,0.5,96,96);
close all;
for itime=10:40
  % subplot(2,1,1)
  % imagesc(movoutorig(20:30,1:10,10+itime));
   %caxis([min(movoutorig(:)),max(movoutorig(:))])
   %colormap gray
   %axis image
   %subplot(2,1,2)
   imagesc(movoutrecons(:,:,itime));
   %caxis([min(movoutrecons(:)),max(movoutrecons(:))])
   colormap gray
   axis image
   pause
end