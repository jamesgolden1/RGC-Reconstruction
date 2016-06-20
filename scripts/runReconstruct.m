clear all;
close all;

%%run Linear reconstruction for on, off, and joint on/off training 
windowsize = 30;
disp('Loading ON parasol spike responses...')
matfON = matfile('../dat/spikeResp_may26_on.mat');
disp('Loading OFF parasol spike responses')
matfOFF = matfile('../dat/spikeResp_may26_off.mat');
%load ON parasols
disp('Loading stimulus movie...')
load ../dat/movie_may26.mat

movielength = size(stim,2);
numOnCells = size(matfON.spikeResp,1);
numOffCells = size(matfOFF.spikeResp,1);


disp(['Total Movie Length in Frames: ' num2str(movielength)]);
disp(['Number of ON Cells: ' num2str(numOnCells)])
disp(['Number of OFF Cells: ' num2str(numOffCells)])


disp('Training ON only filters...')
fileext = 'may26_on';
[recons_test, filters]=linearReconstruct(stim,matfON.spikeResp,fileext, windowsize);

disp('Training OFF only filters...')
%load OFF parasols
fileext = 'may26_off';
[recons_test, filters]=linearReconstruct(stim,matfOFF.spikeResp,fileext, windowsize);

disp('Training joint ON/OFF filters...')
%create ON/OFF spikeResp
fileext = 'may26_on_off';
spikeResp1 = vertcat(matfON.spikeResp, matfOFF.spikeResp);
[recons_test, filters]=linearReconstruct(stim,spikeResp1,fileext, windowsize);


disp('Calculating STAs for comparison...')
%%calculate STAs
spikeResp = spikeResp1;
sta = cell(size(spikeResp,1),1);
for i = 1:size(spikeResp,1)
  disp(['Cell: ' num2str(i)])
  evIdx = find(spikeResp(i,:));
  avg = zeros(size(stim,1),windowsize);
  sumspike = 0;
  for w = 1:20000
      if evIdx(w) > windowsize
          sumspike = sumspike+spikeResp(i,evIdx(w));
          wIdx = evIdx(w) - windowsize : evIdx(w)-1;
          stimdoub = double(stim(:,wIdx));
          avg = avg+double(stimdoub)*spikeResp(i,evIdx(w));
      end
  end
  avg = avg./sumspike;
  sta{i} = avg-0.5;
end
save(strcat('../output/staFull_',fileext),'sta');

%%apply spatial and temporal filtering to white noise movies
% load staFull_may26_on_off.mat
% staFrames = cell(size(spikeResp,1),1);
% ttf = zeros(30,1);
% ttfMat = zeros(30,size(spikeResp,1));
% for icell = 1:size(spikeResp,1)
%     [m,a] = max(sta{icell}(:,26));
%     ttfMat(:,icell) = sta{icell}(a,:)';
%     staFrames{icell} = reshape(sta{icell}(:,26),96,96);
% end
% 
% ttf = mean(ttfMat,2);
% clear ttfMat
% % 
% %mov = reshape(stim,96,96,size(stim,2));
% %movrecons = reshape(recons_stim,96,96,size(recons_stim,2));
% sizeStim = size(stim,2);
% trainSize = 0.6;
% validSize = 0.2;
% numbins=30;
% i1 = floor(trainSize*sizeStim);
% i2 = floor((trainSize+validSize)*sizeStim);
% testTimes = i1+1:i2-numbins;

% movoutrecons = temporal_filter(double(recons_test(:,1:10000)), ttf,0.7);
% movoutorig = temporal_filter(double(stim(:,testTimes(1:10000))), ttf, 0.7);

% movoutrecons = reshape(movoutrecons,96,96,10000);
% movoutorig = reshape(movoutorig,96,96,10000);

% movoutorig = spatial_filter(movoutorig, staFrames, 0.1,0.5,96,96);
% movoutrecons = spatial_filter(movoutrecons, staFrames, 0.1,0.5,96,96);
% save('../output/movout_may30_filt_off', 'movoutorig');
% save('../output/movoutrecons_may30_filt_off', 'movoutrecons');

%%visualize white noise reconstructions
% for itime=1000:1030
%    subplot(2,1,1)
%    imagesc(movoutorig(:,:,itime));
%    %imagesc(reshape(stim(:,itime),[96,96]));
%    %caxis([min(movoutorig(:)),max(movoutorig(:))])
%    colormap gray
%    axis image
%    subplot(2,1,2)
%    imagesc(movoutrecons(:,:,itime));
%    %imagesc(reshape(recons_test(:,itime),[96,96]));
%    %caxis([min(movoutrecons(:)),max(movoutrecons(:))])
%    colormap gray
%    axis image
%    pause
% end
% % 

%%plot 2d correlation of reconstruction and true stimulus
% corr2d =zeros(96,96);
% for ix=1:96
%     for iy=1:96
%     corr2d(ix,iy) = corr(squeeze(movoutorig(ix,iy,:)),squeeze(movoutrecons(ix,iy,:)));
%     end
% end
% figure;imagesc(corr2d);colormap gray


