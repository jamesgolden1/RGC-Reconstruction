clear all;
close all;

name = 'bar';
%set which testing stimulus to use
if strcmp(name,'grating')
    matfOn = matfile('../dat/WNstim_response_OnParasol_36_grating_june10.mat');
    matfOff = matfile('../dat/WNstim_response_OffParasol_64_grating_june10.mat');
    %load ../output/reconstructedSTA/onParasolGratingSTARecon.mat
    %load ../output/reconstructedSTA/offParasolGratingSTARecon.mat
elseif strcmp(name,'bar')
    matfOff = matfile('../dat/WNstim_response_OffParasol_offBig2_64_barGay_june10.mat');
    matfOn = matfile('../dat/WNstim_response_OffParasol_onBig2_36_barGray_june10.mat');
    %load ../output/reconstructedSTA/onParasolBarSTARecon.mat
    %load ../output/reconstructedSTA/offParasolBarSTARecon.mat
end



%test moving bar with off parasols
filtMat_off = matfile('../output/filters_may26_off.mat');
%

%create response matrix for moving bar (first load spikes like in
%loadspikes.m
blocklength = 152;
numcells= 64;

spikesout = double(matfOff.spikesoutsm);
spikeRespOff = downSampResp(spikesout, numcells, blocklength);
stim = double(reshape(matfOff.whiteNoiseSmall,96*96,blocklength));
recons_stim_off = reconsFromFilt(filtMat_off.filterMat, spikeRespOff);


mov = reshape(stim,96,96,size(stim,2));
movrecons_off = reshape(recons_stim_off,96,96,size(recons_stim_off,2));



%reconstruct with ON parasols
filtMat_on = matfile('../output/filters_may26_on.mat');

%create response matrix for moving bar (first load spikes like in
%loadspikes.m

numcells= 36;
spikesout = double(matfOn.spikesoutsm);
spikeRespOn = downSampResp(spikesout, numcells, blocklength);
recons_stim_on = reconsFromFilt(filtMat_on.filterMat, spikeRespOn);


mov = reshape(stim,96,96,size(stim,2));
movrecons_on = reshape(recons_stim_on,96,96,size(recons_stim_on,2));


%reconstruct with on and off
filtMat_on_off = matfile('../output/filters_may26_on_off.mat');
spikeRespOnOff =vertcat(spikeRespOn,spikeRespOff);
recons_stim_on_off = reconsFromFilt(filtMat_on_off.filterMat, spikeRespOnOff);

mov = reshape(stim,96,96,size(stim,2));
movrecons_on_off = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));

%reconstruct with on trained with on/off
filtMatJoint = filtMat_on_off.filterMat;
recons_stim_on_joint = reconsFromFilt(filtMatJoint(1:size(filtMat_on.filterMat,1),:),spikeRespOn);
movrecons_on_joint = reshape(recons_stim_on_joint,96,96,size(recons_stim_on_joint,2));

%reconstruct with off trained with on/off
recons_stim_off_joint = reconsFromFilt(filtMatJoint([1,size(filtMat_on.filterMat,1)+1:end],:),spikeRespOff);
movrecons_off_joint = reshape(recons_stim_off_joint,96,96,size(recons_stim_off_joint,2));


%plot reconstructions with the joint training
fig = figure('position',[100 100 1024 1024]);
set(gcf,'visible', 'off')
for itime=1:121

   subplot(2,4,1)
   imagesc(mov(:,:,itime)); 
   caxis([min(movrecons_off_joint(:)),max(mov(:))])
   title('True Stimulus');
   colormap gray
   axis image


   subplot(2,4,2)
   imagesc(movrecons_off_joint(:,:,itime));
   caxis([min(movrecons_off_joint(:)),max(mov(:))])
   title('Train: ON/OFF, Decode: OFF')
   axis image 

   subplot(2,4,3)
   imagesc(movrecons_on_joint(:,:,itime));
   caxis([min(movrecons_off_joint(:)),max(mov(:))])
   title('Train: ON/OFF, Decode: ON')
   colormap gray
   axis image



   subplot(2,4,4)
   imagesc(movrecons_on_off(:,:,itime));
   caxis([min(movrecons_off_joint(:)),max(mov(:))])
   title('Train: ON/OFF, Decode: ON/OFF');
   colormap gray
   axis image

   M(itime) = getframe(fig);

end
close all;
[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'position', [150 150 w h]);
axis off
save(strcat('../output/presentation/movie_may26_',name),'M');

close all;
%Compare reconstructions trained with only ON and OFF with the joint
%training
fig = figure('position',[100 100 1024 1024]);
set(gcf,'visible', 'off')
for itime=1:121

   subplot(1,5,1)
   imagesc(mov(:,:,itime));
   caxis([min(movrecons_off(:)),max(mov(:))])
   title('True Stimulus');
   colormap gray
   axis image

   subplot(1,5,2)
   imagesc(movrecons_off(:,:,itime));
   caxis([min(movrecons_off(:)),max(mov(:))])
   title('Train: OFF, Decode: OFF');
   colormap gray
   axis image

   subplot(1,5,3)
   imagesc(movrecons_off_joint(:,:,itime));
   caxis([min(movrecons_off(:)),max(mov(:))])
   title('Train: ON/OFF, Decode: OFF')
   axis image 

   subplot(1,5,4)
   imagesc(movrecons_on(:,:,itime));
   caxis([min(movrecons_off(:)),max(mov(:))])
   title('Train: ON, Decode: ON')
   colormap gray
   axis image

   subplot(1,5,5)
   imagesc(movrecons_on_joint(:,:,itime));
   caxis([min(movrecons_off(:)), max(mov(:))])
   title('Train: ON/OFF, Decode: ON')
   colormap gray
   axis image

   M(itime) = getframe(fig);

end
close all;
[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'position', [150 150 w h]);
axis off
save(strcat('../output/presentation/movie_compare_may26_',name),'M');



%%TODO: get correct STA reconstruction from manuallly computed STA and compare with the linear Reconstruction
% old STA reconstruction was using ISETBIO model which was not a fair comparison
%Compare with STA reconstruction

%close all;
%fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
%for itime = 1:121
%    subplot(1,3,1)
%    imagesc(mov(:,:,itime));
%   %plot(mov(:,:,itime));
%    caxis([min(movrecons_on_off(:)),max(mov(:))])
%    title('True Stimulus');
%    colormap gray
%    axis image
%
%    subplot(1,3,2)
%    imagesc(movrecons_on_off(:,:,itime));
%    caxis([min(movrecons_on_off(:)), max(mov(:))])
%    title('Train: ON/OFF, Decode: ON/OFF')
%    colormap gray
%    axis image
%
%    subplot(1,3,3)
%    if strcmp(name, 'grating')
%        imagesc(onParasolGratingSTARecon(1:96,1:96,itime)+offParasolGratingSTARecon(1:96,1:96,itime));
%    elseif strcmp(name, 'bar')
%        imagesc(onParasolBarSTARecon(1:96,1:96,itime)+offParasolBarSTARecon(1:96,1:96,itime));
%    end
%    caxis([min(movrecons_on_off(:)), max(mov(:))])
%    title('Summed STA Reconstruction (ON/OFF)');
%    colormap gray
%    axis image
%    M(itime) = getframe(fig);
%end
%
%close all;
%[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
%hf = figure; 
%% resize figure based on frame's w x h, and place at (150, 150)
%set(hf, 'position', [150 150 w h]);
%axis off
%save(strcat('../output/presentation/movie_compareSTA_may26_',name),'M'); 
