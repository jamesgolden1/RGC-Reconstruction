
% clear all;
% close all;
disp('Loading testing stimulus data...')
name = 'bar2';
% name = 'grating2';
%set which testing stimulus to use
if strcmp(name,'grating')
    matfOn = matfile('../dat/WNstim_response_OnParasol_36_grating_june10.mat');
    matfOff = matfile('../dat/WNstim_response_OffParasol_64_grating_june10.mat');
    %load ../output/reconstructedSTA/onParasolGratingSTARecon.mat
    %load ../output/reconstructedSTA/offParasolGratingSTARecon.mat
elseif strcmp(name,'bar')
    matfOff = matfile('../dat/WNstim_response_OffParasol_offBig2_64_barGay_june10.mat');
    matfOn = matfile('../dat/WNstim_response_OnParasol_onBig2_36_barGray_june10.mat');
    %load ../output/reconstructedSTA/onParasolBarSTARecon.mat
    %load ../output/reconstructedSTA/offParasolBarSTARecon.mat
elseif strcmp(name,'grating2')
    load('gratingsON.mat')
    matfOn.whiteNoiseSmall = whiteNoiseSmall;
    matfOn.spikesoutsm = spikesoutsm;
    load('gratingsOFF.mat')
    matfOff.whiteNoiseSmall = whiteNoiseSmall;
    matfOff.spikesoutsm = spikesoutsm;
elseif strcmp(name,'bar2')
    load('barON.mat')
    matfOn.whiteNoiseSmall = whiteNoiseSmall;
    matfOn.spikesoutsm = spikesoutsm;
    load('barOFF.mat')
    matfOff.whiteNoiseSmall = whiteNoiseSmall;
    matfOff.spikesoutsm = spikesoutsm;
end



disp('Reconstructing with off parasols...')
%test with off parasols
filtMat_off = matfile('../output/filters_may26_off.mat');
%

%create response matrix for test stimulus

blocklength = 152;
numcells= 64;

spikesout = double(matfOff.spikesoutsm);
spikeRespOff = downSampResp(spikesout, numcells, blocklength);
stim = double(reshape(matfOff.whiteNoiseSmall,96*96,blocklength));
recons_stim_off = reconsFromFilt(filtMat_off.filterMat, spikeRespOff);


mov = reshape(stim,96,96,size(stim,2));
movrecons_off = reshape(recons_stim_off,96,96,size(recons_stim_off,2));


disp('Reconstructing with on parasols...')
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

disp('Reconstructing with both on/off parasols...')
%reconstruct with on and off
spikeRespOnOff =vertcat(spikeRespOn,spikeRespOff);
if exist('../output/filters_may26_on_off1.mat','file') ~= 0
    filtMat_on_off = matfile('../output/filters_may26_on_off.mat');
    recons_stim_on_off = reconsFromFilt(filtMat_on_off.filterMat, spikeRespOnOff);
    filtMatJoint = filtMat_on_off.filterMat;
else
    filtMat_on_off1 = matfile('../output/filters_may26_on_off1.mat');
    filtMat_on_off2 = matfile('../output/filters_may26_on_off2.mat');
    
    filterMatCombined = [filtMat_on_off1.filterMat; filtMat_on_off2.filterMat];
    recons_stim_on_off = reconsFromFilt(filterMatCombined, spikeRespOnOff);
    filtMatJoint = filterMatCombined;
end

mov = reshape(stim,96,96,size(stim,2));
movrecons_on_off = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));


%reconstruct with on trained with on/off

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
save(strcat('../output/movie_may26_',name),'M');

close all;
%%
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
% save(strcat('../output/movie_compare_may26_',name),'M');

figure; hold on;
for it = 1:121
    imagesc(movrecons_on_off(:,:,it));
    
   caxis([min(movrecons_on(:)),max(mov(:))])
   colormap gray
    drawnow
end


%%

load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_may26_on_svd_200.mat')

numcells= 36;
spikesout = double(matfOn.spikesoutsm);
spikeRespOn = downSampResp(spikesout, numcells, blocklength);
% recons_stim_on = reconsFromFilt(filtMat_on.filterMat, spikeRespOn);
recons_stim_on = reconsFromFilt(filterMat, spikeRespOn);

mov = reshape(stim,96,96,size(stim,2));
movrecons_on = reshape(recons_stim_on,96,96,size(recons_stim_on,2));

figure; hold on;
for it = 1:121
    imagesc(movrecons_on(:,:,it));
    
%    caxis([min(movrecons_on(:)),max(mov(:))])
   colormap gray
    drawnow
end


%%

% matfOn

figure; hold on;
for it = 1:121
    imagesc(matfOn.whiteNoiseSmall(:,:,it));
    
%    caxis([min(movrecons_on(:)),max(mov(:))])
   colormap gray
    drawnow
end

%% on only

clear reconsError cc coeffdet
for ishift = 0%:5;
barmov = matfOn.whiteNoiseSmall;
barmovShort = double(barmov(:,:,ishift+[1:122]))/2+.5;

icArray = 1000%[100:100:1000];
for icind = 1:length(icArray)
    
    
     %load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_may26_on_svd_' num2str(icArray(icind)) '.mat'])
        load(['/Users/vision/Documents/MATLAB/james/RGC-Reconstruction/output/svd_reconstruct/on/filters_may26_on_svd_' num2str(icArray(icind)) '.mat'])
    
    numcells= 36;
    spikesout = double(matfOn.spikesoutsm);
    spikeRespOn = downSampResp(spikesout, numcells, blocklength);
    % recons_stim_on = reconsFromFilt(filtMat_on.filterMat, spikeRespOn);
    recons_stim_on = reconsFromFilt(filterMat, spikeRespOn);
    
    mov = reshape(stim,96,96,size(stim,2));
    movrecons_on = reshape(recons_stim_on,96,96,size(recons_stim_on,2));
    
    
    err1 = (barmovShort(:) - movrecons_on(:));
%     figure; hist(err1(:),40)
      reconsError(icind) = sqrt(mean((barmovShort(:) - movrecons_on(:)).^2));
 
   % reconsError(icind) = sum(((barmovShort(:)-mean(barmovShort(:))) - (movrecons_on(:)-mean(movrecons_on(:)))).^2);
    cctemp = corrcoef(barmovShort(:), movrecons_on(:));
    cc(icind) = cctemp(2,1);
    
    coeffdet(icind) = 1 - sum( (barmovShort(:) - movrecons_on(:)).^2 ) / sum ((barmovShort(:) - mean(barmovShort(:))).^2);
    
end


% err1 = (barmovShort(:) - movrecons_on(:));
% figure; hist(err1(:),40)

%
% % sznorm = 96*96*122;
% sznorm = 1;
% % figure; 
% hold on;
% plot(icArray./1000,sqrt(reconsError)/(sznorm),'-o');
% 
% % plot(fracsvd(icArray),sqrt(reconsError)/(sznorm),'-o');
% title('Reconstruction Error for Moving Bar');
% xlabel('Fraction of singular values');
% ylabel(sprintf('Reconstruction error ( sqrt(sum(y - y_0)^2) )'));
% grid on; set(gca,'fontsize',14)


%%

[fracval,fracind] = sort(icArray);

figure; 

set(gcf,'position',[463        1078        1097         260]);
subplot(131);
hold on;
plot(fracsvd(icArray(fracind)),reconsError(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('Mean Squared Error'));
grid on;
set(gca,'fontsize',14);

% figure;
subplot(132);
hold on;
plot(fracsvd(icArray(fracind)),cc(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
% title('Contrast Reversing Grating');
title('Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('Correlation Coefficient'));
grid on;
set(gca,'fontsize',14);

% figure;
subplot(133);
hold on;
plot(fracsvd(icArray(fracind)),coeffdet(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('R^2'));
grid on;
set(gca,'fontsize',14);

end;
%%


    errmov = reshape(recons_stim_on,96,96,122);
    
    
figure; hold on;
for it = 1:121
    imagesc(errmov(:,:,it));
    
%    caxis([min(movrecons_on(:)),max(mov(:))])
   colormap gray
    drawnow
end


%% off only
% figure;
clear reconsError cc coeffdet spikesout
for ishift = 0%:5;
barmov = matfOn.whiteNoiseSmall;
barmovShort = double(barmov(:,:,ishift+[1:122]))/2+.5;

icArray = 1800%[200:200:1800];
for icind = 1:length(icArray)
    
    
     %load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_may26_on_svd_' num2str(icArray(icind)) '.mat'])
        load(['/Users/vision/Documents/MATLAB/james/RGC-Reconstruction/output/svd_reconstruct/off/filters_may26_off_svd_' num2str(icArray(icind)) '.mat'])
    
    numcells= 64;
    spikesout = double(matfOff.spikesoutsm);
    spikeRespOff = downSampResp(spikesout, numcells, blocklength);
    % recons_stim_on = reconsFromFilt(filtMat_on.filterMat, spikeRespOn);
    recons_stim_off = reconsFromFilt(filterMat, spikeRespOff);
    
    mov = reshape(stim,96,96,size(stim,2));
    movrecons_off = reshape(recons_stim_off,96,96,size(recons_stim_off,2));
    
    
    err1 = (barmovShort(:) - movrecons_off(:));
%     figure; hist(err1(:),40)
    
    % reconsError(icind) = sum(((barmovShort(:)-mean(barmovShort(:))) - (movrecons_off(:)-mean(movrecons_off(:)))).^2);
    
  reconsError(icind) = sqrt(mean((barmovShort(:) - movrecons_off(:)).^2));
    cctemp = corrcoef(barmovShort(:), movrecons_off(:));
    cc(icind) = cctemp(2,1);
    
    coeffdet(icind) = 1 - sum( (barmovShort(:) - movrecons_off(:)).^2 ) / sum ((barmovShort(:) - mean(barmovShort(:))).^2);
    
end


% err1 = (barmovShort(:) - movrecons_on(:));
% figure; hist(err1(:),40)

%
% % sznorm = 96*96*122;
% sznorm = 1;
% % figure; 
% hold on;
% plot(icArray./1000,sqrt(reconsError)/(sznorm),'-o');
% 
% % plot(fracsvd(icArray),sqrt(reconsError)/(sznorm),'-o');
% title('Reconstruction Error for Moving Bar');
% xlabel('Fraction of singular values');
% ylabel(sprintf('Reconstruction error ( sqrt(sum(y - y_0)^2) )'));
% grid on; set(gca,'fontsize',14)


%%


[fracval,fracind] = sort(icArray);

figure; 

set(gcf,'position',[463        1078        1097         260]);
subplot(131);
hold on;
plot(fracsvd(icArray(fracind)),reconsError(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('Mean Squared Error'));
grid on;
set(gca,'fontsize',14);

% figure;
subplot(132);
hold on;
plot(fracsvd(icArray(fracind)),cc(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
% title('Contrast Reversing Grating');
title('Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('Correlation Coefficient'));
grid on;
set(gca,'fontsize',14);

% figure;
subplot(133);
hold on;
plot(fracsvd(icArray(fracind)),coeffdet(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('R^2'));
grid on;
set(gca,'fontsize',14);

end;

%% on and off 
figure;

if ~exist('Strain','var')
    load('on_off_covar_svd_strain.mat');
   
    ds = diag(Strain);
    fracsvd = (cumsum(ds./sum(ds(:))));
end

for ishift = 0%:5;
barmov = matfOn.whiteNoiseSmall;
% barmovShort = double(barmov(:,:,ishift+[1:122]));%/2+.5;
barmovShort = double(barmov(:,:,ishift+[1:122]));%/2+.5;

% icArray = [100:100:3000];
% icArray = [100:100:3000 500:25:700 575 584 612 637 663]

% icArray = [100:100:3000 212 225 237 250 262 270 275 254 258 255 256 257 252 253 150];% 275 288 312 325 337 350 363]; 
% icArray = [100:100:300 600:100:3000 212 225 237 250 262 270 275 254 258 255 256 257 252 253 150];% 275 288 312 325 337 350 363]; 
icArray = [100:100:3000]
for icind = 1:length(icArray)
    
    
%     load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_may26_on_off_svd_' num2str(icArray(icind)) '.mat'])
    load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct_hoylab/filters_may26_on_off_svd_' num2str(icArray(icind)) '.mat'])
    numcells= 36;
%     spikesout = double(matfOn.spikesoutsm);
%     spikeRespOn = downSampResp(spikesout, numcells, blocklength);
%     % recons_stim_on = reconsFromFilt(filtMat_on.filterMat, spikeRespOn);
    recons_stim_on_off = reconsFromFilt(filterMat, spikeRespOnOff);
    
    mov = reshape(stim,96,96,size(stim,2));
%     movrecons_on_off = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));
        movrecons_on_off_full = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));
    movrecons_on_off = .5*movrecons_on_off_full./mean(movrecons_on_off_full(:));
%         movrecons_on_off = zeros(size(movrecons_on_off_full));
%     movrecons_on_off(movrecons_on_off_full>0.5) = 1; movrecons_on_off(movrecons_on_off_full<0.5) = 0;
    
    cctemp = corrcoef(barmovShort(:), movrecons_on_off(:));
    cc(icind) = cctemp(2,1);
    
    coeffdet(icind) = 1 - sum( (barmovShort(:) - movrecons_on_off(:)).^2 ) / sum ((barmovShort(:) - mean(barmovShort(:))).^2);
    
    err1 = (barmovShort(:) - movrecons_on_off(:));
%     figure; hist(err1(:),80)
    
  reconsError(icind) = sqrt(mean((barmovShort(:) - movrecons_on_off(:)).^2));
% reconsError(icind) = sum((((barmovShort(:)-mean(barmovShort(:)))./max(1)) - ((movrecons_on_off(:)-mean(movrecons_on_off(:)))./max(1))).^2);
     %reconsError(icind) = sum((((barmovShort(:)-mean(barmovShort(:)))./max(barmovShort(:))) - ((movrecons_on(:)-mean(movrecons_on_off(:)))./max(movrecons_on_off(:)))).^2);
end


% err1 = (barmovShort(:) - movrecons_on(:));
% figure; hist(err1(:),40)

%%

[fracval,fracind] = sort(icArray);

figure; 

% set(gcf,'position',[463        1078        1097         260]);
  set(gcf,'position',[344         913        1068         425]);
subplot(121);
hold on;
plot(fracsvd(icArray(fracind)),reconsError(fracind),'-o');
% plot((icArray(fracind)),reconsError(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('Mean Squared Error'));
grid on;
set(gca,'fontsize',14);

% figure;
subplot(122);
hold on;
plot(fracsvd(icArray(fracind)),cc(fracind),'-o');
% plot((icArray(fracind)),cc(fracind),'-o');
% title('Reconstruction Error for Moving Bar');
title('Contrast Reversing Grating');
% title('Moving Bar');
xlabel('Fraction of singular values');
ylabel(sprintf('Correlation Coefficient'));
grid on;
set(gca,'fontsize',14);

% % figure;
% subplot(133);
% hold on;
% plot(fracsvd(icArray(fracind)),coeffdet(fracind),'-o');
% % plot((icArray(fracind)),coeffdet(fracind),'-o');
% % title('Reconstruction Error for Moving Bar');
% xlabel('Fraction of singular values');
% ylabel(sprintf('R^2'));
% grid on;
% set(gca,'fontsize',14);
end;
%%


    errmov = reshape(recons_stim_on,96,96,122);
    
    
figure; hold on;
for it = 1:121
    imagesc(movrecons_on_off(:,:,it));
    
   caxis([min(movrecons_on_off(:)),max(mov(:))])
   colormap gray
    drawnow
end

%%


figure; plot(diag(Strain)./max(diag(Strain)),'linewidth',3)
hold on; plot(fracsvd,'linewidth',3)

grid on
xlabel('Singular Value Index'); ylabel('Normalized Mangitude');
legend('SVD','Cumulative SVD')
set(gca,'fontsize',14)