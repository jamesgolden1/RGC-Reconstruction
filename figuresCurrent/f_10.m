
% f_10
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 10: plot frames from reconstructed movies.
%
% 2018 JRG (c) isetbio team

%%
clear
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/mov_hallway_mosaic0.mat')

figure; %set(gcf,'position',[440   453   961   345]);
set(gcf,'position',[440    62   620   736]);
% ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% ha = tight_subplot(4,3,[.01 .03],[.1 .01],[.01 .01]);
% subplot(131);
% axes(ha(1));

fi = 0;
figure;
imagesc(reshape(stim(:,150),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi) '.pdf']);

% subplot(132);
% axes(ha(2));
figure;
imagesc(reshape(stim(:,250),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+1) '.pdf']);

% subplot(133);
% axes(ha(3));
figure;
imagesc(reshape(stim(:,400),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+2) '.pdf']);

%%
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/hallway_movRecon_oct10_healthy.mat')
fi = 3;
% figure; set(gcf,'position',[440   453   961   345]);
% ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% ha = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]);
% subplot(131);
% axes(ha(1+fi));
figure;
imagesc(reshape(movRecon(:,16+150),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi) '.pdf']);

% subplot(132);
% axes(ha(2+fi));
figure;
imagesc(reshape(movRecon(:,16+250),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+1) '.pdf']);

% subplot(133);
% axes(ha(3+fi));
figure;
imagesc(reshape(movRecon(:,16+400),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+2) '.pdf']);

%%
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/hallway_movRecon_oct10_learning.mat')
fi = 6;
% figure; set(gcf,'position',[440   453   961   345]);
% ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% ha = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]);
% subplot(131);
% axes(ha(1+fi));
figure;
imagesc(reshape(movRecon(:,5+150),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi) '.pdf']);

% subplot(132);
% axes(ha(2+fi));
figure;
imagesc(reshape(movRecon(:,5+250),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+1) '.pdf']);

% subplot(133);
% axes(ha(3+fi));
figure;
imagesc(reshape(movRecon(:,5+400),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+2) '.pdf']);



%%
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/hallway_movRecon_oct10_nolearning_on.mat')
fi = 9;
% figure; set(gcf,'position',[440   453   961   345]);
% ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% ha = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]);
% subplot(131);
% axes(ha(1+fi));
figure;
imagesc(reshape(movRecon(:,5+150),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi) '.pdf']);

% subplot(132);
% axes(ha(2+fi));
figure;
imagesc(reshape(movRecon(:,5+250),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+1) '.pdf']);

% subplot(133);
% axes(ha(3+fi));
figure;
imagesc(reshape(movRecon(:,5+400),[100 100]))
colormap gray; axis equal;
axis off;
print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+2) '.pdf']);

% print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames_jan30.pdf']);

%%
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/_hallway_aug22_.mat')
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/_rhallway_aug22_.mat')
filterMat = load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/prosthesis_70_training_aug13/filtersmosaic0_sv 5_w1_sh3_dr30_pitch_70_decay_2.mat')

spikeAug(1,:) = ones(1,size(spikeResp,2));
spikeAug(2:9716+1,:) = spikeResp;
% load('filters__mosaic0.mat')


if isstruct(filterMat)
    droputIndices = filterMat.dropoutIndices;
    filterMatSm = filterMat.filterMat;
    filterMat=zeros(9717,10000);
    filterMat(droputIndices,:)=filterMatSm;
    filterMatSm=[];
    
    %         spikeAugFull = spikeAug;
    %         spikeAug = zeros(size(spikeAug));
    %         spikeAug(droputIndices,:) = spikeAugFull(droputIndices,:);
end
% filterMat2 = filterMat; filterMat = [];

movRecon = filterMat'*spikeAug;
%%
fi = 12;
% figure; set(gcf,'position',[440   453   961   345]);
% ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% ha = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]);
% subplot(131);
% axes(ha(1+fi));
figure;
imagesc(reshape(movRecon(:,5+150),[100 100]))
colormap gray; axis equal;
axis off;
% print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi) '.pdf']);
print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/figures_aug2018/dynamic_frames' num2str(fi) '.pdf']);

% subplot(132);
% axes(ha(2+fi));
figure;
imagesc(reshape(movRecon(:,5+250),[100 100]))
colormap gray; axis equal;
axis off;
% print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+1) '.pdf']);

print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/figures_aug2018/dynamic_frames' num2str(fi+1) '.pdf']);
%
% subplot(133);
% axes(ha(3+fi));
figure;
imagesc(reshape(movRecon(:,5+400),[100 100]))
colormap gray; axis equal;
axis off;
% print(gcf, '-dpdf',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/oct_movies/jan30/dynamic_frames' num2str(fi+2) '.pdf']);

print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/figures_aug2018/dynamic_frames' num2str(fi+2) '.pdf']);