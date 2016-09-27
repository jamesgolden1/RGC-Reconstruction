
disp('Loading filters...')

%%%%%%%%%%%%%%%%%
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_off_parasol_prosthesis_svd_600.mat')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_prosthesis_single_svd_1021.mat')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_prosthesis_single_transpose2_svd_1000.mat')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_short3_svd_1000.mat')
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_short6_svd_800.mat')
numcells=48;

blocklength = 152;
% frstart = 1;
% frend = frstart + blocklength-1;
% spikesout = double(matfOff.spikesoutsm);
% stim = double(reshape(matfOff.whiteNoiseSmall,40*80,blocklength));
% spikeRespOff = downSampRespPhys(spikesout, numcells, blocklength);
% recons_stim_off = reconsFromFilt(filterMat, spikeRespOff);
%   mov_off = reshape(recons_stim_off,[80 40 122]);
%   figure; ieMovie(mov_off);
 %%%%%%%%%%%%%%%%%

% disp('Loading movie...')
% load ../dat/movie_may26.mat
% disp('Loading response matrix and raw spike responses...')
% % load ../output/respTest_may26_on_off.mat
% load ../output/on_off_respTest.mat
% matfON = matfile('../dat/spikeResp_may26_on.mat');
% matfOFF = matfile('../dat/spikeResp_may26_off.mat');
% spikeResp = vertcat(matfON.spikeResp, matfOFF.spikeResp);

% disp('Loading sta...')
% matfSTA = matfile('../output/staFull_may26_on_off.mat');
% sta = matfSTA.sta;
% % sta{1} = zeros(96*96,30);
% 
% disp('Compute temporal filters for selected subset of cells...')
numbins = 20;
% %compare temporal filters
% %ttfMat = zeros(numbins,size(8,1));
% ttfOn = zeros(numbins,size(4,1));
% ttfOff = zeros(numbins,size(4,1));
% ttfOnOff = zeros(numbins, size(8,1));
% staFrames = cell(size(spikeResp,1),1);
% ttfMat = zeros(numbins,size(spikeResp,1));
% for icell = 1:size(spikeResp,1)
%     [m,a] = max(abs(sta{icell}(:,26)));
%     ttfMat(:,icell) = sta{icell}(a,:)';
%     staFrames{icell} = reshape(sta{icell}(:,26),96,96);
% end
% ttf = mean(ttfMat,2);

filtMatON = filterMat;%filtMat_on.filterMat;
filtMatONSVD = filterMat;%filtMat_on.filterMat;
% filtMatOFF = filtMat_off.filterMat;
% filtMatONOFF = filtMat_on_off.filterMat;

icelllist = [1:48];%,57,64,73,82];
for icell = 1:48
     cellval = icelllist(icell);
%      if icell <= 4
        tc = abs(fliplr(filtMatON((cellval-1)*numbins+2:cellval*numbins+1,:)'));
        [m2,a2] = max(tc(:,26));
        ttfOn(:,icell) = fliplr(filtMatON((cellval-1)*numbins+2:(cellval)*numbins+1,a2)');
%      end
%      else
%         tc = abs(fliplr(filtMatOFF((cellval-37)*numbins+2:(cellval-36)*numbins+1,:)'));
%         [m2, a2] = max(tc(:,26));
%         ttfOff(:,icell-4) = fliplr(filtMatOFF((cellval-37)*numbins+2:(cellval-36)*numbins+1,a2)');
% 
%      end
%      tc = abs(fliplr(filtMatONOFF((cellval-1)*numbins+2:cellval*numbins+1,:)'));
%      [m3, a3] = max(tc(:,26));
%      ttfOnOff(:,icell) = fliplr(filtMatONOFF((cellval-1)*numbins+2:(cellval)*numbins+1,a3)');
end
figure; plot(ttfOn)
%%
% icOn = 300; icOff = 400; icOnOff = 400;
% load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_may26_on_svd_' num2str(icOn) '.mat'])
% filtMatONSVD = filterMat;
% load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_may26_off_svd_' num2str(icOff) '.mat'])
% filtMatOFFSVD = filterMat;
% load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct_hoylab/filters_may26_on_off_svd_' num2str(icOnOff) '.mat'])
% filtMatONOFFSVD = filterMat;
% 
% icelllist = [15,20,23,34,57,64,73,82];
% for icell = 1:8
%      cellval = icelllist(icell);
%      if icell <= 4
%         tc = abs(fliplr(filtMatONSVD((cellval-1)*numbins+2:cellval*numbins+1,:)'));
%         [m2,a2] = max(tc(:,26));
%         ttfOnSVD(:,icell) = fliplr(filtMatONSVD((cellval-1)*numbins+2:(cellval)*numbins+1,a2)');
% 
%      else
%         tc = abs(fliplr(filtMatOFFSVD((cellval-37)*numbins+2:(cellval-36)*numbins+1,:)'));
%         [m2, a2] = max(tc(:,26));
%         ttfOffSVD(:,icell-4) = fliplr(filtMatOFFSVD((cellval-37)*numbins+2:(cellval-36)*numbins+1,a2)');
% 
%      end
%      tc = abs(fliplr(filtMatONOFFSVD((cellval-1)*numbins+2:cellval*numbins+1,:)'));
%      [m3, a3] = max(tc(:,26));
%      ttfOnOffSVD(:,icell) = fliplr(filtMatONOFFSVD((cellval-1)*numbins+2:(cellval)*numbins+1,a3)');
% end

%%

disp('Plotting spatial filters...')
%plot spatial filters
% fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
figure;

icelllist = [1:14];%,57,64,73,82];

    for icell = 1:34%:4
        for i = 1:numbins-3
%         subplot(3,4,icell)
        index = icell+0;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONSVD(rangecols(i),:),40,80));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        drawnow
%         subplot(3,4,icell+4)
%         index = icell+0;
%         %for off cells use icell+4 as index, for on cells use index = icell
%         cellval = icelllist(index);
%         rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
%         imagesc(reshape(filtMatONOFF(rangecols(i),:),96,96));
%         %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
%         %caxis([0,1]);
%         colormap gray
%         axis image
%         
%         subplot(3,4,icell+8)
%         cellval =icelllist(index);
%         imagesc(reshape(sta{cellval}(:,numbins-i),96,96));
%         %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
%         colormap gray
%         axis image
    end
%     subtitle('Example ON Cell Learned Spatial Filters vs. STA');
%     M(i) = getframe(fig);
end
%close all;
%[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
%hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
%set(hf, 'position', [150 150 w h]);
%axis off

% save('../output/presentation/movie_STA_may26_on','M')


%%
figure;
    for icell = 1%1:5%:4
        
        sta = reshape(filterMat((icell-1)*numbins+2:(icell)*numbins+1,:)',40,80,numbins);
%         maxsta = max(sta(:));
        ieMovie(sta);
    end
%     figure; plot(RGB2XWFormat(sta)');

%%


    %%
    
    figure; 
    for ci = 1:48
        subplot(7,7,ci)
    
    imagesc(reshape(filterMat((ci-1)*numbins+2,:),40,80));% caxis([-.002 .002])
    end
%%

disp('Plotting spatial filters...')
%plot spatial filters
fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
for i = 1:numbins-3
    for icell = 1:4
        subplot(3,4,icell)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFFSVD(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+4)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFF(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+8)
        cellval =icelllist(index);
        imagesc(reshape(sta{cellval}(:,numbins-i),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        colormap gray
        axis image
    end
    subtitle('Example ON Cell Learned Spatial Filters vs. STA');
    M(i) = getframe(fig);
end
%close all;
%[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
%hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
%set(hf, 'position', [150 150 w h]);
%axis off

% save('../output/presentation/movie_STA_may26_on','M')

%%
disp('Plotting example temporal filters...')
%plot temporal filters
figure;
for i = 1:4
%     cellval = icelllist(i);
%     rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
%     subplot(3,4,i)
%     imagesc(reshape(filtMatONOFF(rangecols(4),:),96,96));
%     colormap gray
%     axis image

    subplot(3,4,i)
    plot(8*ttfOnSVD(:,i),'b')
    hold on; plot(8*ttfOnOffSVD(:,i),'--r')    
    legend('Train On', 'Train On/Off', 'Location', 'Best')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([-0.2,0.3])
    
    subplot(3,4,i+8)
    plot(ttfMat(:,icelllist(i)),'k')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([-0.2,0.3])

    subplot(3,4,i+4)  
    plot(ttfOn(:,i),'b')
    hold on; plot(ttfOnOff(:,i),'--r')
    legend('Train On', 'Train On/Off', 'Location', 'Best')
    ylim([-0.02,0.06])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);


end
suptitle('ON Cell Temporal Filters vs. STA')
% savefig('../output/presentation/temp_filter_may26_on')
% close all;


%%
disp('Plotting example temporal filters...')
%plot temporal filters
figure;
for i = 1:4
%     cellval = icelllist(i);
%     rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
%     subplot(3,4,i)
%     imagesc(reshape(filtMatONOFF(rangecols(4),:),96,96));
%     colormap gray
%     axis image

    subplot(3,4,i)
    plot(8*ttfOffSVD(:,i),'b')
    hold on; plot(8*ttfOnOffSVD(:,i),'--r')    
    legend('Train On', 'Train On/Off', 'Location', 'Best')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([-0.2,0.3])
    
    subplot(3,4,i+8)
    plot(ttfMat(:,icelllist(i)),'k')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([-0.2,0.3])

    subplot(3,4,i+4)  
    plot(ttfOff(:,i),'b')
    hold on; plot(ttfOnOff(:,i),'--r')
    legend('Train On', 'Train On/Off', 'Location', 'Best')
    ylim([-0.02,0.06])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);


end
suptitle('OFF Cell Temporal Filters vs. STA')
% savefig('../output/presentation/temp_filter_may26_on')
% close all;

%%
figure;
for i = 5:8
%     subplot(3,4,i-4)
%     cellval = icelllist(i);
%     rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
%     imagesc(reshape(filtMatONOFFSVD(rangecols(4),:),96,96));
%     colormap gray
%     axis image
    
    subplot(3,4,i+4)
    plot(ttfMat(:,icelllist(i)),'k')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([-0.4,0.2])
    
    subplot(3,4,i)
    plot(ttfOff(:,i-4),'b')
    hold on; plot(ttfOnOff(:,i),'--r')
    legend('Train Off', 'Train On/Off', 'Location','Best')
    ylim([-0.08,0.05])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
        subplot(3,4,i-4)
    plot(1*ttfOffSVD(:,i-4),'b')
    hold on; plot(2*ttfOnOffSVD(:,i),'--r')
    legend('Train Off', 'Train On/Off', 'Location','Best')
    ylim([-0.08,0.05])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
end
suptitle('OFF Cell Temporal Filters vs. STA')
% savefig('../output/presentation/temp_filter_may26_off')
% 

%%
%%%%%%%%%%%%%%%%%%%

%%
disp('Plotting spatial filters...')
%plot spatial filters
fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
for i = 1:numbins-3
    for icell = 1:4
        subplot(3,4,icell)
        index = icell+0;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONSVD(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+4)
        index = icell+0;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatON(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+8)
        cellval =icelllist(index);
        imagesc(reshape(sta{cellval}(:,numbins-i),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        colormap gray
        axis image
    end
    subtitle('Example ON Cell Learned Spatial Filters vs. STA');
    M(i) = getframe(fig);
end


%%
disp('Plotting spatial filters...')
%plot spatial filters
fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
for i = 1:numbins-3
    for icell = 1:4
        subplot(3,4,icell)
        index = icell+0;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFFSVD(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+4)
        index = icell+0;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONSVD(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+8)
        cellval =icelllist(index);
        imagesc(reshape(sta{cellval}(:,numbins-i),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        colormap gray
        axis image
    end
    subtitle('Example ON Cell Learned Spatial Filters vs. STA');
    M(i) = getframe(fig);
end

%%
disp('Plotting spatial filters...')
%plot spatial filters
fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
for i = 1:numbins-3
    for icell = 1:4
        subplot(3,4,icell)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
%         cellval = icelllist(index)-36;
%         rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
%         imagesc(reshape(filtMatOFFSVD(rangecols(i),:),96,96));
        cellval = icelllist(index)-0;
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFFSVD(rangecols(i),:),96,96));

        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+4)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index)-36;
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatOFFSVD(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+8)
        cellval =icelllist(index+0)
        imagesc(reshape(sta{cellval}(:,numbins-i),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        colormap gray
        axis image
    end
    subtitle('Example ON Cell Learned Spatial Filters vs. STA');
    M(i) = getframe(fig);
end
%close all;
%[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
%hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
%set(hf, 'position', [150 150 w h]);
%axis off

% save('../output/presentation/movie_STA_may26_on','M')

%% Flip through
disp('Plotting spatial filters...')
%plot spatial filters
% fig = figure('position',[100 100 1024 1024]);
hold on;
%set(gcf,'visible', 'off')

icelllist = [15,20,23,34,57,64,73,82];
% 
% % for listind = 1:16
% % listind = 0;
listind = listind+1;
    icelllist = 36 + 4*(listind-1) + [1:8];
%     icelllist(5:8) = [56 57 58 63];
%     icelllist(5:8) = [64 66 72 73];
    icelllist(5:8) = [74 79 83 89];
% 56 57 58 63 64 66 72 73 74 83 89
for i = 1%:numbins-3
    for icell = 1:4
        subplot(3,4,icell)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
%         cellval = icelllist(index)-36;
%         rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
%         imagesc(reshape(filtMatOFFSVD(rangecols(i),:),96,96));
        cellval = icelllist(index)-0;
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFFSVD(rangecols(i),:),96,96));

        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+4)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index)-36;
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatOFFSVD(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+8)
        cellval =icelllist(index+0);
        imagesc(reshape(sta{cellval}(:,numbins-i),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        colormap gray
        axis image
    end
    suptitle(sprintf('%d %d %d %d',icelllist(5:8)));
%     subtitle('Example ON Cell Learned Spatial Filters vs. STA');
%     M(i) = getframe(fig);
end
% end
%%

name_str = 'on_filters.mp4';
path_str = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/'];
    
vObj = VideoWriter([path_str name_str],'MPEG-4');
vObj.FrameRate = 30;
vObj.Quality = 100;
open(vObj);

disp('Plotting spatial filters...')
%plot spatial filters
fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
for i = 1:numbins-3
    for icell = 1:4
        subplot(3,4,icell)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFFSVD(rangecols(i),:),96,96)./max(filtMatONOFFSVD(:)));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+4)
        index = icell+4;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFF(rangecols(i),:),96,96)./max(filtMatONOFF(:)));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        
        subplot(3,4,icell+8)
        cellval =icelllist(index);
        imagesc(reshape(sta{cellval}(:,numbins-i),96,96)./max(sta{cellval}(:)));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        colormap gray
        axis image
    end
%     suptitle('Example ON Cell Learned Spatial Filters vs. STA');
%     M(i) = getframe(fig);
    
    
    F = getframe(fig);
    writeVideo(vObj,F);
    
end

close(vObj);
%close all;
%[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
%hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
%set(hf, 'position', [150 150 w h]);
%axis off

% save('../output/presentation/movie_STA_may26_on','M')

