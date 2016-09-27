% clear all;
% close all;
disp('Loading filters...')
filtMat_off = matfile('../output/filters_may26_off.mat');
filtMat_on = matfile('../output/filters_may26_on.mat');
filtMat_on_off = matfile('../output/filters_may26_on_off.mat');

disp('Loading movie...')
load ../dat/movie_may26.mat
disp('Loading response matrix and raw spike responses...')
% load ../output/respTest_may26_on_off.mat
load ../output/on_off_respTest.mat
matfON = matfile('../dat/spikeResp_may26_on.mat');
matfOFF = matfile('../dat/spikeResp_may26_off.mat');
spikeResp = vertcat(matfON.spikeResp, matfOFF.spikeResp);

% disp('Loading sta...')
% matfSTA = matfile('../output/staFull_may26_on_off.mat');
% sta = matfSTA.sta;
sta{1} = zeros(96*96,30);

disp('Compute temporal filters for selected subset of cells...')
numbins = 30;
%compare temporal filters
%ttfMat = zeros(numbins,size(8,1));
ttfOn = zeros(numbins,size(4,1));
ttfOff = zeros(numbins,size(4,1));
ttfOnOff = zeros(numbins, size(8,1));
staFrames = cell(size(spikeResp,1),1);
ttfMat = zeros(numbins,size(spikeResp,1));
for icell = 1:size(spikeResp,1)
    icell = 1;
    [m,a] = max(abs(sta{icell}(:,26)));
    ttfMat(:,icell) = sta{icell}(a,:)';
    staFrames{icell} = reshape(sta{icell}(:,26),96,96);
end
ttf = mean(ttfMat,2);

filtMatON = filtMat_on.filterMat;
filtMatOFF = filtMat_off.filterMat;
filtMatONOFF = filtMat_on_off.filterMat;

icelllist = [15,20,23,34,57,64,73,82];
for icell = 1:8
     cellval = icelllist(icell);
     if icell <= 4
        tc = abs(fliplr(filtMatON((cellval-1)*numbins+2:cellval*numbins+1,:)'));
        [m2,a2] = max(tc(:,26));
        ttfOn(:,icell) = fliplr(filtMatON((cellval-1)*numbins+2:(cellval)*numbins+1,a2)');

     else
        tc = abs(fliplr(filtMatOFF((cellval-37)*numbins+2:(cellval-36)*numbins+1,:)'));
        [m2, a2] = max(tc(:,26));
        ttfOff(:,icell-4) = fliplr(filtMatOFF((cellval-37)*numbins+2:(cellval-36)*numbins+1,a2)');

     end
     tc = abs(fliplr(filtMatONOFF((cellval-1)*numbins+2:cellval*numbins+1,:)'));
     [m3, a3] = max(tc(:,26));
     ttfOnOff(:,icell) = fliplr(filtMatONOFF((cellval-1)*numbins+2:(cellval)*numbins+1,a3)');
end




disp('Plotting spatial filters...')
%plot spatial filters
fig = figure('position',[100 100 1024 1024]);
%set(gcf,'visible', 'off')
for i = 1:numbins-3
    for icell = 1:4
        subplot(2,4,icell)
        index = icell;
        %for off cells use icell+4 as index, for on cells use index = icell
        cellval = icelllist(index);
        rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
        imagesc(reshape(filtMatONOFF(rangecols(i),:),96,96));
        %caxis([min(sta{cellval}(:,numbins-i)),max(sta{cellval}(:,numbins-i))]);
        %caxis([0,1]);
        colormap gray
        axis image
        subplot(2,4,icell+4)
        cellval =1;
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
save('../output/presentation/movie_STA_may26_on','M')

disp('Plotting example temporal filters...')
%plot temporal filters
figure;
for i = 1:4
    cellval = icelllist(i);
    rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
    subplot(3,4,i)
    imagesc(reshape(filtMatONOFF(rangecols(4),:),96,96));
    colormap gray
    axis image
    
    subplot(3,4,i+4)
    plot(ttfMat(:,icelllist(i)),'k')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([-0.2,0.3])

    subplot(3,4,i+8)
  
    plot(ttfOn(:,i),'b')
    hold on; plot(ttfOnOff(:,i),'--r')
    legend('Train On', 'Train On/Off', 'Location', 'Best')
    ylim([-0.02,0.06])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);


end
subtitle('ON Cell Temporal Filters vs. STA')
savefig('../output/presentation/temp_filter_may26_on')
close all;
figure;
for i = 5:8
    subplot(3,4,i-4)
    cellval = icelllist(i);
    rangecols = [(cellval-1)*numbins+2 : (cellval)*numbins+1];
    imagesc(reshape(filtMatONOFF(rangecols(4),:),96,96));
    colormap gray
    axis image
    subplot(3,4,i)
    plot(ttfMat(:,icelllist(i)),'k')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([-0.4,0.2])
    subplot(3,4,i+4)
    plot(ttfOff(:,i-4),'b')
    hold on; plot(ttfOnOff(:,i),'--r')
    legend('Train Off', 'Train On/Off', 'Location','Best')
    ylim([-0.08,0.05])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
end
subtitle('OFF Cell Temporal Filters vs. STA')
savefig('../output/presentation/temp_filter_may26_off')
% 
