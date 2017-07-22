

% addpath('/Volumes/Lab/Users/james/matlab/private/colleen/Motion Computation Code/Figure Code/'); % tight_subplot

% load('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/june5wn/filters_mosaic0_sv50_w1_sh6.mat')
% filterWN = filterMat; clear filterMat

% load('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/june5wn/sta.mat')

load('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/may22/filters_mosaic0_sv60_w1_sh4.mat')
filterNS = filterMat; clear filterMat

% load('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/may22/sta.mat')

load('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/june16prima/filters_wmean/filters_mosaic0_sv 5_w1_sh2_tr80.mat')
filterPros = filterMat;

%% Get max location of each filter for zooming
[mgr,mgc] = meshgrid(1:100,1:100);

[cmax,cind] = max(abs(filterNS),[],2);
[fmaxc,fmaxr] = ind2sub([100 100],cind);

mgrmat = mgr(:)*ones(1,size(fmaxr,1));
fmaxrmat = ones(size(mgrmat,1),1)*fmaxr';
mgrd = ((mgrmat - fmaxrmat)').^2;

mgcmat = mgc(:)*ones(1,size(fmaxc,1));
fmaxcmat = ones(size(mgcmat,1),1)*fmaxc';
mgcd = ((mgcmat - fmaxcmat)').^2;

dp = sqrt(mgrd+mgcd);

%%
% figure;
% set(gcf,'position',[         1000         163        1348        1153]);

% set(gcf,'position',[   1566         968         653         335]);
% ha = tight_subplot(8,4,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); imagesc(rand(10,10)); axis([1 5 1 5]); axis off; end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% moff = [0 26*26 26*26+30*30+60*60 26*26+30*30 51*51];


%  27*31+31*35+54*62+63*72
onParasolInd = 1+27*31;
offParasolInd = onParasolInd+31*35;
onMidgetInd = offParasolInd+54*62;
offMidgetInd = onMidgetInd+63*72;
moff = [0 onParasolInd offParasolInd onMidgetInd offMidgetInd];

% spikeNoLearn(onParasolInd+1:offParasolInd,:) = 0;
% spikeNoLearn(onMidgetInd+1:end,:) = 0;
% spikeNoLearn(1:onParasolInd,:) = 0;
% spikeNoLearn(offParasolInd+1:onMidgetInd,:) = 0;


cellArr = [122 322 874 1134];

pltCtr = 1;

titleStr{1} = 'on parasol';
titleStr{2} = 'off parasol';
titleStr{3} = 'off midget';
titleStr{4} = 'on midget';

windsize=20;
signValRF =[1 1 1 1];
signValWN = [1 1 1 1];
signVal = [1 1 1 1];
mindCtr = 1;
% 
% ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);

for mosind =4%[1 2 3 4]
    figure;
    
    ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);

    nCells = 4;
    
signValRF =ones(1,nCells);
signValWN =ones(1,nCells);
signVal = ones(1,nCells);
%     cellArr = round(3000*rand(nCells,1));
% cellArr = round((moff(mosind+1)-moff(mosind))*rand(nCells,1));

if mosind == 1
    cellArr = [229 156+137 371 314];
elseif mosind == 2
    cellArr = [293 609 283+20 845];
    
elseif mosind == 3
    cellArr = [247+499 1177 1540 1267];
    
elseif mosind == 4
    cellArr = [2310 3326 1514+100 -100+2639];
% else
%     
%     cellArr = round((moff(mosind+1)-moff(mosind))*rand(nCells,1));
end

% cellArr(:,4) = [
%         1434
%         1438
%         1726
%         2831];
% ha = tight_subplot(nCells,4,[.01 .03],[.1 .01],[.01 .01]);
hold on;
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%     imc = zeros(nCells*(windSize+1),4*(windSize+1));
    for mind = 1:nCells
    
    paraInd = cellArr(mind);
    paraIndPlus = cellArr(mind)+moff(mosind)
    mindCtr = mindCtr+1;
        
    
     
    
%     tight_subplot(4,4,pltCtr); 
    axes(ha(pltCtr));
    pltCtr = pltCtr+1;
    filtIm = reshape(filterNS(1+paraIndPlus,:),[100 100]);
    imagesc(signVal(mind)*filtIm);
    mabs = max(abs(filtIm(:)));
    minabs = min((filtIm(:)));
    caxis([1*minabs 1*mabs]); colormap gray;
    axis image;
%     title(sprintf([titleStr{mind} '\nCell %d NS Decoding'], paraInd),'fontsize',14); % set(gca,'fontsize',14);
%     if mod(mindCtr,2)
    axis([fmaxr(paraIndPlus+1)-windsize , fmaxr(paraIndPlus+1)+windsize,fmaxc(paraIndPlus+1)-windsize , fmaxc(paraIndPlus+1)+windsize]);
%     end
    
    axis off
%     
%     
% %     tight_subplot(4,4,pltCtr); 
%     axes(ha(pltCtr));
%     pltCtr = pltCtr+1;
%     filtIm = reshape(filterPros(1+paraIndPlus,:),[100 100])';
%     imagesc(signVal(mind)*filtIm);
%     mabs = max(abs(filtIm(:)));
%     minabs = min((filtIm(:)));
% %     caxis([minabs 1*mabs]); 
% %     if mosind == 1 || mosind == 3
% %     caxis([-mabs -.5*mabs]);
% %     else
% %     caxis([.5*mabs mabs]);
% %     end
% colormap gray;
%     axis image;
% %     title(sprintf([titleStr{mind} '\nCell %d NS Decoding'], paraInd),'fontsize',14); % set(gca,'fontsize',14);
% %     if mod(mindCtr,2)
%     axis([fmaxr(paraIndPlus+1)-windsize , fmaxr(paraIndPlus+1)+windsize,fmaxc(paraIndPlus+1)-windsize , fmaxc(paraIndPlus+1)+windsize]);
% %     end
%     
%     axis off
    
%     subplot(4,4,pltCtr); pltCtr = pltCtr+1;
%     %     irf = innerRetina.mosaic{mind}.tCenter{1};
%     %     plot(.008*(1:length(irf)),irf); axis square
%     plot([1:size(tFiltInt,1)]/121,tFiltInt(:,mind)); hold on;
%     plot(tAxis,tFilt(:,mind)+.6*(-.5+rand(size(tFilt,1),1)),'r');
%     
%     %     plot([1:size(tFiltInt,1)]/121,tFiltInt(:,mind)+.02*(-.5+rand(size(tFiltInt,1),1)),'-r'); hold on;
%     axis([0 0.4 -.5 1.1]);
%     xlabel('Time (sec)');
%     ylabel('AU');
%     grid on;
%     if mind ==1
%         legend('WN/NS STA','Decoding');delay:end
%     end;
%     set(gca,'fontsize',12);
    end
end
% savefig(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/healthy' num2str(mosind) '.fig']);
% savefig(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/filters/healthy' num2str(mosind) '.fig']);
savefig(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/filters/july14_healthy' num2str(mosind) '.fig']);
% print(gcf,'-deps',['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/filters/healthy' num2str(mosind) '.eps']);

% print(gcf,'-dsvg',['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/filters/healthy' num2str(mosind) '.svg'])
%%
pltCtr = 1;
mindCtr = 1;
for mosind =4%[1 2 3 4]
    figure;
    
    ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);

    nCells = 4;
    
signValRF =ones(1,nCells);
signValWN =ones(1,nCells);
signVal = ones(1,nCells);
%     cellArr = round(3000*rand(nCells,1));
% cellArr = round((moff(mosind+1)-moff(mosind))*rand(nCells,1));

if mosind == 1
    cellArr = [229 156+137 371 314];
    
elseif mosind == 2
    cellArr = [293 609 283+20 845];
    
elseif mosind == 3
    cellArr = [247+499 1177 1540 1267];
    
    
elseif mosind == 4
    cellArr = [2310 3326 1514+100 -100+2639];
    
else
    
    cellArr = round((moff(mosind+1)-moff(mosind))*rand(nCells,1));
end

% cellArr(:,4) = [
%         1434
%         1438
%         1726
%         2831];
% ha = tight_subplot(nCells,4,[.01 .03],[.1 .01],[.01 .01]);
hold on;
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%     imc = zeros(nCells*(windSize+1),4*(windSize+1));
    for mind = 1:nCells
    
    paraInd = cellArr(mind);
    paraIndPlus = cellArr(mind)+moff(mosind)
    mindCtr = mindCtr+1;
        
    
     
%     
% %     tight_subplot(4,4,pltCtr); 
%     axes(ha(pltCtr));
%     pltCtr = pltCtr+1;
%     filtIm = reshape(filterNS(1+paraIndPlus,:),[100 100]);
%     imagesc(signVal(mind)*filtIm);
%     mabs = max(abs(filtIm(:)));
%     minabs = min((filtIm(:)));
%     caxis([minabs 1*mabs]); colormap gray;
%     axis image;
% %     title(sprintf([titleStr{mind} '\nCell %d NS Decoding'], paraInd),'fontsize',14); % set(gca,'fontsize',14);
% %     if mod(mindCtr,2)
%     axis([fmaxr(paraIndPlus+1)-windsize , fmaxr(paraIndPlus+1)+windsize,fmaxc(paraIndPlus+1)-windsize , fmaxc(paraIndPlus+1)+windsize]);
% %     end
%     
    axis off
%     
    
%     tight_subplot(4,4,pltCtr); 
    axes(ha(pltCtr));
    pltCtr = pltCtr+1;
    filtIm = reshape(filterPros(1+paraIndPlus,:),[100 100])';
    imagesc(signVal(mind)*filtIm);
    mabs = max(abs(filtIm(:)));
    minabs = min((filtIm(:)));
%     caxis([minabs 1*mabs]); 
%     if mosind == 1 || mosind == 3
%     caxis([-mabs -.5*mabs]);
%     else
%     caxis([.5*mabs mabs]);
%     end
colormap gray;
    axis image;
%     title(sprintf([titleStr{mind} '\nCell %d NS Decoding'], paraInd),'fontsize',14); % set(gca,'fontsize',14);
%     if mod(mindCtr,2)
    axis([fmaxr(paraIndPlus+1)-windsize , fmaxr(paraIndPlus+1)+windsize,fmaxc(paraIndPlus+1)-windsize , fmaxc(paraIndPlus+1)+windsize]);
%     end
    
    axis off
    
%     subplot(4,4,pltCtr); pltCtr = pltCtr+1;
%     %     irf = innerRetina.mosaic{mind}.tCenter{1};
%     %     plot(.008*(1:length(irf)),irf); axis square
%     plot([1:size(tFiltInt,1)]/121,tFiltInt(:,mind)); hold on;
%     plot(tAxis,tFilt(:,mind)+.6*(-.5+rand(size(tFilt,1),1)),'r');
%     
%     %     plot([1:size(tFiltInt,1)]/121,tFiltInt(:,mind)+.02*(-.5+rand(size(tFiltInt,1),1)),'-r'); hold on;
%     axis([0 0.4 -.5 1.1]);
%     xlabel('Time (sec)');
%     ylabel('AU');
%     grid on;
%     if mind ==1
%         legend('WN/NS STA','Decoding');delay:end
%     end;
%     set(gca,'fontsize',12);
    end
end
savefig(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/filters/july14_pros_' num2str(mosind) '.fig']);
%%

% sta = stimzm(:,9:end)*single(spikeResp(:,1:end-8)');
% figure; staex = sta(:,323); staim = reshape(staex,[100 100]); stamax = max(abs(staex(:))); imagesc(staim); caxis([-stamax stamax]);
% figure; staex = sta(:,6323); staim = reshape(staex,[100 100]); stamax = max(abs(staex(:))); imagesc(staim); caxis([-stamax stamax]);

% stim = RGB2XWFormat(testmovieshort);
% delay = 6;
% stimwnzm = (single(stim)-(ones(size(stim,1),1)*mean(stim,1)));
% stawn = stimwnzm(:,1:end-(delay-1))*single(spikeResp(:,delay:end)');

% figure; staex = stawn(:,6323); staim = reshape(staex,[100 100]); stamax = max(abs(staex(:))); imagesc(staim); caxis([-stamax stamax]);
