% f_8
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 6 reconstruction filters for healthy vs. prosthesis
%
% Compares the reconstruction filters for the healthy and prosthesis
% repsonses.
%
% 2018 JRG (c) isetbio team
% [formerlhy fitlersAndSTAs_paper_august]

%% Load data

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0.mat');
% filterNS = filterMat; clear filterMat
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv 5_w1_sh4_dr0_pitch_70_decay_2.mat')
% filterPros = filterMat;

rd = RdtClient('isetbio');
rd.crp('/resources/data/reconstruction');
filterFile = 'filtersmosaic0_sv50_w1_sh15_dr0_aug27.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterNS = data.filterMat; clear data;

rd = RdtClient('isetbio');
rd.crp('/resources/data/reconstruction');
filterFile = 'filtersmosaic0_sv05_w1_sh4_dr0_pitch_70_decay_2_aug29.mat';
% filterFile = 'filtersmosaic0_sv5_w1_sh4_dr0_pitch_70_decay_2_aug27.mat'
% filterFile = 'filtersmosaic0_sv10_w1_sh4_dr0_pitch_70_decay_2.mat';
% rd.crp('/resources/data/istim');
% filterFile = 'filters_mosaic0_sv10_w1_sh2_dr0.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterPros = data.filterMat; clear data;
% filterPros = filterNS;

%% Get max location of each filter for zooming
[mgr,mgc] = meshgrid(1:100,1:100);

% Find row and col location of max value
[cmax,cind] = max(abs(filterNS),[],2);
[fmaxc,fmaxr] = ind2sub([100 100],cind);

mgrmat = mgr(:)*ones(1,size(fmaxr,1));
fmaxrmat = ones(size(mgrmat,1),1)*fmaxr';
mgrd = ((mgrmat - fmaxrmat)').^2;

mgcmat = mgc(:)*ones(1,size(fmaxc,1));
fmaxcmat = ones(size(mgcmat,1),1)*fmaxc';
mgcd = ((mgcmat - fmaxcmat)').^2;

dp = sqrt(mgrd+mgcd);

%% Choose cells to look at

% Total number of cells in each mosaic
% 28*32+31*35+1:28*32+31*35+55*63+1

onParasolInd = 1+28*32;
offParasolInd = onParasolInd+31*35;
onMidgetInd = offParasolInd+55*63;
offMidgetInd = onMidgetInd+4371;
moff = [0 onParasolInd offParasolInd onMidgetInd offMidgetInd];

cellArr = [122 322 874 1134];

titleStr{1} = 'on parasol';
titleStr{2} = 'off parasol';
titleStr{3} = 'off midget';
titleStr{4} = 'on midget';

windsize=24;
signValRF =[1 1 1 1];
signValWN = [1 1 1 1];
signVal = [1 1 1 1];
mindCtr = 1;

for mosind =1:2%[1 2 3 4]
    figure;
    set(gcf,'position',[440   367   385   431]);
    ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);
    
    nCells = 4;
    
    signValRF =ones(1,nCells);
    signValWN =ones(1,nCells);
    signVal = ones(1,nCells);
    
    if mosind == 1
        cellArr = [229+13 156+137+59 371 229+13];
    elseif mosind == 2
        cellArr = [293 609+33 283+41 845];
    elseif mosind == 3
        cellArr = [247+499 1177 1540 1267];
        
    elseif mosind == 4
        cellArr = [2310 3326 1514+100 -100+2639];
    end
    
    hold on;
    
    pltCtr = 1;
    
    for mind = 1:nCells
        
        paraInd = cellArr(mind);
        paraIndPlus = cellArr(mind)+moff(mosind);
        mindCtr = mindCtr+1;
        
        axes(ha(pltCtr));
        pltCtr = pltCtr+1;
        filtIm = reshape(filterNS(1+paraIndPlus,:),[100 100]);
        imagesc(signVal(mind)*filtIm);
        mabs = max(abs(filtIm(:)));
        minabs = min((filtIm(:)));
        
        colormap gray;
        axis image;
        
        axis([fmaxr(paraIndPlus+1)-windsize , fmaxr(paraIndPlus+1)+windsize,fmaxc(paraIndPlus+1)-windsize , fmaxc(paraIndPlus+1)+windsize]);
        axis off
    end
end

set(gcf,'PaperPositionMode','auto')
% print(gcf,'-dsvg',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/healthy' num2str(mosind) '.svg'])

%% Prosthesis filters

mindCtr = 1;
for mosind =1:2%[1 2 3 4]
    figure;
    
    set(gcf,'position',[440   367   385   431]);
    ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);
    
    nCells = 4;
    
    signValRF =ones(1,nCells);
    signValWN =ones(1,nCells);
    signVal = ones(1,nCells);
    
    if mosind == 1
        cellArr = [229+13 156+137+59 371 229+13];
        
    elseif mosind == 2
        cellArr = [293 609+33 283+41 845];
        
    elseif mosind == 3
        cellArr = [247+1499 1177 1540 1267];
        
        
    elseif mosind == 4
        cellArr = [2310 3326 1514+100 -100+2639];
        
    else
        
        cellArr = round((moff(mosind+1)-moff(mosind))*rand(nCells,1));
    end
    
    hold on;
    pltCtr = 1;
    for mind = 1:nCells
        
        paraInd = cellArr(mind);
        paraIndPlus = cellArr(mind)+moff(mosind);
        mindCtr = mindCtr+1;
        
        axis off
        axes(ha(pltCtr));
        pltCtr = pltCtr+1;
        filtIm = reshape(filterPros(1+paraIndPlus,:),[100 100]);
        imagesc(signVal(mind)*filtIm);
        mabs = max(abs(filtIm(:)));
        minabs = min((filtIm(:)));
        colormap gray;
        axis image;
        axis([fmaxr(paraIndPlus+1)-windsize , fmaxr(paraIndPlus+1)+windsize,fmaxc(paraIndPlus+1)-windsize , fmaxc(paraIndPlus+1)+windsize]);
        axis off
        
    end
end
% savefig(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/filters/july14_pros_' num2str(mosind) '.fig']);

set(gcf,'PaperPositionMode','auto')
% print(gcf,'-dsvg',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/pros' num2str(mosind) '.svg'])
%%
%%

% sta = stimzm(:,9:end)*single(spikeResp(:,1:end-8)');
% figure; staex = sta(:,323); staim = reshape(staex,[100 100]); stamax = max(abs(staex(:))); imagesc(staim); caxis([-stamax stamax]);
% figure; staex = sta(:,6323); staim = reshape(staex,[100 100]); stamax = max(abs(staex(:))); imagesc(staim); caxis([-stamax stamax]);

% stim = RGB2XWFormat(testmovieshort);
% delay = 6;
% stimwnzm = (single(stim)-(ones(size(stim,1),1)*mean(stim,1)));
% stawn = stimwnzm(:,1:end-(delay-1))*single(spikeResp(:,delay:end)');

% figure; staex = stawn(:,6323); staim = reshape(staex,[100 100]); stamax = max(abs(staex(:))); imagesc(staim); caxis([-stamax stamax]);
