% f_2a
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 2a RGC spatial RF
%
% Generates an RGC mosaic with four types of cells and subretinal
% prosthesis model and shows 2D spatial RF profiles.
%
% 2017 JRG (c) isetbio team

%% Generate mosaic and prosthesis array

reconHealthy = recon();

[reconHealthy,primaRecon] = reconHealthy.buildPrimaMosaic();

% Or load from .mat file
% load([reconstructionRootPath '/dat/figures_data/primaRecon.mat'],'primaRecon');

innerRetina = primaRecon.innerRetina;

%% Generate plot of
figure;
% set(gcf,'position',[ 440   359   977   439]);

titlestr{1} = 'isetbio on parasol RGC RF';
titlestr{2} = 'isetbio off parasol RGC RF';
titlestr{3} = 'isetbio on midget RGC RF';
titlestr{4} = 'isetbio off midget RGC RF';

bipolarMosaicSizeMicrons = 1e6*primaRecon.bpMosaic.mosaic{1}.patchSize(1);
bipolarNumber = size(primaRecon.bpMosaic.mosaic{1}.cellLocation,1);
bipolarCellSizeMicrons = bipolarMosaicSizeMicrons/bipolarNumber; % = 2

for mosaicNumber =1:4;
    clear s2
    
    s2=(innerRetina.mosaic{mosaicNumber}.sRFcenter{1,1}-...
        innerRetina.mosaic{mosaicNumber}.sRFsurround{1,1}); 
    
    subplot(2,2,mosaicNumber);
    
    [mv,mi1]=max(max(s2));
    [mv,mi2]=max(max(s2'));
    
    plot(-size(s2,1)-1+bipolarCellSizeMicrons*(1:length(s2(mi1,:))),s2(mi1,:));
    
    axis([-50 50 -.02 1.05*max(s2(:))]);
    xlabel(sprintf('\\mum'));
    set(gca,'fontsize',14);
    grid on;
    
    title(titlestr{mosaicNumber});
end
