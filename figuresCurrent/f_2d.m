% f_2d
% 
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 2d RGC mosaic with prosthesis array
%
% Generates an RGC mosaic with four types of cells and subretinal
% prosthesis model and shows the prosthesis array overlaid on the RGC
% mosaic of different types.
%
% 2017 JRG (c) isetbio team
% 
% [formerly pixel_array_mosaic.m]

%% Load mosaic
reconHealthy = recon();

[reconHealthy,primaRecon] = reconHealthy.buildPrimaMosaic();

% Or load from .mat file
% load([reconstructionRootPath '/dat/figures_data/primaRecon.mat'],'primaRecon');

innerRetina = primaRecon.innerRetina;

%% Plot mosaic and prothesis array
for typeInd = 1%:4;
    innerRetina.mosaic{typeInd}.plot('mosaic')
    hold on;
    primaRecon.plot('array')
    axis off
end
