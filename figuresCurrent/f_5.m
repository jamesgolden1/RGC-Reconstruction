% f_5
% 
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 4 prosthesis RGC reconstruction of leaf image
%
% Generates the different types of reconstructed images for the leaf image
% stimulus with prosthesis stimulation as well as activaiton maps for the
% RGC types.
%
% 2017 JRG (c) isetbio team
% 

%%
pRecon.degenFlag = 1;
reconHealthy = recon(pRecon);
reconHealthy.buildPrimaLeaf();%'pixelWidth',35);
% reconHealthy.buildPrimaWheel('pixelWidth',70);