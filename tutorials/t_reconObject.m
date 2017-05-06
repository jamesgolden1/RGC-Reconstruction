% t_reconObject
% 
% The recon object is used to reconstruct a movie stimulus from RGC spikes.
% 
% In order to use the recon object, one must be running isetbio from
% vision@bertha.stanford.edu in order to access the imageNet stimulus
% blocks.
% 
% The recon object is specified primarily by file names. Unfortunately, for
% a large RGC mosaic (>1000 cells), the variables becomes so large that is
% easier to store them on disk and call them when necessary.
% 
% The recon object is built by specifying a number of parameters. The
% object is then created. The 'build' method calls on isetbio to simulate
% the retina through to RGCs in response to a large database of imageNet
% patches.
% 
% The 'train' method takes the RGC spike responses and stimulus images and
% uses them to learn decoding filters.
% 
% The 'plot' method visualizes the decoding filters, while the 'test'
% method runs the reconstruction for a subset of data and evaluates the
% accuracy of the reconstructions.
% 
% JRG NP (Nikhil Parthsarathy) 4/2017

%%

% New mosaic
mosaicFile = '_test4_mosaic_may5';%_mosaicAll_23282';
movieFile = 'test4/mov_may5'; 
spikesFile = 'test4/sp_may5';


% Old mosaic
% mosaicFile = '_mosaicAll_1246640' ;
% movieFile = 'ns100_r2_10/ns100_jan1_mov3'; 
% spikesFile = 'ns100_r2_10/ns100_jan1_sp3';

% Set filter file name
% evArr = [.2 .1 .3 .4];
% trainFraction = [.16 .5 .66 .83 1];
% evInd = 4; trainFractionInd = 5;
% filterFile = ['test4/filters_ns100_may1_sh9_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];

filterFile = ['test4/filters_may5'  mosaicFile];

clear pRecon
pRecon.buildFile = 'test4/build_may5';
pRecon.stimFile = movieFile;
pRecon.respFile = spikesFile;
pRecon.filterFile = filterFile;
pRecon.mosaicFile = mosaicFile;
pRecon.windowSize = 1;
pRecon.percentSV = .5;
% pRecon.stimType = 'wn';

reconHealthy = recon(pRecon);

blockIn = 1;

% pool = parpool(16);
% parfor blockIn = 1:16
%     reconHealthy.build(pRecon,'blockIn',blockIn);
% end
% delete(pool);

reconHealthy.train(pRecon,'shifttime',15);
reconHealthy.plot('filters');
% reconHealthy.test(pRecon);
