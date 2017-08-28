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
clear

pRecon.pixelWidth = 70/1;
folderName = 'aug27test';%prima70';
% folderName = 'aug122prima18test';

% folderName = 'aug8';
mosaicFile = 'mosaic0';

movieFile  = fullfile(folderName, 'mov');
spikesFile = fullfile(folderName, 'sp');
buildFile  = fullfile(folderName, 'raw','build');

% movieFile  = fullfile(reconstructionRootPath, 'dat', folderName, 'mov');
% spikesFile = fullfile(reconstructionRootPath, 'dat', folderName, 'sp');
% buildFile  = fullfile(reconstructionRootPath, 'dat', folderName, 'build');

windowSize = 1;
percentSV = .5;%.5;
% shifttime = 2;
shifttime =15;%4;%15;%17;%5;
dropout = 0;

filterFile  = fullfile(folderName,...    
    ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);
% filterFile  = fullfile(reconstructionRootPath, 'dat', folderName,...    
%     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);
%     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime)]);

pRecon.buildFile = buildFile;
pRecon.stimFile = movieFile;
pRecon.respFile = spikesFile;
pRecon.filterFile = filterFile;
pRecon.mosaicFile = mosaicFile;
pRecon.windowSize = windowSize;
pRecon.percentSV = percentSV;
pRecon.dropout = dropout;
pRecon.stimTypeBuild = 'ns';
% pRecon.stimTypeBuild = 'wn'; 

reconHealthy = recon(pRecon);
 
% blockIn = 200;
% nCores = 12;
% pool = parpool(nCores);
% for ii = 1:39%27%:floor(576/nCores)%:floor(576/nCores)]
% parfor blockIn = [nCores*(ii-1)+1:nCores*(ii)]
%     reconHealthy.build(pRecon,'blockIn',blockIn);
%     reconHealthy.buildPrima(pRecon,'blockIn',blockIn);
    reconHealthy.buildHallway(pRecon);
%     reconHealthy.buildPrimaHallway(pRecon);
% end
% delete(pool);
% end
% delete(gcp);
% reconHealthy.train(pRecon,'shifttime',shifttime);
% reconHealthy.plot('filters');
% reconHealthy.test(pRecon);
% reconHealthy.movie();
