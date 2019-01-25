% t_reconObject
% 
% The recon object is used to reconstruct a movie stimulus from RGC spikes.
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
% In order to train the reconstruction filter, one must be running isetbio
% from the appropriate Stanford server in order to access the imageNet
% stimulus blocks.
% 
% The 'train' method takes the RGC spike responses and stimulus images and
% uses them to learn decoding filters.
% 
% The 'plot' method visualizes the decoding filters, while the 'test'
% method runs the reconstruction for a subset of data and evaluates the
% accuracy of the reconstructions.
% 
% 
% INSTRUCTIONS: see line 75 if parpool is available to cut down on training
% time. The default is a standard for loop, 
% 
% TOOLBOX DEPENDENCIES - these must be downloaded and added to the matlab
%                               path with subfolders.
%       isetbio:            http://github.com/isetbio/isetbio [branch:phosphenePaper]
%       RGC-Reconstruction: https://github.com/Chichilnisky-Lab/RGC-Reconstruction
%       EJLPhosphene:       https://github.com/isetbio/EJLPhosphene
%       RemoteDataToolbox:  https://github.com/isetbio/RemoteDataToolbox
% 
% (c) isetbio team 2017 JRG  NP (Nikhil Parthsarathy)
%%

isetbioLocalPath = '/Volumes/Lab/Users/james/';
addpath(genpath([isetbioLocalPath 'isetbio']));

rdtLocalPath = '/Volumes/Lab/Users/james/';
addpath(genpath([rdtLocalPath 'RemoteDataToolbox']));

phospheneLocalPath = '/Volumes/Lab/Users/james/current';
addpath(genpath([phospheneLocalPath 'EJLPhosphene']));

reconLocalPath = '/Volumes/Lab/Users/james/current';
addpath(genpath([reconLocalPath 'RGC-Reconstruction/']));

%%

% Each image block is 500 different images
% Total available is 360*500 images
% Save 20% for test set
numberImageBlocks = 360;

folderName = 'healthy'; 
mosaicFile = 'mosaic0';
pRecon.testFlag = 0;

movieFile  = fullfile(folderName, 'mov');
spikesFile = fullfile(folderName, 'sp');
buildFile  = fullfile(folderName, 'raw','build');

windowSize = 1;
percentSV = .5;
shifttime = 15;
dropout = 0;

filterFile  = fullfile(folderName,...    
    ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);

pRecon.buildFile = buildFile;
pRecon.stimFile = movieFile;
pRecon.respFile = spikesFile;
pRecon.filterFile = filterFile;
pRecon.mosaicFile = mosaicFile;
pRecon.windowSize = windowSize;
pRecon.percentSV = percentSV;
pRecon.dropout = dropout;
% pRecon.stimTypeBuild = 'ns500';
pRecon.stimTypeBuild = 'ns';
% pRecon.stimTypeBuild = 'wn'; 

reconHealthy = recon(pRecon);

% % For parallel pool, uncomment below:
% blockIn = 1; 
% nCores = 18;
% pool = parpool(nCores);
% for ii = 1:floor(numberImageBlocks/nCores)
% parfor blockIn = [nCores*(ii-1)+1:nCores*(ii)]
for blockIn = 1:numberImageBlocks
    % Healthy training
    reconHealthy.build(pRecon,'blockIn',blockIn);
    
    % Prosthesis training
%     pRecon.pixelWidth = 70/1;
%     reconHealthy.buildPrima(pRecon,'blockIn',blockIn);
end
% delete(pool);
% end
% end
% delete(gcp);
reconHealthy.loadSpikes(pRecon);
reconHealthy.train(pRecon,'shifttime',shifttime);
