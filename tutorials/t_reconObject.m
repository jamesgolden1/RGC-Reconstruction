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
% INSTRUCTIONS:
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
clear
% 

% folderName = 'aug30';

% folderName = 'sep20hall100'; pRecon.testFlag = 1;
% folderName = 'sep20prima35test'; pRecon.testFlag = 1;
% folderName = 'sep20prima35hall'; pRecon.testFlag = 1;
% folderName = 'sep20prima35test'; pRecon.numTest = 20;
% folderName = 'nov_results/healthyTest2/'; pRecon.testFlag = 0;
folderName = 'nov_results/primaTestLearning2/'; pRecon.testFlag = 0;
% folderName = 'aug8';
mosaicFile = 'mosaic0';

movieFile  = fullfile(folderName, 'mov');
spikesFile = fullfile(folderName, 'sp');
buildFile  = fullfile(folderName, 'raw','build');

% movieFile  = fullfile(reconstructionRootPath, 'dat', folderName, 'mov');
% spikesFile = fullfile(reconstructionRootPath, 'dat', folderName, 'sp');
% buildFile  = fullfile(reconstructionRootPath, 'dat', folderName, 'build');

windowSize = 1;
percentSV = .25;%.5;%.05;%.5;%.5;
% shifttime = 2;
shifttime = 3;%15;%4;%15;%17;%5;
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
pRecon.stimTypeBuild = 'ns500';
% pRecon.stimTypeBuild = 'wn'; 

reconHealthy = recon(pRecon);
 
blockIn = 1; startInd = 1;
% nCores = 18;
% pool = parpool(nCores);
% for startInd = 5:20
% for ii = 27:1:32%28:2:32%1:26%27%:floor(576/nCores)%:floor(576/nCores)]
% parfor blockIn = [nCores*(ii-1)+1:nCores*(ii)]
%     reconHealthy.build(pRecon,'blockIn',blockIn,'startInd',startInd);
    
    pRecon.pixelWidth = 70/1;
    reconHealthy.buildPrima(pRecon,'blockIn',blockIn,'startInd',startInd);

%     reconHealthy.buildHallway(pRecon);
%     reconHealthy.buildPrimaHallway(pRecon);
% end
% delete(pool);
% end
% end
% delete(gcp);
reconHealthy.loadSpikes(pRecon);
% reconHealthy.loadSpikes
% reconHealthy.train(pRecon,'shifttime',shifttime);
% c
% reconHealthy.test(pRecon);
% reconHealthy.movie();



movieFile  = fullfile(folderName);
spikesFile = fullfile(folderName);
buildFile  = fullfile(folderName, 'raw','build');

pRecon.buildFile = buildFile;
pRecon.stimFile = movieFile;
pRecon.respFile = spikesFile;

pRecon.spatialFilterLambda = 0.001;

    pRecon.testShift = 0;
    
    reconHealthy = recon(pRecon);
    
    mse1 = reconHealthy.testImagenet(pRecon);