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
clear


% addpath(genpath('/Volumes/Lab/Users/james/EJLPhosphene'))
% addpath(genpath('/Volumes/Lab/Users/james/current/RGC-Reconstruction/'))
% addpath(genpath('/Volumes/Lab/Users/james/current/isetbio'))

% Each image block is 500 different images
% Total available is 576*500 images
% Save 20% for test set
numberImageBlocks = 576;%460;

% folderName = 'healthy_training'; 
drctr = 0;
mse35 =[];
cc35=[];

% pixelWidth = 70/2;

%  for dropout =  .99:-.01:.91%.9:-.1:.1
    

% folderName = 'aug29prima35';
folderName = 'prosthesis_35_training_aug13'; 
% folderName = 'prosthesis_70_real_testing_aug6'; 
% folderName = 'prosthesis_70_testing_july16'; 
% folderName = 'prosthesis_35_training_aug7'; 
% folderName = 'aug27';
mosaicFile = 'mosaic0';
% pRecon.testFlag = 0;
% 
movieFile  = fullfile(folderName, 'mov');
spikesFile = fullfile(folderName, 'sp');
buildFile  = fullfile(folderName, 'raw','build');

% % healthy
% windowSize = 1;
% percentSV = .5;
% shifttime = 15;
% dropout = 0;
% 
% filterFile  = fullfile(folderName,...    
%     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);% '_pitch' sprintf('_%d',pixelWidth) '_decay_2']);

pixelWidth = 70/2;
% prosthesis
windowSize = 1;
percentSV = .05;
shifttime =3;
dropout = 0;%.3;


% pRecon.pixelWidth = pixelWidth;

% filterFile  = fullfile(folderName,...    
%     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%.0f',100*dropout) '_pitch' sprintf('_%d',pixelWidth) '_decay_2']);
filterFile  = fullfile(folderName,...    
    ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%.0f',100*0) '_pitch' sprintf('_%d',pixelWidth) '_decay_2']);

% filterFile  = fullfile(folderName,...    
%     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);


for dropout = [0 ]% [.9:-.1:.1]%.9:-.1:.7 .5 .3:-.1:.1]%.9:-.1:.1
    
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


pRecon.testFlag = 1;

reconHealthy = recon(pRecon);
%%

% reconHealthy.loadSpikes(pRecon);
% reconHealthy.train(pRecon,'shifttime',shifttime);
%%
% % % For parallel pool, uncomment below:
% % % blockIn = 1; 
% nCores = 12;
% % pool = parpool(nCores);
% % pRecon.pixelWidth = 70/2;
% % for ii = 27:floor(numberImageBlocks/nCores) % 18 core
% for ii = 40:floor(numberImageBlocks/nCores) % 12 cores
%     
% for startInd =1% 6:20
% for blockIn = [nCores*(ii-1)+1:nCores*(ii)]
% % for blockIn = 1:numberImageBlocks
%     % Healthy training
% %     reconHealthy.build(pRecon,'blockIn',blockIn);
%     
%     % Prosthesis training
%     reconHealthy.buildPrima(pRecon,'blockIn',blockIn,'startInd',startInd);
% % end
% end
% % 
% delete(gcp);
% end
% end
% delete(gcp);
% % reconHealthy.loadSpikes(pRecon);
% % % reconHealthy.train(pRecon,'shifttime',shifttime);
% 
%%
% % folderName = 'testHealthy'; 
% % folderName = 'testPros70_july14';
% folderName = 'prosthesis_35_testing_july16';
% folderName = 'prosthesis_testing_70_aug22_noLearn'; pRecon.pixelWidth = 70/1;

folderName = 'prosthesis_35_testing_aug13'; pRecon.pixelWidth = 70/2;
% folderName = 'aug29test';
% folderName = 'healthy_testing_aug12'
% folderName = 'healthy_testing_aug12'
% Stimulus is 25 natural scene images presented for 20 frames each at 0.001 sec
% pRecon.stimTypeBuild = 'ns500'; 

% Set destination files for movie and spike data
pRecon.stimFile  = fullfile(folderName,'mov');
pRecon.respFile  = fullfile(folderName,'sp');
pRecon.buildFile = fullfile(folderName, 'raw','build');

pRecon.testFlag = 1;
pRecon.onlyOnFlag = 0;
pRecon.noLearnFlag = 0;
% pRecon.stimTypeBuild = 'ns500';
% reconHealthy.build(pRecon);

pRecon.dropoutFlag = 1;

% reconHealthy.loadSpikes(pRecon);

pRecon.numTest = 1;
reconHealthy = recon(pRecon);
% reconHealthy.loadSpikes(pRecon);
    
pRecon.stimFile  = fullfile(folderName);%,'mov');
pRecon.respFile  = fullfile(folderName);%,'sp');

drctr = drctr+1;
tsctr = 1;    
for testShift = 1%1%:5%:20
    
%     pRecon.testShift = testShift;
%     pRecon.startInd = startInd;

    pRecon.testShift = testShift;
    
    % reconHealthy.train(pRecon,'shifttime',shifttime);

%     testShift=1;
    pRecon.spatialFilterLambda =.001; %pros
    
%     pRecon.spatialFilterLambda = .06; % healthy
    [mse1, cc] = reconHealthy.testImagenet(pRecon); 
    cc35(drctr,tsctr) = cc;
    mse35(drctr,tsctr) = mse1;
    tsctr = tsctr+1;
end
% disp('Healthy: 0.0952    0.8305');   

end