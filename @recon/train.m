function obj = train(obj, varargin)
%TRAIN - Run the reconstruction algorithm to generate the filters

% Must run loadSpikesAll.m to build stimFile and respFile

p = inputParser;
p.addParameter('stimFile',[],@ischar);
p.addParameter('respFile',[],@ischar);
p.addParameter('filterFile',[],@ischar);
p.addParameter('windowSize',8,@isnumeric);
p.addParameter('percentSV',0.5,@isnumeric);
p.addParameter('mosaicFile',[],@ischar);
p.addParameter('trainFraction',1,@isnumeric);
p.addParameter('shiftTime',0,@isnumeric);
p.addParameter('stimType','ns',@ischar);
p.parse(varargin{:});
filterFile = p.Results.filterFile;
stimFileName = p.Results.stimFile;
respFileName = p.Results.respFile;
windowSize = p.Results.windowSize;
percentSV = p.Results.percentSV;
mosaicFile = p.Results.mosaicFile;
trainSizeArray = p.Results.trainFraction;
shiftTime = p.Results.shiftTime;
stimType = p.Results.stimType;

if isempty(filterFile)
    filterFile = ['filters_' num2str(round(cputime*100))];
end

%% run Linear reconstruction for on, off, and joint on/off training 

disp('Loading spike responses...')

% stimFileName = 'NSmovie_spikeResp_overlap0';
% respFileName = 'NSspikeResp_overlap0';

% stimFileName = 'NSmovie_40reps_ns0';
% respFileName = 'NSspikeResp_40reps_ns0';

% matfON = matfile([reconstructionRootPath '\dat\' respFileName]);

% matfON = matfile([reconstructionRootPath '/dat/' respFileName mosaicFile]);

matfON = matfile('/Volumes/Lab/Users/james/RGC-Reconstruction/dat/ns100_r2_10/ns100_jan1_sp3_mosaicAll_1246640');

% movielength = 1*240000;%size(stim,2);
% disp(['Total Movie Length in Frames: ' num2str(movielength)]);

fileext = 'mosaic_ns_all_mult';
% trainSizeArray = 1;%[.6/8:.6/8:.6]; 
trainInd = 1;
% includedComponentsArray = 1000;
includedComponentsArray = percentSV;

srON = matfON.spikeResp;
% srOFF = matfOFF.spikeResp;
% clear matfON matfOFF
% spikeResp1 = vertcat(srON(:,1:12*12000));%, srOFF(:,1:240000), matfONP.spikeResp, matfOFFP.spikeResp);
spikeResp1 = srON;

% for ri = 1:2:size(spikeResp1,1)
%     spikeResp1(ri,:) = zeros(1,size(spikeResp1,2));
% end

% scov = spikeResp1*spikeResp1';
% figure; imagesc(scov); colormap parula;
stimFileName = [stimFileName mosaicFile];
for incInd = 1%:length(includedComponentsArray)
    filterMat = linearReconstructSVD_short_midgets_both(stimFileName,spikeResp1,fileext, windowSize,includedComponentsArray(incInd),trainSizeArray(trainInd),shiftTime,stimType);
end

save([reconstructionRootPath '/dat/' filterFile],'filterMat','-v7.3');
