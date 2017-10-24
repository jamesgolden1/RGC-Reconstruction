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
p.addParameter('dropout',0,@isnumeric);
p.addParameter('stimType','ns',@ischar);
p.addParameter('buildFile',[],@ischar);

p.addParameter('pixelWidth',70,@isnumeric);
p.addParameter('currentDecay',2,@isnumeric);

p.addParameter('numTest',1,@isnumeric);
p.KeepUnmatched = true;
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
buildFile = p.Results.buildFile;
dropout = p.Results.dropout;

pixelWidth = p.Results.pixelWidth;
currentDecay = p.Results.currentDecay;
numTest = p.Results.numTest;

if isempty(filterFile)
    filterFile = ['filters_' num2str(round(cputime*100))];
end
% % 

stimFileTest = [ stimFileName '_' mosaicFile];

% stimFileTest = [ stimFileName '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile];

if exist([stimFileTest '.mat'],'file') ~= 2
    loadSpikes(buildFile,stimFileName,respFileName,mosaicFile,pixelWidth,currentDecay,numTest)
end


% if exist([stimFileName mosaicFile '.mat'],'file') ~= 2
%     loadSpikes(buildFile,stimFileName,respFileName,mosaicFile)
% end
%% run Linear reconstruction for on, off, and joint on/off training 

disp('Loading spike responses...')

% stimFileName = 'NSmovie_spikeResp_overlap0';
% respFileName = 'NSspikeResp_overlap0';

% stimFileName = 'NSmovie_40reps_ns0';
% respFileName = 'NSspikeResp_40reps_ns0';

% matfON = matfile([reconstructionRootPath '\dat\' respFileName]);

matfON = matfile([reconstructionRootPath '/dat/' respFileName '_' mosaicFile]);
% matfON = matfile([reconstructionRootPath '/dat/' respFileName '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat']);

% matfON = matfile('/Volumes/Lab/Users/james/RGC-Reconstruction/dat/ns100_r2_10/ns100_jan1_sp3_mosaicAll_1246640');

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
% stimFileName = [stimFileName mosaicFile];

% stimFileName = [ stimFileName '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile];

stimFileName = [ stimFileName '_' mosaicFile];
for incInd = 1%:length(includedComponentsArray)
    [filterMat,dropoutIndices] = linearReconstructSVD_short_midgets_both(stimFileName,spikeResp1,respFileName, windowSize,includedComponentsArray(incInd),trainSizeArray(trainInd),shiftTime,stimType,dropout);
end

save([reconstructionRootPath '/dat/' filterFile '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay)],'filterMat','dropoutIndices','-v7.3');


function loadSpikes(buildFile,stimFile,respFile,mosaicFile,pixelWidth,currentDecay,numTest)
loadFile = buildFile;

%% loadSpikesAll %%%%5
% % % load('test0_block_1_mosaicAll_23282.mat')
%
if isempty(stimFile)
    stimFile = ['movie_' num2str(round(cputime*100))];
end

if isempty(respFile)
    respFile = ['spikes_' num2str(round(cputime*100))];
end

loadFile = buildFile;

for fi = 1:numTest
 
if ismac || isunix
%     dNames = (dir([reconstructionRootPath '/dat/' loadFile '*block_*' mosaicFile '.mat']));
    dNames = (dir([reconstructionRootPath '/dat/' loadFile '*block_*' '_start_' num2str(fi) '_' mosaicFile '.mat']));
%  dNames = (dir([reconstructionRootPath '/dat/' loadFile '*block_*' '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat']));
else
%     dNames = (dir([phospheneRootPath '\dat\' loadFile '*block_*' mosaicFile '.mat']));
 dNames = (dir([reconstructionRootPath '\dat\' loadFile '*block_*' mosaicFile '.mat']));
end
% blocklength = 12000;
numReps = length(dNames)
% numCells= 36+64+169+225;
% spikeResp = zeros(numCells, blocklength*numReps);

if isempty(strfind(dNames(1).name,'469'))
    blockStart = 0;
else
    blockStart = 468;
end

  % filename1 = [reconstructionRootPath '/dat/' buildFile '_block_' num2str(1) '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat'];
if isunix || ismac
%     filename1 = [phospheneRootPath '/dat/' loadFile '_block_' num2str(1) mosaicFile '.mat'];
%     filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(blockStart+1) '_' mosaicFile '.mat'];
filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(blockStart+1) '_start_' num2str(fi) '_' mosaicFile '.mat'];

%     filename1 = [reconstructionRootPath '/dat/' loadFile dNames(1).name];
else
%     filename1 = [phospheneRootPath '\dat\' loadFile '_block_' num2str(1) mosaicFile '.mat'];
filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(1) mosaicFile '.mat'];
end

matf = matfile(filename1);
szMov = size(matf.whiteNoiseSmall);
blocklength = szMov(3);
stim = zeros(szMov(1)*szMov(2),blocklength*numReps,'uint8');
clear matf
blockNum = 0+blockStart;
for blockNumInd =[1:length(dNames) ]
% for blockNumInd =[1:12 21:50]
    blockNum = blockNum+1
    % filename1 = [reconstructionRootPath '\dat\WNstim_response_stx2_block_' num2str(blockNum) '.mat'];    
    % filename1 = [reconstructionRootPath '\dat\WNstim_response_block_' num2str(blockNumInd) '.mat'];    
    % filename1 = [reconstructionRootPath '/dat/nsResponses/NSstim_response_betha_ns0_block_' num2str(blockNum) '.mat'];    
    % filename1 = [reconstructionRootPath '\dat\NSstim_response_overlap0_block_' num2str(blockNum) '.mat'];

%   filename1 = [reconstructionRootPath '/dat/' buildFile '_block_' num2str(blockNum) '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat'];
    if isunix || ismac
%         filename1 = [phospheneRootPath '/dat/' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
%  filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(blockNum) '_' mosaicFile '.mat'];
   filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(blockNum) '_start_' num2str(fi) '_' mosaicFile '.mat'];
 
    else
%         filename1 = [phospheneRootPath '\dat\' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
 filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(blockNum) mosaicFile '.mat'];
    end
    
    matf = matfile(filename1);
    spikesoutsm = matf.spikesoutsm;
    % Spikes in this variable for each block
    spikesout = double(spikesoutsm);
    pointer = (blockNum-1-blockStart)*blocklength;
    
    for i = 1:blocklength
        blocksize = 10;
        endval = i*blocksize;
        if endval > size(spikesout,2)
            endval = size(spikesout,2);
        end
        startval = (i-1)*blocksize + 1;
        spikeResp(:,pointer+i) = sum(spikesout(:,startval:endval),2);
    end
    clear spikesout
    clear spikesoutsm
    szMov = size(matf.whiteNoiseSmall);
    stimtmp = reshape(matf.whiteNoiseSmall,szMov(1)*szMov(2),blocklength);
%     stimtmp = uint8(128+double(stimtmp) - ones(size(stimtmp,1),1)*mean(stimtmp,1));
    stim(:,(blockNum-1-blockStart)*blocklength +1 : (blockNum-blockStart)*blocklength) = stimtmp;

    % Stimulus here
    % whiteNoiseSmall;
end

% stim = (single(stim')-(ones(size(stim,2),1)*mean(stim,2)'));
% stim = stim';

% stim = (single(stim)-(ones(size(stim,2),1)*mean(stim,2)')');
% stimzm = (single(stim)-(ones(size(stim,2),1)*mean(stim,2)')');
% clear stim; stim = stimzm;
    

%     save([reconstructionRootPath '/dat/' respFile '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat'],'spikeResp','-v7.3');
%     save([reconstructionRootPath '/dat/' stimFile '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '_' mosaicFile '.mat'],'stim','-v7.3')

if isunix || ismac
%     
    save([reconstructionRootPath '/dat/' respFile '_' mosaicFile '.mat'],'spikeResp','-v7.3');
    save([reconstructionRootPath '/dat/' stimFile '_' mosaicFile '.mat'],'stim','-v7.3')
    
%     save([reconstructionRootPath '/dat/' respFile '_' num2str(fi) '_' mosaicFile '.mat'],'spikeResp','-v7.3');
%     save([reconstructionRootPath '/dat/' stimFile '_' num2str(fi) '_' mosaicFile '.mat'],'stim','-v7.3')
else
    
    save([reconstructionRootPath '\dat\' respFile mosaicFile],'spikeResp','-v7.3');
    save([reconstructionRootPath '\dat\' stimFile mosaicFile],'stim','-v7.3')
end

end
