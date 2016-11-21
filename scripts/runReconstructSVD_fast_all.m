function filterFile = runReconstructSVD_fast_all(varargin)
% Run the reconstruction algorithm to generate the filters

% Must run loadSpikesAll.m to build stimFile and spikesFile

p = inputParser;
p.addParameter('movieFile',[],@ischar);
p.addParameter('spikesFile',[],@ischar);
p.addParameter('filterFile',[],@ischar);
p.parse(varargin{:});
loadFile = p.Results.loadFile;
movieFile = p.Results.movieFile;
spikesFile = p.Results.spikesFile;

if isempty(filterFile)
    filterFile = ['filters_' num2str(round(cputime*100))];
end

%% run Linear reconstruction for on, off, and joint on/off training 
windowsize = 8;
disp('Loading spike responses...')

% stimFileName = 'NSmovie_spikeResp_overlap0';
% spikesFileName = 'NSspikeResp_overlap0';

stimFileName = 'NSmovie_40reps_ns0';
spikesFileName = 'NSspikeResp_40reps_ns0';
% matfON = matfile([reconstructionRootPath '\dat\' spikesFileName]);
matfON = matfile([reconstructionRootPath '/dat/' spikesFileName]);

movielength = 1*240000;%size(stim,2);
disp(['Total Movie Length in Frames: ' num2str(movielength)]);

fileext = 'mosaic_ns_all_40reps_ns0_rz';
trainSizeArray = 1;%[.6/8:.6/8:.6]; 
trainInd = 1;
includedComponentsArray = 1000;

srON = matfON.spikeResp;
% srOFF = matfOFF.spikeResp;
% clear matfON matfOFF
% spikeResp1 = vertcat(srON(:,1:12*12000));%, srOFF(:,1:240000), matfONP.spikeResp, matfOFFP.spikeResp);
spikeResp1 = srON;

% scov = spikeResp1*spikeResp1';
% figure; imagesc(scov); colormap parula;

for incInd = 1%:length(includedComponentsArray)
    filterMat = linearReconstructSVD_short_midgets_both(stimFileName,spikeResp1(:,:),fileext, windowsize,includedComponentsArray(incInd),trainSizeArray(trainInd));
end

save([reconstructionRootPath '/dat/' filterFile],'filterMat','-v7.3');

% figure; imagesc(reshape(filterMat(2+1*8,:),96,96));

%%
cellNumber = 40;
figure;
for fr = 1:8
    starf(:,:,fr) = reshape(filterMat(fr+cellNumber*8,:),96,96);
    subplot(2,4,fr); imagesc(starf(:,:,fr)); 
%     caxis([-3e-3 3e-3]);
    caxis([-1e-4 1e-4]);
    colormap parula;
end
% figure; ieMovie(starf);

%% on parasol

% load([reconstructionRootPath '\dat\'  'filters_mosaic_all_nat_onParasol_svd_288_len_100']);

% load([reconstructionRootPath '\dat\'  'filters_mosaic_all_nat_onParasol_svd_1000_len_100']);
% load([reconstructionRootPath '\dat\'  'filters_mosaic_all_nat_onParasol_svd_3952_len_100']);

figure;
for cellNumber = 1:36
    subplot(6,6,cellNumber);
    imagesc(reshape(filterMat(2+cellNumber*8,:),96,96)); 
%     caxis([-.5e-3 .5e-3]); 
%     caxis([-1e-3 1e-3]);
%     caxis([-9e-6 9e-6]);
caxis([-.35 .35])
axis off
    colormap parula
end
%% off parasol
figure;
for cellNumber = 37:37+64    
% for cellNumber = 1:64
    subplot(8,8,cellNumber-36);
    imagesc(reshape(filterMat(2+cellNumber*8,:),96,96)); 
%     caxis([-.5e-3 .5e-3]); 
%     caxis([-2e-3 2e-3]);

% caxis([-8e-6 8e-6]);
axis off
    colormap parula
end
%% off midget
figure;
startCell = 65;
totalCells = 8;
for cellNumber = startCell+36+64+1:startCell+36+64+1+totalCells^2%225
%     subplot(15,15,cellNumber-(36+64));
    subplot(totalCells,totalCells,cellNumber-(startCell+36+64));
    imagesc(reshape(filterMat(2+cellNumber*8,:),96,96)); 
%     caxis([-.125e-3 .125e-3]); 
%     caxis([-4e-3 4e-3]);
caxis([-6e-6 6e-6]);
% caxis([-2e-3 2e-3]);
%  caxis([-1e-3 1e-3]);
    colormap parula; %axis square
    axis off
%     drawnow;
end
%% on midget
figure;
startCell = 1;
totalCells = 8;
cMax = max(max(filterMat(2+startCell+36+64+225+1:startCell+36+64+1+225+totalCells^2*8,:)));
for cellNumber = startCell+36+64+225+1:startCell+36+64+1+225+totalCells^2%225
%     subplot(15,15,cellNumber-(36+64));
    subplot(totalCells,totalCells,cellNumber-(startCell+36+64+225));
    imagesc(reshape(filterMat(2+cellNumber*8,:),96,96)); 
%     caxis([-.125e-3 .125e-3]); 
%     caxis([-1e-3 1e-3]);
caxis([-cMax cMacx]);
% caxis([-9e-6 9e-6]);
    colormap parula; %axis square
%     drawnow;
axis off
end