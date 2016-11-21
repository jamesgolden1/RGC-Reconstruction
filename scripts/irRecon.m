function [movrecons_on_off_full, movrecons_on_off_dropout] = irRecon(innerRetina, filterName, varargin)
% Converts an RGC mosaic spike train to a reconstructed movie.
% 
% testReconAll loads an RGC mosaic as well as a movie input and computes 
% the spike response. The spike response is passed to this function, where
% the spikes are pulled out of the object, rearranged and linearly decoded
% to produce a reconstructed movie.
% 
% 10/2016 JRG (c) isetbio team
% 
% See also testReconAll, runReconstruct, irOptimalRecon
%%

p = inputParser;
p.addRequired('innerRetina',@isa('innerRetina'));
p.addRequired('filterFile');
p.addParameter('percentDead',0,@isnumeric);
p.parse(innerRetina,filterFile,varargin{:});
innerRetina = p.Results.innerRetina;
filterFile = p.Results.filterFile;
mosaicFile = p.Results.filterFile;
percentDead = p.Results.percentDead;

%% Get spikes from mosaic objects

y = cell(length(innerRetina.mosaic));

% Loop over the mosaic objects
for mosaicInd = 1:length(innerRetina.mosaic)
    
    cellCtr=0; dt = .01;
    maxTrials = innerRetina.mosaic{1}.numberTrials;
    nCells = size(innerRetina.mosaic{mosaicInd}.responseSpikes);
    
    y{mosaicInd} = zeros(nCells(1)*nCells(2),10000);
    
    for ycell = 1:nCells(1)        
        for xcell = 1:nCells(2)
            cellCtr = cellCtr+1;            
            for trial = 1%:maxTrials                
                yind =  innerRetina.mosaic{mosaicInd}.responseSpikes{xcell,ycell,trial,1};                
                y{mosaicInd}(cellCtr,ceil(yind./dt))=1;                
            end
        end
    end    
end


%% Load the reconstruction filters, locally or from RDT

% all four mosaics, nov 4
% load('/Users/james/Downloads/filters_mosaic_ns_all_overlap0_svd_1000_len_100.mat')
load('/Users/james/Downloads/filters_mosaic_wn_all_42reps_overlap0_svd_1000_len_100.mat')
% rdt = RdtClient('isetbio');
% rdt.crp('/resources/data/reconstruction');
% data = rdt.readArtifact('filters_may26_parasol_midget_combined_svd_3000_len_100', 'type', 'mat');
% filterMat = data.filterMat; clear data;

% Reshape spikes and decode

recons_stim_on_off = cell(length(innerRetina.mosaic),1);
movrecons_on_off = cell(length(innerRetina.mosaic),1);

recons_stim_on_off_dropout = cell(length(innerRetina.mosaic),1);
movrecons_on_off_dropout = cell(length(innerRetina.mosaic),1);


blocklength = 100;
sizeSpikes = floor([ size(y{1},2) size(y{2},2) size(y{3},2) size(y{4},2)  ]/blocklength);

% on Parasol
spikesout = (y{1});%double(matfOff.spikesoutsm);
numcells1 = 36;

spikeRespOn = zeros(size(spikesout,1), max(sizeSpikes));

spikeRespOn(:,1:floor(size(spikesout,2)/blocklength)) = downSampResp(spikesout, numcells1, floor(size(spikesout,2)/blocklength));
%     spikeRespOn= zeros(size(downSampResp(spikesout, numcells1, blocklength)));

% off Parasol
numcells2 = 64;
spikesout2 = (y{2});%double(matfOff.spikesoutsm);
spikeRespOff = zeros(size(spikesout2,1), max(sizeSpikes));
spikeRespOff(:,1:floor(size(spikesout2,2)/blocklength)) = ((downSampResp(spikesout2, numcells2, floor(size(spikesout2,2)/blocklength))));
%     spikeRespOff = zeros(size(downSampResp(spikesout2, numcells2, blocklength)));

% off Midget
numcells3 = 225;
spikesout3 = (y{3});%double(matfOff.spikesoutsm);

spikeRespOffM = zeros(size(spikesout3,1), max(sizeSpikes));
spikeRespOffM(:,1:floor(size(spikesout3,2)/blocklength)) = ((downSampResp(spikesout3, numcells3, floor(size(spikesout3,2)/blocklength))));
%     spikeRespOffM = zeros(size(downSampResp(spikesout3, numcells3, blocklength)));

% on Midget
numcells4 = 169;
spikesout4 = (y{4});%double(matfOff.spikesoutsm);

spikeRespOnM = zeros(size(spikesout4,1), max(sizeSpikes));
spikeRespOnM(:,1:floor(size(spikesout4,2)/blocklength)) = ((downSampResp(spikesout4, numcells4, floor(size(spikesout4,2)/blocklength))));
%     spikeRespOnM = zeros(size(downSampResp(spikesout4, numcells4, blocklength)));

% spikesout = vertcat(onSR(:,1:15000), offSR(:,1:15000), onPSR(:,1:15000), offPSR(:,1:15000));

% whole mosaic
spikeRespOnOff =vertcat(spikeRespOn,spikeRespOff, spikeRespOffM,spikeRespOnM);
% only parasols
%     spikeRespOnOff =vertcat(spikeRespOn,spikeRespOn, zeros(size(spikeRespOffM)),zeros(size(spikeRespOnM)));
% only midgets
%     spikeRespOnOff =vertcat(zeros(size(spikeRespOn)),zeros(size(spikeRespOff)), ((spikeRespOffM)),((spikeRespOnM)));

numbins = 8;
recons_stim_on_off{mosaicInd,1} = reconsFromFiltLen(filterMat, spikeRespOnOff, numbins);
% filterMatInd = find(abs(filterMat)<0.001); filterMat2 = filterMat; filterMat2(filterMatInd)=0;
% recons_stim_on_off{mosaicInd,1} = reconsFromFiltLen(filterMat2, spikeRespOnOff, numbins);
movrecons_on_off_full = reshape(recons_stim_on_off{mosaicInd,1},96,96,size(recons_stim_on_off{mosaicInd,1},2));

%% Reconstruction with a random subset of spiking cells zeroed out
deadIndicesAll = randperm(size(spikeRespOnOff,1));
numberDead = round(percentDead*(size(spikeRespOnOff,1)));
deadIndices = deadIndicesAll(1:numberDead);
spikeRespOnOffDead = spikeRespOnOff;
spikeRespOnOffDead(deadIndices,:) = zeros(length(deadIndices),size(spikeRespOnOff,2));

recons_stim_on_off_dropout{mosaicInd,1} = reconsFromFiltLen(filterMat, spikeRespOnOffDead, numbins);

movrecons_on_off_dropout = reshape(recons_stim_on_off_dropout{mosaicInd,1},96,96,size(recons_stim_on_off{mosaicInd,1},2));
% movrecons_on_off = .5*movrecons_on_off_full./mean(movrecons_on_off_full(:));
