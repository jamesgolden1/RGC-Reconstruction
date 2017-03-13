% t_expDataReconstruction
% 
% Load data from a physiology experiment and put it through the
% linear reconstruction algorithm.
% 
% Training of decoder on line 167
% 
% 2/2017 JRG

%% Initialize

clear;
addpath(genpath(reconstructionRootPath));
% javaaddpath /home/vision/Nishal/Java/Java/vision7/bin/

% % % % 
% Description: The training data and testing data are interleaved in blocks
% during the experiment. The even blocks are a different, unique movie
% every time (the training data), and the odd blocks are a repeated movie
% (the testing data). The STA and the movies are saved [pixels, pixels,
% frames]. This folder contains two matfiles, one with spikes (CellData)
% and one with the movie (StimData), for each experiment and stimulus type:
% White Noise (WN) or Natural Scenes with Eye Movements (NSEM).
% 
% MAT files: CellData Spikes: Each cell has its spike times in seconds,
% separated into blocks. The even blocks correspond to the FitMovie, and
% the odd blocks correspond to the TestMovie. STA: Spike triggered average
% calculated from white noise StimData FitMovie: contains one cell for each
% fitting block TestMovie: contains one movie which is repeated every
% testing block
% 
% Note: The frame rate for our monitor is NOT 120Hz exactly. The actual
% time each frame is displayed for is: 0.00832750 s.

% piece = 'WN-2012-09-27-3';

% piece = 'WN-2013-08-19-6';
piece = 'NSEM-2013-08-19-6';
load(['/Volumes/Lab/Users/Nora/ShareData/CarlosData/' piece '-CellData.mat']);
load(['/Volumes/Lab/Users/Nora/ShareData/CarlosData/' piece '-StimData.mat']);

% load('/Volumes/Lab/Users/Nora/ShareData/CarlosData/NSEM-2013-08-19-6-CellData.mat');
% load('/Volumes/Lab/Users/Nora/ShareData/CarlosData/NSEM-2013-08-19-6-StimData.mat');

% 2013-08-19-6
% Isolated, used in Heitman and ICLR papers
% dt = .00832750
% NS: 59 1-minute training movies interleaved with 59 iterations of 30
% second test movie, skip first two iterations: 114 30-second movies: 91
% 30-second training movies, 23 30-second validation movies, 57 iterations
% of 30-second test movie
	
% WN: 60 30-second training movies interleaved with 60 iterations of
% 10-second test movie, skip first three iterations: 57 30-second movies:
% 46 30-second training movies, 11 30-second validation movies, 57
% iterations of 10-second test movie

if strcmpi(piece(1:2),'WN'); NSEMStimData = WNStimData; NSEMCellData = WNCellData; clear WNStimData; clear WNCellData; end;
movieLength= size(NSEMStimData.FitMovie{1},3);
for blockNum = 1:59
    fitMovie(:,:,(blockNum-1)*movieLength+1:blockNum*movieLength) = NSEMStimData.FitMovie{blockNum};
end

names = fieldnames(NSEMCellData);
cellSpikes = cell(length(names),1);
for nameInd = 1:length(names)
    blockNumCtr = 0;
    for blockNum = 2:2:118
        blockNumCtr = blockNumCtr+1;
        eval(['cellSpikesTemp = NSEMCellData.' names{nameInd} '.Spikes{blockNum};']);
        cellSpikes{nameInd} = [cellSpikes{nameInd}; (blockNumCtr-1)*movieLength*.00832750+cellSpikesTemp];
    end
end
%% Get STA

tstim = .00832750;%2/119.5172;
STA_length = 30;
movie_size = size(fitMovie);
STA = zeros(movie_size(1)*movie_size(2),STA_length,100);
fitframes = movie_size(3);

fitMovieRS = reshape(fitMovie,[movie_size(1)*movie_size(2),movie_size(3)]);
    
for cellind = 10%:length(cellSpikes)
    sp_frame = floor(cellSpikes{cellind}(:)/tstim);    
%     sp_rel = find((sp_frame>STA_length)&(sp_frame<fitframes));
    sp_rel = find(sp_frame>STA_length);
    

    for i = 1:STA_length
%         STA(:,i,cellind) = sum((fitMovieRSzm(:,(sp_frame(sp_rel)-STA_length+1)+i)),3);  

%     fitMovieMean = single(mean(fitMovieRS(:,sp_frame(sp_rel)-STA_length+1+i),2));
    fitMovieRSzm =-.5+ single(fitMovieRS(:,sp_frame(sp_rel)-STA_length+1+i));% - fitMovieMean*ones(1,size(sp_rel,1));
        STA(:,i,cellind) = sum(fitMovieRSzm,2);
    end
end

t0 = squeeze((sum(sum(STA,1),3)));
figure; plot(t0);
%%
STAview = reshape(STA(:,:,cellind),[movie_size(1),movie_size(2),STA_length]);
%  % Show STA movie
figure;
for i = 1:STA_length
   imagesc(STAview(:,:,i)')
%    imagesc(NSEMCellData.OFFPar_1006.STA(:,:,i)');
%    caxis([-14000 14000]);
   colormap gray
   axis image
   title('You should see an STA here')
   pause(0.1)
end
% figure; imagesc(NSEMCellData.OFFPar_1006.STA(:,:,5)')

%% Pull out spikes for each cell

% For only spatial decoding filters, get movie frame from peak of STA temporal response
% Loop over every cell
spikemat = zeros(length(cellSpikes),fitframes);
zshift =0;
for cellind = 1:length(cellSpikes)
for i = 1:length(cellSpikes{cellind})
    sp_frame = floor(cellSpikes{cellind}(i)/tstim);
    if sp_frame > STA_length && sp_frame<fitframes
        % STA = STA+double(fitmovie{1}(:,:,(sp_frame-STA_length+1):sp_frame));
        spikemat(cellind,(sp_frame-1)) = 1 + spikemat(cellind,(sp_frame-zshift));
    end
end
end

%%
% % Visualize spike covariance matrix
% figure; imagesc(spikemat*spikemat');
% % spikeResp = spikemat;
spikematmean = mean(spikemat,2);
spikematzm = spikemat - spikematmean*ones(1,fitframes);
% figure; imagesc(spikematzm*spikematzm');
spikematnorm = spikematzm./((std(spikematzm')'*ones(1,fitframes)));
spikematnormcov = spikematnorm*spikematnorm';
figure; imagesc(spikematnormcov - max(spikematnormcov(:))*eye(size(spikematnormcov,1))); colormap parula

%% Compute decoding filters

if ~isdir([reconstructionRootPath '/dat/' piece])
    mkdir(([reconstructionRootPath '/dat/' piece]));
    addpath([reconstructionRootPath '/dat/' piece]);
end

mosaicFile = '' ;
movieFileSave = [reconstructionRootPath '/dat/' piece '/movie']; 
spikesFileSave = [reconstructionRootPath '/dat/' piece '/spikes'];

stim = (reshape(fitMovie,40*80,size(fitMovie,3)));
spikeResp = spikemat;%(:,1:length(stim));

% % Visualize RF to make sure everything is still lined up
rf1 = (single(stim(:,1:end-5)-.5)*spikeResp(109,6:end)');
figure; imagesc(reshape(rf1,[80 40]))

% % Save for loading in future
save(movieFileSave,'stim');
save(spikesFileSave,'spikeResp');


movieFile = ['/' piece '/movie'];
spikesFile = ['/' piece '/spikes'];

% % % % 
pRecon.movieFile = [movieFile];
pRecon.spikesFile = [spikesFile];
pRecon.mosaicFile = mosaicFile;

% windowSize = 1 gives spatial decoding filters only
% set to value > 1 for spatiotemporal decoding filters
pRecon.windowSize = 10; 

evArr = [.01 .05 .1 .2 .4 .6 .8 .99];
trainFraction = [0.2 0.4 0.6 0.8];
zshift = 0;
cell_type = piece;

for evInd = 5%length(evArr)
    for trainFractionInd = length(trainFraction)
        
        if pRecon.windowSize == 1
            filterFile = ['' piece '/' cell_type '_sh' sprintf('%d',zshift) '_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd)) mosaicFile];
        else
            filterFile = ['' piece '/' cell_type '_sh' sprintf('%d',zshift) '_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd)) '_wind' sprintf('%d',pRecon.windowSize) mosaicFile];
        end
        
        pRecon.filterFile = filterFile;
        
        pRecon.percentSV = evArr(evInd);
        pRecon.trainFraction = trainFraction(trainFractionInd);
        pRecon.shiftTime = zshift;
        pRecon.stimType = 'wn';
        % RECONSTRUCTION CALL HERE
        [filterFile] = runReconstructSVD_fast_all(pRecon);
        % end
    end
end

% Load and visualize decoding filters
load(filterFile);

% % Visualize sum(abs(filters)) to see spatial extent
% % figure; imagesc(reshape(sum(abs(filterMat)),[80 40])')

numFilters = size(filterMat,1);
figure; for fr = 1:64; subplot(8,8,fr); imagesc(reshape(filterMat(00+fr,:),[80 40])'); colormap parula;  end;

figure; for fr = 1:25; subplot(5,5,fr); imagesc(reshape(filterMat(150+fr,:),[80 40])'); colormap parula;  end;

% numberCells = [103 112 117 182 12]; % for each cell type in '2016-02-17-6/data025'

%% Test decoding accuracy

% Uncomment to run only test
% cell_type_ind = 1;
% zshift=3;

% cell_str = {'on parasol','off parasol','on midget','off midget','on smooth','off smooth','all','all5'};
% cell_type = cell_str{cell_type_ind};

% For spatial only decoder, set tempFlag = 1 to add temporal STA course as temporal decoding filter timecourse
tempFlag = 0; 

mosaicFile = '' ;
% movieFileSave = [reconstructionRootPath '/dat/' piece '/movie']; 
% spikesFileSave = [reconstructionRootPath '/dat/' piece '/spikes'];

% Uncomment to run only test
% load(movieFileSave); 
% load(spikesFileSave);

mse = [];
% pRecon.windowSize = 1; % 1 = .4429, 4 = .4433

figure;

% Percent of eigenvalues/SVs to use
evArr = [.01 .05 .1 .2 .4 .6 .8 .99];
% Percent of total dataset for filters to have been trained on
trainFraction = [0.2 0.4 0.6 0.8];


stimTest = stim(:,ceil(0.8*length(stim)):length(stim));

spikeRespTest = zeros(1+size(spikeResp,1),length(ceil(0.8*length(stim)):length(stim)));

spikeRespTest = ones(1,length(ceil(0.8*length(stim)):length(stim)));
spikeRespTest(2:(1+size(spikeResp,1)),:) = spikeResp(:,ceil(0.8*length(stim)):length(stim));

% for evInd = length(evArr)
%     for trainFractionInd = length(trainFraction)
%         [evInd trainFractionInd]
                
        
        if pRecon.windowSize == 1
            filterFile = ['' piece '/' cell_type '_sh' sprintf('%d',zshift) '_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd)) mosaicFile];
        else
            filterFile = ['' piece '/' cell_type '_sh' sprintf('%d',zshift) '_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd)) '_wind' sprintf('%d',pRecon.windowSize) mosaicFile];
        end
        
        load(filterFile);
        
        if tempFlag 
             pRecon.windowSize = 12;
            load('t02.mat');
            filterMatSp = filterMat;
            filterMat = zeros(12*(size(filterMatSp,1)-1)+1,size(filterMatSp,2));
            
            cellTypeNumber = cumsum(2+numberCells);
            for ci = 2:size(filterMatSp,1)
                tind = min(find(cellTypeNumber>ci));
                filterTemp = filterMatSp(ci,:)'*(t0{1}-mean(t0{tind}))'./sum(abs(t0{1}-mean(t0{1})));
                filterMat(ci:ci+11,:) = filterTemp(:,19:30)';
            end
        end
        
        % stimRecon = -filterMat'*spikeRespTest;
        stimRecon =  reconsFromFiltLen(filterMat, spikeRespTest(2:end,:), pRecon.windowSize);
        % stimRecon = -filterMat(randperm(size(filterMat,1)),randperm(size(filterMat,2)))'*spikeRespTest;
        % stimRecon = -randn(size(filterMat))'*spikeRespTest;
        
        % figure; subplot(121); imagesc(stimTest(:,1:1000)); subplot(122); imagesc(stimRecon(:,1:1000))
        
        mc1rz = single(stimTest(:,1:end-pRecon.windowSize,1)) - ones(size(stimTest,1),1)*mean(stimTest(:,1:end-pRecon.windowSize,1));
        mc2rz = single(stimRecon) - ones(size(stimRecon,1),1)*mean(stimRecon,1);
         figure; subplot(121); imagesc(mc1rz(:,1:1000)); subplot(122); imagesc(mc2rz(:,1:1000))
    
        mc1rs = reshape(mc1rz,[80 40 size(mc1rz,2)]);
        mc2rs = reshape(mc2rz,[80 40  size(mc2rz,2)]);
        
%         mc1red = reshape(mc1rs(11:30,6:15,:),[10*20,size(stimTest,2)]); 
%         mc2red = reshape(mc2rs(11:30,6:15,:),[10*20,size(stimRecon,2)]); 
%         figure; imagesc(mc2rs(:,:,1));

        
        shiftTime = 0;
        errmov =(mc1rz(:,1+shiftTime:end))+1*(mc2rz(:,1:end-shiftTime));
        % errmov =(mc1rz(:,1+shiftTime:end-(pRecon.windowSize)))-1*(mc2rz(:,1:end-shiftTime));        
%         errmov =(mc1rs(:,1+shiftTime:end-(pRecon.windowSize)))-1*(mc2rs(:,1:end-shiftTime));
%         errmov =(mc1red(:,1+shiftTime:end-(pRecon.windowSize)))-1*(mc2red(:,1:end-shiftTime));
        errtot = ((errmov.^2));
        
        figure; subplot(121); imagesc(mc1rs(:,:,1)); subplot(122); imagesc(mc2rs(:,:,1))
         imv = 4; figure; subplot(121); imagesc(mc1rs(:,:,imv*80+1)); subplot(122); imagesc(mc2rs(:,:,imv*80+20));colormap gray %imagesc(sum(mc2rs(:,:,imv*80+[1:80]),3)); colormap gray
%         figure; subplot(131); hist(mc1rz(:),40); subplot(132); hist(mc2rz(:),40); subplot(133); hist(errmov(:),40);
        
        mse(evInd,trainFractionInd) = sqrt(mean(errtot(:)));
        mss(evInd,trainFractionInd) = sqrt(var(errtot(:)));
        
        trainSizeMat(evInd,trainFractionInd) = trainFraction(trainFractionInd);
        evMat(evInd,trainFractionInd) = evArr(evInd);
        
        clear mc1 mc2 errmov errtot movieRecon
        
%     end
% end

% 
% figure; 
% plot(1e0*trainSizeMat',mse','-x','linewidth',4)
% % hold on;
% % plot(1e0*trainSizeMat',mseh','-x','linewidth',4)
% grid on; % axis([0.2 1 .479 .485])
% xlabel('Training Set Size','fontsize',15); 
% ylabel('Mean Square Error - Reconstruction','fontsize',15);
% set(gca,'fontsize',13)
% legend('1% SVs','5% SVs','10% SVs','20% SVs','40% SVs','60% SVs','80% SVs','100% SVs');
% 
% fr = 50;
% figure; subplot(121); imagesc(reshape(stimTest(:,fr),40,20)); subplot(122);  imagesc(reshape(stimRecon(:,fr),40,20));
