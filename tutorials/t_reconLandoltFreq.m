

% scp vision@bertha.stanford.edu:/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/healthy_gratings_aug16_cont/cont/* '/Users/james/Documents/Matlab/RGC-Reconstruction/dat/healthy_gratings_aug16_cont/'


%%

clear

% addpath(genpath('/Users/james/Documents/MATLAB/RemoteDataToolbox'))
% addpath(genpath('/Users/james/Documents/MATLAB/EJLPhosphene'))
% addpath(genpath('/Users/james/Documents/MATLAB/isetbio'))
% addpath(genpath('/Users/james/Documents/MATLAB/RGC-Reconstruction'))


% degRotArr = [180 178 175 170 165 160 155 150 145 140 135 130 125 120 100 90 45 20 10 0];

% degRotArr = [180 179 178 177 176.5 176 175.75 175.5 175.25 175 172];
degRotArr = 0%[180 179.5 175:-2:164 162:-3:150];
% pixelRadArr = [5:5:20 22:28 30 35 43 50];
% pixelRadArr = [ 10 20 23 25 27 28 29 30 40 50];
pixelRadArr = 50;%[5 17 23 24 25 26 28 29 45 ];
noiseVal = 0;

nBasis = 200;

% [xi,yi] = meshgrid(-49:50,-49:50);
% radVal = (sqrt(xi.^2+yi.^2));
% radInd = find(radVal>45);
% radVal(radInd) = 0;
% figure; imagesc(radVal)

freqArr =2;%[.05 .1 .2 .5 1 2 4 5 8 10 16];

% sep 14
% contrastArr = 1*[0 .00002 .0002 .0004 .001 .002:.001:.006 .008 .02:.01:.06];
% sep 15
% contrastArr = 1*[0 .0002 .0004 .001 .002:.002:.01 .02:.02:.1 .15 .2];

contrastArr = .04;%1*[0 .0001 .001 .008 .01 .0125 .015 .0175 .02 .0225 .025 .0275 .03 .05 .1 .15 ];
freqArr = sqrt([ .01 .1 .2 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .875 1]);
freqArr = [ freqArr sqrt([ .001 .005 .008 .02 .025 .03 .04 .06 .08 .09 .11 .12 .15])];% .2 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .875 1]);

freqArr = sort(freqArr,'ascend');
% dirname = 'gratings_sep17/prima_sep17/';
% dirname1 = 'august_cont1/healthy_landolt_aug16_nogap';
% dirname2 = 'august_cont1/healthy_landolt_aug16_gap';
dirname1 = 'healthy_landolt_aug16_nogap_long';
dirname2 = 'healthy_landolt_aug16_gap_long';
numReps =4 ;

maxTrial = 150;%size(framesOut,3);

% contrastArr = [.002 .003 .004 .005 .006 .007];
% 
% dirname = 'gratings_sep17/prima_sep17add/';
% numReps =16 ;
% maxTrial = 100;

plotSingle = 0;

% contrastArr = [ .0008 .0009 .003 .005 .007 .009  ];
% dirname = 'oct3/oct3/prima_sep19add';
figure; 


% For only on, use healthy with lambda = .05
% % load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0_zero05.mat')
% % load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/prosthesis_35_training_aug13/filtersmosaic0_sv 5_w1_sh3_dr30_pitch_35_decay_2.mat')
% % load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/prosthesis_35_training_aug13/filtersmosaic0_sv 5_w1_sh3_dr30_pitch_35_decay_2.mat')
% % load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0_zero10.mat')
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0_zero05.mat')
% % load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0.mat'); filterMat2 = filterMat; clear filterMat;
% % lambda = .06;
% filterMat = load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/prosthesis_70_training_aug13/filtersmosaic0_sv 5_w1_sh3_dr30_pitch_70_decay_2.mat')

% if isstruct(filterMat)
%     droputIndices = filterMat.dropoutIndices;
%     filterMatSm = filterMat.filterMat;
%     filterMat=zeros(9717,10000);
%     filterMat(droputIndices,:)=filterMatSm;
%     filterMatSm=[];        
% 
% %         spikeAugFull = spikeAug;
% %         spikeAug = zeros(size(spikeAug));
% %         spikeAug(droputIndices,:) = spikeAugFull(droputIndices,:);
% end
% filterMat2 = filterMat; filterMat = [];
% % filterMat2 = zeroFilter(filterMat,lambda);
figure;
for contrastInd =1:length(contrastArr)%-8:length(contrastArr)%length(contrastArr)-15:length(contrastArr)-3
for freqInd =  1:length(freqArr)
% for degInd = 1:length(degRotArr)
for degInd = length(pixelRadArr)

%     d0 = dir('/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/prima10/* 0.mat');
% 
% load(['/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/prima10/' d0(1).name]);
% d1 = dir('/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/prima10/*.mat');
% 
% load(['/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/prima10/' d1(ord1(1)).name])

% load(fullfile('/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/primaFreqH/',[sprintf('primaRecon_freqH%4d.mat',100*freqArr(freqInd))]),'trialReconPlay');

% load(fullfile('/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/primaFreqH/',[sprintf('primaRecon_freqH%4d.mat',100*freqArr(freqInd))]),'trialReconPlay');
% %     save(fullfile(reconstructionRootPath,'dat',folderNameTest,'cont',[sprintf('primaRecon_freqH%4d_cont%6d.mat',100*gratingSpFreq,1e5*contrast)]),'trialReconPlay','contrast');

% load(fullfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/gratings/gratings_sep134_100/',[sprintf('recon_freq%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(contrastInd))]),'trialReconPlay','contrast');

% load(fullfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/landolt_sep/sep15/healthy/',[sprintf('recon_freq%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(contrastInd))]),'trialReconPlay','contrast');
% load(fullfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/landolt_sep/sep15/healthy/',[sprintf('recon_freq%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(contrastInd))]),'trialReconPlay','contrast','spikeTrials');

load(fullfile(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/' dirname1],[sprintf('recon_freqH%4d_cont%6d_rs%2d.mat',100*25,1e5*contrastArr(contrastInd),100*freqArr(freqInd))]),'trialReconPlay','contrast','spikeTrials');

% trialReconPlayXW = filterMat2'*spikeTrials;
% trialReconPlay2 = reshape(trialReconPlayXW,[100 100 size(trialReconPlayXW,2)]);


spikeOn = spikeTrials;
% spikeOn = zeros(size(spikeTrials));
% spikeOn(1:1+28*32,:) = spikeTrials(1:1+28*32,:);
% spikeOn(28*32+31*35+1:28*32+31*35+55*63+1,:) = spikeTrials(28*32+31*35+1:28*32+31*35+55*63+1,:);

trialReconPlayXW = filterMat2'*spikeOn;

framesOut = reshape(trialReconPlayXW,[100 100 size(trialReconPlayXW,2)]);

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/primaRecon_one.mat');
% framesOut = trialReconPlay;%2;
% dataRaw = single(zeros(size(framesOut,1)*size(framesOut,2), 2*size(framesOut,3)));
% dataRaw(1:size(framesOut,1)*size(framesOut,2),1:size(framesOut,3)) = RGB2XWFormat(framesOut);


% centerValsR = [1:size(framesOut,1)];
% centerValsC = [1:size(framesOut,2)];

centerValsR = [50-pixelRadArr(degInd)+1:50+pixelRadArr(degInd)];
centerValsC = [50-pixelRadArr(degInd)+1:50+pixelRadArr(degInd)];

dataRaw = single(zeros(length(centerValsR)*length(centerValsC), 2*maxTrial));
% dataRaw(1:size(framesOut,1)*size(framesOut,2),1:size(framesOut,3)) = RGB2XWFormat(framesOut);


%%

for fr = 1:maxTrial
    imr1 = framesOut(centerValsR,centerValsC,fr);
%     imr1 = imrotate(framesOut(:,:,fr),noiseVal*(-1+2*rand(1,1)),'crop');
%     zInd = find(imr1==0);
%     imr1mid = imr1(41:60,41:60);
%     imr1(zInd) = mean(imr1mid(:))*ones(length(zInd),1) + std(imr1mid(:))*randn(length(zInd),1);
    
    % dataRot(zInd,fr) = mean(imr1mid(:))*ones(length(zInd),1) + std(imr1mid(:))*randn(length(zInd),1);
    % figure; imagesc(imr1);
    dataRot(:,fr) = RGB2XWFormat(imr1); 
end
% dataRot(radInd,:) = 0;

% figure; imagesc(reshape(dataRot(:,10),[100 100]))

% dataRaw(1:size(framesOut,1)*size(framesOut,2),[1:size(framesOut,3)]) = dataRot;

dataRaw(1:length(centerValsR)*length(centerValsC),[1:maxTrial]) = dataRot(:,1:maxTrial);

clear framesOut trialReconPlay trialReconPlay2 spikeTrials


%%

% d1 = dir('/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/prima10/*.mat');
% 
% load(['/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/prima10/' d1(ord1(freqInd)).name])

% load(fullfile('/Users/james/Documents/MATLAB/EJLPhosphene/local/paper/july17_gratings/primaFreq/',[sprintf('primaRecon_freq%4d.mat',100*freqArr(freqInd))]),'trialReconPlay');

% load(fullfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/gratings/gratings_sep134_100/',[sprintf('recon_freqH%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(contrastInd))]),'trialReconPlay','contrast');
% load(fullfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/landolt_sep/sep15/healthy/',[sprintf('recon_freqH%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(contrastInd))]),'trialReconPlay','contrast','spikeTrials');
% load(fullfile(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/' dirname],[sprintf('primaRecon_freqH%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(contrastInd))]),'trialReconPlay','contrast','spikeTrials');
% load(fullfile(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/' dirname2],[sprintf('primaRecon_freq%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(contrastInd))]),'trialReconPlay','contrast','spikeTrials');
load(fullfile(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/' dirname2],[sprintf('recon_freq%4d_cont%6d_rs%2d.mat',100*25,1e5*contrastArr(contrastInd),100*freqArr(freqInd))]),'trialReconPlay','contrast','spikeTrials');

% load(fullfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/gratings/gratings_sep13/',[sprintf('recon_freq%4d_cont%6d.mat',100*freqArr(freqInd),1e5*contrastArr(1))]),'trialReconPlay','contrast');

spikeOn = spikeTrials;
% spikeOn = zeros(size(spikeTrials));
% spikeOn(1:1+28*32,:) = spikeTrials(1:1+28*32,:);
% spikeOn(28*32+31*35+1:28*32+31*35+55*63+1,:) = spikeTrials(28*32+31*35+1:28*32+31*35+55*63+1,:);

trialReconPlayXW = filterMat2'*spikeOn;

framesOut = reshape(trialReconPlayXW,[100 100 size(trialReconPlayXW,2)]);
% trialReconPlayXW = filterMat2'*spikeTrials;
% trialReconPlay2 = reshape(trialReconPlayXW,[100 100 size(trialReconPlayXW,2)]);

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/primaRecon_zero.mat');
% framesOut = trialReconPlay;%2;
% dataRaw(1:size(framesOut,1)*size(framesOut,2),size(framesOut,3)+[1:size(framesOut,3)]) = RGB2XWFormat(framesOut);


% imr1 = imrotate(framesOut(:,:,1),45,'crop');
% zInd = find(imr1==0); 
% imr1mid = imr1(41:60,41:60);
% imr1(zInd) = mean(imr1mid(:))*ones(length(zInd),1) + std(imr1mid(:))*randn(length(zInd),1); 
% figure; imagesc(imr1);

%%

% degRot = 178;

for fr = 1:maxTrial
    
    imr1 = framesOut(centerValsR,centerValsC,fr);
%     imr1 = imrotate(flipud(framesOut(:,:,fr)),degRotArr(degInd) + noiseVal*(-1+2*rand(1,1)),'crop');
%     imr1 = imrotate(fliplr(flipud((framesOut(:,:,fr)))),degRotArr(degInd) + noiseVal*(-1+2*rand(1,1)),'crop');
%     zInd = find(imr1==0);
%     imr1mid = imr1(41:60,41:60);
%     imr1(zInd) = 0;%mean(imr1mid(:))*ones(length(zInd),1) + std(imr1mid(:))*randn(length(zInd),1);
    
%     dataRot(zInd,fr) = mean(imr1mid(:))*ones(length(zInd),1) + std(imr1mid(:))*randn(length(zInd),1);
    % figure; imagesc(imr1);
    dataRot(:,fr) = RGB2XWFormat(imr1); 
end
% dataRot(radInd,:) = 0;
% figure; imagesc(reshape(dataRot(:,10),[100 100]))

% dataRaw(1:size(framesOut,1)*size(framesOut,2),size(framesOut,3)+[1:size(framesOut,3)]) = dataRot;

dataRaw(1:length(centerValsR)*length(centerValsC),maxTrial+[1:maxTrial]) = dataRot(:,1:maxTrial);

clear framesOut trialReconPlay
% figure; subplot(121);imagesc(reshape(mean(dataRaw(:,1:200),2),[100 100])); 
% subplot(122); imagesc(reshape(mean(dataRaw(:,200+[1:200]),2),[100 100])); 
% figure; subplot(121); imagesc(reshape(dataRaw(:,10),[100 100])); colormap gray
% subplot(122); imagesc(reshape(dataRaw(:,200+10),[100 100])); colormap gray;
%
% % dataCov = dataRaw*dataRaw'; figure; imagesc(dataCov(1:500,1:500));

% [~,~,V] = svd(dataRaw','econ');
% dataBasis = V(:,1:nBasis);
% 
% % figure;
% % for i = 1:25
% %     subplot(3,3,i)
% %     imagesc(reshape(dataBasis(:,i),[100 100])); colormap gray;
% % end
% 
% % Time series of weights
% data  = dataRaw' * dataBasis(:,[1:nBasis]);

data = dataRaw';
nTrials = size(data,1)/2;

% figure; hold on; grid on;
% scatter3(data(1:nTrials,2),data(1:nTrials,3),data(1:nTrials,7),'r');
% scatter3(data(nTrials+[1:nTrials],2),data(nTrials+[1:nTrials],3),data(nTrials+[1:nTrials],7),'b');
% % title(sprintf('PC Space, %d Degrees Separation', 180-degRotArr(degInd)));
% title(sprintf('%1.4f',contrast));
% xlabel('PC 2'); ylabel('PC 3'); zlabel('PC 4');
% set(gca,'fontsize',14);

%% Start classification training

for rep = 1:numReps

% Put the weights from each trial into the rows of a matrix
% Each row is another trial
nWeights = size(data,2);
% data = zeros(2*nTrials,nWeights*tSamples);
% for ii = 1 : (2*nTrials)
%     start = (ii-1)*tSamples + 1;
%     thisTrial = weightSeries(start:(start+tSamples - 1),:);
%     data(ii,:) = thisTrial(:)';
% end
label = [ones(nTrials, 1); -ones(nTrials, 1)];

% Select some of the data (80%) as the training set.
train_index = zeros(nTrials, 1);
train_index(randperm(nTrials, round(0.8*nTrials))) = 1;
train_index = train_index > 0;

% The aligned and offset trials are still matched
train_index = repmat(train_index, 2, 1);

% Fit the SVM model.
mdl = fitcsvm(data(train_index, :), label(train_index), ...
    'KernelFunction', 'linear');

% predict the data not in the training set.
yp = predict(mdl, data(~train_index, :));
classLoss = sum(label(~train_index) ~= yp) / length(yp);

% X(bb) = barOffset(bb);
P(degInd,rep) = (1-classLoss) * 100;
clear label train_index mdl yp classLoss
end
[freqArr(freqInd) mean(P')]

% P
imagesc(dataRaw); drawnow;
clear dataRaw data dataBasis V dataRot mdl label yp
end


% Pbig(:,:,freqInd) = P;

Pbig(:,:,freqInd) = P;
% 
% figure;
% imagesc(reshape(dataBasis(:,1:3)*mdl.Beta(1:3),[100 100]))
% 
% figure; 
% subplot(121);
% % imagesc(reshape(dataRaw(:,12),[100 100]))
% imagesc(reshape(mean(dataRaw(:,1:500),2),[100 100]))
% subplot(122);
% % imagesc(reshape(dataRaw(:,12+500),[100 100]))
% imagesc(reshape(mean(dataRaw(:,500+[1:500]),2),[100 100]))
% % 
% figure; imagesc(reshape(mean(dataRaw(:,1:500),2),[100 100])-reshape(mean(dataRaw(:,500+[1:500]),2),[100 100]));
% figure; imagesc(reshape(mean(dataRaw(:,1:500),2),[100 100]),[100 100]);
% figure; imagesc(reshape(mean(dataRaw(:,1:500),2),[length(centerValsR) length(centerValsC)]));
% figure; imagesc(reshape(mean(dataRaw(:,1:500),2),[length(centerValsR) length(centerValsC)])-reshape(mean(dataRaw(:,500+[1:500]),2),[length(centerValsR) length(centerValsC)]));
% 
% save([isetbioRootPath '/local/' sprintf('landoltC_noise_%d_res.mat',noiseVal)],'P');
% 
% 
% 
% figure; 
% plot(180-degRotArr,mean(P'),'x'); 
% set(gca,'xscale','log'); 
% grid on


%

% Fit psychometric curve to thresholds as a function of contrast
%     [xData, yData] = prepareCurveData( maxContrast-.0, rocArea);
    
    % Set up fittype and options.
    ft = fittype( '1 - 0.5*exp(-(x/a)^b + c)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.323369521886293 0.976303691832645 0];
    opts.lower = [1 -1 0]; opts.upper = [180 50 0];
    
    % Fit a curve between contrast level (x) and probability of correction
    % detection.
    % bootWeibullFit(stimLevels, nCorrect, nTrials, varargin)
    % in computationaleyebrain/simulations/Pixel Visibility/ ...
    
%     xData = (180-degRotArr)'*ones(1,size(P,2)); 
    xData = (pixelRadArr)'*ones(1,size(P,2)); 
    yData = P/100;
%     figure; scatter(xData(:),yData(:))
    
%     [fitresult, gof] = fit( xData(:), yData(:), ft, opts );
    
%     [fitresult, gof] = fit( mean(xData'),mean(yData'), ft, opts );
    %%
    if plotSingle == 1
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
%     h = plot( fitresult, mean(xData'),mean(yData') );
    hold on; scatter(xData(:),yData(:),40,'filled');
%     hold on; scatter(mean(xData'),mean(yData'),40,'filled');
    set(gca,'xscale','log')
    legend( h, 'data', 'fitted curve', 'Location', 'NorthWest' );
    % Label axes
    xlabel('Angle (degrees)');
    ylabel('p(Correct)');
    grid on
    thresh1 = fitresult.a;
    title(sprintf('Detection, \\alpha = %1.2f pixels',(thresh1)));
    set(gca,'fontsize',16)
%     axis([0 1 0.5 1]);
% axis([0.8620  243.8354    0.5454    0.9783])
%     end
    %%
    
    contrastData = freqArr'*ones(1,size(P,2)); 
    % Plot fit with data.
%     figure( 'Name', 'untitled fit 1' ); %hold on;
%     h = plot( fitresult, xData, yData );

%     pData = plot(fitresult); pDataCopy{freqInd} = pData;
%     plot3(pData.XData,freqArr(freqInd)*ones(length(pData.XData),1),pData.YData,'linewidth',4); hold on;
% %     scatter3(xData(:), freqArr(freqInd)*ones(length(xData(:)),1),yData(:));
% %     h = plot( fitresult, mean(xData'),mean(yData') );
%     hold on; scatter3(xData(:),freqArr(freqInd)*ones(length(P(:)),1), yData(:),40,'filled');
    
%     hold on; plot3(mean(xData'),freqArr(freqInd)*ones(length(mean(yData')),1),mean(yData'),'-o','linewidth',6);%,40,'filled');
    hold on; plot(mean(xData'),mean(yData'),'-o','linewidth',6);%,40,'filled');
%     hold on; scatter(mean(xData'),mean(yData'),40,'filled');
    set(gca,'xscale','log')
%     legend( h, 'data', 'fitted curve', 'Location', 'NorthWest' );
    % Label axes
    xlabel('Radius (pixels)');
    ylabel('Fraction Correct');
    zlabel('Fraction Correct');
    grid on
    thresh1 = fitresult.a;
    title(sprintf('Detection, \\alpha = %1.2f pixels',(thresh1)));
    set(gca,'fontsize',16)
%     axis([0 1 0.5 1]);
    % hold off;
    drawnow;
end
%     print(gcf,'-dpng',['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/landoltC_pixelRadius.png']);
end
end
% axis([1 100 0.001 .08 0 1])

axis([1 100 0 1])
clear filterMat2 spikeTrials trialReconPlay2 trialReconPlayXW


%%
% % Pbig(1,:,1:7)=50;
% % Pbig(1,:,23:25)=100;

% load('ws_landolt_healthy_aug17_2018_4reps_out0.mat');
% 
freqArr2(1,1,:) = freqArr;

% % % added from sep17add
% contAdd = ones(1,1,6);
% % % gratings prima
% % % contAdd(1,:,:) = [55.9375   63.4375   62.1875   69.6875   75.0000   89.6875];
% gratings only on
% contAdd(1,:,:) = [57.8125   55.3125   60.3125   63.4375   67.5000   82.1875];
% gratings on and off
% contAdd(1,:,:) = [57.1875   54.3750   53.1250   50.6250   50.0000   50.6250];
% Pbig(1,:,17:22)=repmat(contAdd,[1 4 1]);
% contrastArr(17:22) = [.002:.001:.007];

% xa = (repmat(freqArr2,[length(pixelRadArr) length(freqArr) 1]));
xa = (repmat(freqArr2(1:end),[length(pixelRadArr) numReps 1]));
% pixelArr2(1,1,:) = pixelRadArr;
ya = permute(repmat(pixelRadArr,[numReps 1 length(contrastArr)]),[2 1 3]);
% figure; scatter3(xa(:),ya(:),Pbig(:));
% set(gca,'xScale','log');

radVal=1;
[~,sortInd] = sort(freqArr(1:end),'ascend');
freqArrSort = freqArr((sortInd));

% freqArrSort(1) = 2e-5;
if size(Pbig,2)>1
Pbigsort = squeeze(mean(Pbig(radVal,:,:)));
else

Pbigsort = squeeze((Pbig(radVal,:,:)));
end    
% figure; 
hold on;
for radVal = 1%1:15
    plot(freqArrSort(:),Pbigsort(sortInd),'-ro','linewidth',4);
end

xlabel('Contrast'); ylabel('Fraction Correct');
% title('
grid on;
set(gca,'xScale','log');
set(gca,'fontsize',14)
% axis([3.1676e-07 .2 48 100]);
% axis([2.8515e-05 .02 48 100]);
legend('Healthy','Prosthesis','location','se');
title('Gratings Orientation Discrimination');