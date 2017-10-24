function obj = buildPrimaGratings(obj, varargin)
%BUILD - builds training set for the recon object
% Run natural scenes through the RGC array for the big four RGC cell types.
%
% inputs:
%   mosaicFile - a string that is used to save the mosaic file
%   buildFile - a string that is used to store the spikes and the movie stim
%
% See also: trainAndTest.m
%
p = inputParser;
p.addParameter('mosaicFile',[],@ischar);
p.addParameter('buildFile',[],@ischar);
p.addParameter('stimFile',[],@ischar);
p.addParameter('respFile',[],@ischar);
p.addParameter('blockIn',1,@isnumeric);
p.addParameter('stimTypeBuild','ns',@ischar);
p.addParameter('nTrials',100,@isnumeric);
p.addParameter('contrast',0.1,@isnumeric);

p.addParameter('pixelWidth',70,@isnumeric);
p.addParameter('currentDecay',2,@isnumeric);

p.addParameter('gratingSpFreq',25,@isnumeric);
p.addParameter('horizontalFlag',false,@isnumeric);
p.addParameter('folderNameTrain','',@ischar);
p.addParameter('folderNameTest','',@ischar);
p.addParameter('filterFile','',@ischar);
p.addParameter('windowSize','',@isnumeric);
p.addParameter('percentSV','',@isnumeric);
p.addParameter('shifttime','',@isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
mosaicFile = p.Results.mosaicFile;
buildFile = p.Results.buildFile;
stimFile = p.Results.stimFile;
respFile = p.Results.respFile;
blockIn = p.Results.blockIn;
stimTypeBuild = p.Results.stimTypeBuild;
nTrials = p.Results.nTrials;
contrast = p.Results.contrast;
horizontalFlag = p.Results.horizontalFlag;

gratingSpFreq = p.Results.gratingSpFreq;

pixelWidth = p.Results.pixelWidth;
currentDecay = p.Results.currentDecay;
folderNameTrain = p.Results.folderNameTrain;
folderNameTest = p.Results.folderNameTest;
filterFile = p.Results.filterFile;
windowSize = p.Results.windowSize;
percentSV = p.Results.percentSV;
shifttime  = p.Results.shifttime;
if isempty(buildFile)
    buildFile = ['primaLandolt_training_' num2str(round(cputime*100))];
end

tic

%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 1.8; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.7;% 3.2;

% Stimulus length = nSteps*nBlocks;
nPixels = 100;
nSteps = 16;%000;
% nBlocks = 15;%30;


%%

numCols = 100; numRows = 100; % gratingSpFreq = 50;
% gratingsMovieRow = sin((2*pi/(numCols/gratingSpFreq))*[1:numCols]);
gratingsMovieRow = 127+contrast*127*sin((2*pi/(numCols/gratingSpFreq))*[1:numCols] + pi/2);
if horizontalFlag
    gratingsMovieIm= (ones(numRows,1)*gratingsMovieRow)';%*tempFreq;
else
    
    gratingsMovieIm= (ones(numRows,1)*gratingsMovieRow);%*tempFreq;
end
    
gratingsMovie = repmat(gratingsMovieIm,[1 1 nSteps]);

gratingsMovieContrast = (gratingsMovie(1:100,1:100,:));
gratingsMovieContrast(:,:,1) = 127*zeros(100,100,1);
% gratingsMovieContrast = 128+127*gratingsMovieContrast;


% Need factor of two on contrast for accurate match
% gratingsMovieContrast = (2*contrast)*(.5*gratingsMovie(1:100,1:100,:));
% gratingsMovieContrast(end,end,:) = .5; gratingsMovieContrast(1,end,:) = -.5;
% gratingsMovieContrast(:,:,1) = zeros(100,100);
% gratingsMovieContrast = 128+127*gratingsMovieContrast;
%% Reconstruction filter parameters

%% Reconstruction filter parameters

% mosaicFile = 'mosaic0';
% % windowSize = 1;
% % percentSV = .05;%.5;%.5;%.25;%.12;
% % % shifttime = 2;
% % shifttime = 4;%15;%3;%15;
% dropout = 0;
% 
% filterFile  = fullfile(folderNameTrain,...    
%     ['filters' mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_dr%d',100*dropout)]);

load(filterFile);

% lambda = .01;
% filterMat2 = zeroFilter(filterMat,lambda);
% clear filterMat;

%% Load image

primaParams.pixelWidth = pixelWidth*1e-6; % meters
primaParams.ecc = 1.8;       % mm
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance


primaRecon = primaArray(gratingsMovieContrast,primaParams);

primaRecon.compute(gratingsMovieContrast);

%% Set signal normalization


% %     s1 = primaRecon.bpMosaic.mosaic{m1}.responseCenter-primaRecon.bpMosaic.mosaic{m1}.responseSurround;
% %     v2(m1,:) = [max(s1(:)) mean(s1(:)) min(s1(:))];
% 
% maxBp = 5.2*contrast+3.12;
% for m1 = 1:2
%      st1 = primaRecon.bpMosaic.mosaic{m1}.responseCenter;  
%      st1scale = maxBp*st1./max(st1(:));     
%      primaRecon.bpMosaic.mosaic{m1}.set('responseCenter',st1scale);
% end
% 
% 
% maxBp = 4*contrast+2.93;
% for m1 = 3:4
%      st1 = primaRecon.bpMosaic.mosaic{m1}.responseCenter;  
%      st1scale = maxBp*st1./max(st1(:));     
%      primaRecon.bpMosaic.mosaic{m1}.set('responseCenter',st1scale);
% end

%%
    
for iTrial = 1:nTrials
    iTrial
    tic
    % primaRecon = primaArray(landoltMovieContrast,primaParams);
    
    primaRecon.computeRGC();
    innerRetina = primaRecon.innerRetina;
    
    %% Get spikes and make reconstruction
    
    spikeResp = mosaicSpikes(primaRecon.innerRetina);
    % save('spikeResp_primaLandolt.mat','spikeResp');
    
    % load(filterFile);
    spikeAug(1,:) = ones(1,size(spikeResp,2));
    spikeAug(2:9717,:) = spikeResp;
%     movRecon = filterMat'*spikeAug;
%     toc
%     movTrials(:,iTrial) = -movRecon(:,6);
    
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     figure; ieMovie(movReconPlay);
    % save('hallwayReconMovie.mat','movRecon');
    
    movRecon = filterMat'*spikeAug;
    toc
    movTrials(:,iTrial) = movRecon(:,8);
    spikeTrials(:,iTrial) = spikeAug(:,8);
    clear spikeResp spikeAug


    trialReconPlay = reshape(movTrials,[100 100 size(movTrials,2)]);
    
    if mod(iTrial,100)==0
        % save(fullfile(reconstructionRootPath,'dat','july17gratings','primaH',[sprintf('primaRecon_%4d.mat',10000*contrast)]),'trialReconPlay');
        
        if horizontalFlag
            save(fullfile(reconstructionRootPath,'dat',folderNameTest,'cont',[sprintf('primaRecon_freqH%4d_cont%6d.mat',100*gratingSpFreq,1e5*contrast)]),'trialReconPlay','contrast','spikeTrials');
            
        else
            save(fullfile(reconstructionRootPath,'dat',folderNameTest,'cont',[sprintf('primaRecon_freq%4d_cont%6d.mat',100*gratingSpFreq,1e5*contrast)]),'trialReconPlay','contrast','spikeTrials');
            
        end
        
    end

end
%% Find best frame
gratingsMovieImage = gratingsMovieContrast(:,:,2);
movRecon0 = movRecon - ones(10000,1)*mean(movRecon);
stdRecon0 = std(movRecon0);
std0 = std(single(gratingsMovieImage(:)));
movReconNorm = std0*movRecon0./(ones(10000,1)*stdRecon0);
% movRecon0 = movRecon - ones(10000,1)*me an(movRecon);

errMov = (movReconNorm - (single(gratingsMovieImage(:))*ones(1,16)));

% errMov = (movReconNorm - RGB2XWFormat(gratingsMovieContrast));
size(errMov)
mse1 = mean(errMov.^2);
figure; plot(mse1)
[mv,mi] = min(mse1)

%%
% cont = .5;
%     5.7294    0.8850   -0.3851
%     0.3851   -0.8850   -5.7294
%     4.9324    0.9173   -0.3830
%     0.3753   -0.9244   -4.9345

% cont = .1
%     3.6468    0.6627   -0.2765
%     0.2765   -0.6627   -3.6468
%     3.3315    0.6994   -0.3237
%     0.3237   -0.7054   -3.3315
    
% cont = .01
%     3.1733    0.6572   -0.3110
%     0.3110   -0.6572   -3.1733
%     2.9753    0.6941   -0.3521
%     0.3278   -0.7001   -2.9753

% cont = 0
%     3.1232    0.6574   -0.3228
%     0.3228   -0.6574   -3.1232
%     2.9329    0.6943   -0.3464
%     0.3375   -0.7003   -2.9329
% 
% figure; scatter([0 .01 .1 .5],[3.12 3.17 3.65 5.7]);
% f1 = 1/((-3.12+[3.12 3.17 3.65 5.72])'\[0 .01 .1 .5]')
% % f1 = 5.2; y = 5.2x+3.12
% % hold on; plot(0:.01:.5,3.12+f1*[0:.01:.5]);
% f2 = 1/((-2.93+[2.93 2.975 3.33 4.93])'\[0 .01 .1 .5]') % = 4.00
% 
% for m1 = 1:4
%     s1 = bpL.mosaic{m1}.responseCenter-bpL.mosaic{m1}.responseSurround;
%     v2(m1,:) = [max(s1(:)) mean(s1(:)) min(s1(:))];
% end
% 
% for m1 = 1:4
%     s1 = primaRecon.bpMosaic.mosaic{m1}.responseCenter-primaRecon.bpMosaic.mosaic{m1}.responseSurround;
%     v2(m1,:) = [max(s1(:)) mean(s1(:)) min(s1(:))];
% end