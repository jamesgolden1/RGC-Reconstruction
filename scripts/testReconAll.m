% testReconAll
% 
% Tests the accuracy of a movie stimulus reconstructed from RGC spikes. 
% 
% An outer segment osDisplayRGB object is created in order to provide RGB
% input to an RGC mosaic. The rgbData is set to the movie stimulus, and the
% RGC spikes in response to the movie are computed. The mosaic spikes are
% passed to the irOptimalReconSingle function, which generates the
% reconstructed movie from the spikes and decoding filters.
% 

%% RGB stimulus stored in os

clear

% Make this cone mosaic
os = osCreate('displayrgb');

% params.nSteps = 100;
% params.row = 96;
% params.col = 96;
% params.fov = 1.6;
% iStim = ieStimulusBinaryWhiteNoise(params);
% stimulusRGBdata = iStim.sceneRGB(:,:,:,1);    

% iStim = ieStimulusGratingSubunitFast('barWidth',20,'row',96,'col',96,'nSteps',100,'freq',5);
% stimulusRGBdata = iStim.sceneRGB(:,:,:,1);


% load([ reconstructionRootPath  '\dat\movsm_' num2str(1) '.mat'],'movsm');
% natScenes = movsm(1:96,1:96,randperm(nSteps));

% images1 = dir([reconstructionRootPath '/dat/*.tif']);
% im1 = rgb2gray(imread([reconstructionRootPath '/dat/' images1(1).name]));
im1 = rgb2gray(imread('peppers.png'));
stimulusRGBdata = ieScale(double(repmat(im1(200+[1:96],200+[1:96]),[1 1 100])));
os = osSet(os, 'rgbData', stimulusRGBdata);


%% RGC spikes
% Load mosaic
load('/Users/james/Downloads/mosaic_all_overlap0.mat')
% Compute response
innerRetina.compute(os,'coupling',false);
% Generate reconstruction
[movrecons, ~] = irOptimalReconSingle(innerRetina, 0);
figure; ieMovie(movrecons);
movreconsScale = ieScale(movrecons);
% 
% movreconsThresh = zeros(size(movreconsScale));
% movreconsThresh(movreconsScale>mean(movreconsScale(:))) = 1;
% figure; ieMovie(movreconsThresh);
%%
clear combMovie
combMovie(:,1:96,:) = stimulusRGBdata(1:96,1:96,5:end-4);
combMovie(:,96+[1:96],:) = 1-movreconsScale;
figure; ieMovie(combMovie);
%%
for shiftFrame = 1:22
    mseTmp = ((movreconsScale(:,:,shiftFrame:shiftFrame+70) - stimulusRGBdata(:,:,1:71)).^2);
    reconError(shiftFrame) = sqrt(sum(mseTmp(:)));
    clear mseTmp
end

figure; plot(reconError,'-x')