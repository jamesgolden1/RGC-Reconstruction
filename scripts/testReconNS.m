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

images1 = dir([reconstructionRootPath '/dat/*.tif']);
im1 = (rgb2gray(imread([reconstructionRootPath '/dat/' images1(1).name])));
% im1 = rgb2gray(imread('peppers.png'));
tic
for xblock = 1:4
    xblock
    for yblock = 1:4
        stimulusRGBdata(:,:,21:40) = 2*(double(repmat(im1((xblock-1)*96+[1:96],(yblock-1)*96+[1:96]),[1 1 20])));
        
        os = osCreate('displayrgb');
        os = osSet(os, 'rgbData', stimulusRGBdata);
        % RGC spikes
        % Load mosaic

        load('/Users/james/Downloads/mosaic_all_overlap0.mat')
        % Compute response
        innerRetina.compute(os,'coupling',false);
        % Generate reconstruction
        [movrecons, ~] = irOptimalReconSingle(innerRetina, 0);
        % figure; ieMovie(movrecons);
        movreconsScale = (movrecons-mean(movrecons(:)));
        
        
        imStitchOrig((xblock-1)*96+[1:96],(yblock-1)*96+[1:96],:) =movrecons(:,:,:);
        imStitch((xblock-1)*96+[1:96],(yblock-1)*96+[1:96],:) = movreconsScale(:,:,:);
        imOrig((xblock-1)*96+[1:96],(yblock-1)*96+[1:96],:) = stimulusRGBdata(:,:,:);
        clear innerRetina os stimulusRGBdata
    end
end
toc
%%
figure; 
subplot(121);
imagesc(imOrig(:,:,40)); colormap gray;
subplot(122);
imagesc(1-imStitch(:,:,11)); colormap gray
% subplot(133);
% errIm = ieScale(imStitch(:,:,1))-ieScale(imOrig(:,:,1));
% imagesc(errIm); colormap gray

%%
% clear combMovie
% combMovie(:,1:96,:) = stimulusRGBdata(1:96,1:96,5:end-4);
% combMovie(:,96+[1:96],:) = 1-movreconsScale;
% figure; ieMovie(combMovie);
%%