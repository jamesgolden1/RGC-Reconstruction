function testReconNS(varargin)
% 
% Tests the accuracy of a movie stimulus reconstructed from RGC spikes. 
% 
% An outer segment osDisplayRGB object is created in order to provide RGB
% input to an RGC mosaic. The rgbData is set to the movie stimulus, and the
% RGC spikes in response to the movie are computed. The mosaic spikes are
% passed to the irOptimalReconSingle function, which generates the
% reconstructed movie from the spikes and decoding filters.
% 


p = inputParser;

p.addParameter('mosaicFile',[],@ischar);
p.addParameter('filterFile',[],@ischar);
p.parse(varargin{:});
filterFile = p.Results.filterFile;
mosaicFile = p.Results.mosaicFile;

%% RGB stimulus stored in os

% clear

% Make this cone mosaic
os = osCreate('displayrgb');

if isunix || ismac
    images1 = dir([reconstructionRootPath '/dat/*.tif']);
    im1 = (rgb2gray(imread([reconstructionRootPath '/dat/' images1(1).name])));
else
    images1 = dir([reconstructionRootPath '\dat\FoliageBig\*.tif']);
    im1 = (rgb2gray(imread([reconstructionRootPath '\dat\FoliageBig\' images1(1).name])));
end
% im1 = rgb2gray(imread('peppers.png'));
tic
for xblock = 1:2
    xblock
    for yblock = 1:2
        stimulusRGBdata(:,:,21:40) = 1*(double(repmat(im1((xblock-1)*96+[1:96],(yblock-1)*96+[1:96]),[1 1 20])));
        
        os = osCreate('displayrgb');
        os = osSet(os, 'rgbData', stimulusRGBdata);
        % RGC spikes
        % Load mosaic

        if isunix || ismac
            filenameRGC = [reconstructionRootPath '/dat/' mosaicFile '.mat'];
        else
            filenameRGC = [reconstructionRootPath '\dat\' mosaicFile '.mat'];
        end
        load(filenameRGC);
        % load('/Users/james/Downloads/mosaic_all_overlap0.mat')
        
        % Compute response
        innerRetina.compute(os,'coupling',false);
        % Generate reconstruction
        % [movrecons, ~] = irOptimalReconSingle(innerRetina, 0);
        pRecons.innerRetina = innerRetina; pRecons.filterFile = filterFile; pRecons.numbins =1;
%         [movrecons, ~] = irOptimalReconSingle(pRecons);
        [movrecons, ~] = irOptimalReconSingle(pRecons);
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
imagesc(1-imStitch(:,:,33)); colormap gray
ph=1;
% subplot(133);
% errIm = ieScale(imStitch(:,:,1))-ieScale(imOrig(:,:,1));
% imagesc(errIm); colormap gray

%%
% clear combMovie
% combMovie(:,1:96,:) = stimulusRGBdata(1:96,1:96,5:end-4);
% combMovie(:,96+[1:96],:) = 1-movreconsScale;
% figure; ieMovie(combMovie);
%%