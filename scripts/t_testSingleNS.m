% t_pixiumTestNSSingle
% 
% Reconstructs the hallway navigation movie stimulus with decoder trained
% on the prosthesis simulation.
%
% Stimulating the activations of retinal ganglion cells with an array of
% electrodes. The stimulus image is defined, the electrode array is
% generated, the electrode activations are computed, the RGC mosaics are
% generated, the RGC mosaic responses are computed and the stimulus is
% inferred from the RGC mosaic resposnes using a linear decoder.
%
% Outline of computation:
% 1. Load image/movie
% 2. Outer segment representation
% 3. Build electrode array
% 4. Compute electrode activations from image/movie
% 5. Build RGC array
% 6. Calculate RGC input - only one, change to input from multiple
% 7. Build RGC activation functions
% 8. Compute RGC activations/spikes
% 9. Invert representation to form image/movie
% 10. Tile over big FOV
%
% 12/2016 JRG (c) isetbio
%
%  3bd0152 % nov 19\

%% Initialize
clear;
% ieInit;
tic

%% Parameters to alter
clear electrodeArray

% Electrode size
% Set the size of implant pixels
electrodeArray.width = 30e-6; % meters
% electrodeArray.width = 140e-6; % meters

% Retinal patch eccentricity
patchEccentricity = 4; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2/3;

% % % % % % KNOBS
pulseFreq = 25;           % Hz, electrode pulse frequency
pulseDutyCycle = .2;       % Fraction of cycle pulse is on
irradianceFraction = .5;  % Fraction of maximum irradiance 

% Stimulus length
nSteps = 520;

percentDead = 0;

% %%% Grating subunit stimulus
%
% % params.barWidth = bwL;
% % iStim = ieStimulusGratingSubunit(params);
%
% % iStim = iStimC;
% absorptions = iStim.absorptions;
% % movingBar = iStim;
%
% nSteps = size(iStim.sceneRGB,3);
% params.nSteps = nSteps;
%
% size(iStim.sceneRGB)

params.nSteps = nSteps;
params.row = 96;
params.col = 96;
params.fov = 1.6;

%% Loop for different resizing factors
% The hall movies captures a large field of view. In order to reconstruct
% this large-FOV stimulus with a small RGC mosaic, we must tile the mosaic
% over the stimulus a number of times. The mosaic is necessarily small
% because the training algorithm takes a long time.

for rsFactor = 6%[1 2 3 5 6]
    
    %% Resize the hallway movie stimulus for tiling
    rsFactor
    tic
%     load('C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\stimuli\hallMovie.mat')
%     % load([reconstructionRootPath '\dat\stimuli\hallMovie.mat'])
%     szFrames = size(vidFrame,3);
%     hallMovieResize = zeros(rsFactor*96,rsFactor*96,szFrames);
%     for ii = 1:100%szFrames
%         hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*96,rsFactor*96]);
%     end
%     
%     % Set hallway movie stimulus
%     testmovieshort = (255*ieScale(hallMovieResize)); clear hallMovieResize;
    
    images1 = dir([reconstructionRootPath '\dat\FoliageBig\*.tif']);
    
    %     while pctr <= 12350
    im1 = imread([reconstructionRootPath '\dat\FoliageBig\' images1(2).name]);
    %     imrs = reshape(rgb2gray(im1),size(im0,1),size(im0,2));
    imrs = rgb2gray(im1);
    figure; imagesc(imrs); colormap gray;
    testmovieshort = 128*zeros([size(imrs) 20]);
    
    for k = 1:30
        testmovieshort(:,:,k) = (imrs);
    end
    
    % Stimulus parameters
    paramsStim.nsteps = 1;%nFrames;%size(testmovieshort,3);
    paramsStim.timeInterval = 1/125;%0.001; % secclose all
    
    paramsStim.expTime = 1/125;%0.001; % sec
    nFrames = nSteps;
    
    paramsStim.nsteps = nFrames;
    paramsStim.fov = 8;
    paramsStim.radius = 36/3*1e-6;
    paramsStim.theta = 330;
    paramsStim.side = 'left';
    paramsStim.fov = fov;
    blockctr = 0;
    
    tic
    
    %% Loop over each tiled mosaic
    for iblock = 1:rsFactor
        for jblock = 1:rsFactor
            [iblock jblock]
            % Get iStim structure with os object with desired properties
            blockctr = blockctr+1;
            paramsStim.nsteps = 1;
            iStim = ieStimulusMovie(testmovieshort((iblock-1)*96+[1:96],(jblock-1)*96+[1:96],1:10),paramsStim);
                        
            iStim.absorptions = iStim.sensor;
            clear movingBar
            movingBar = iStim;
            
            %% Outer segment calculation
            % There is no simulated outer segment, this identity outer segment acts as
            % a pass-through holder of the stimulus intensity information.
            
            % Input = RGB
            os = osCreate('displayrgb');
            
            % Get retinal patch properties for electrode array
            sceneSize = sceneGet(movingBar.scene,'size');
            retinalPatchWidth = sensorGet(movingBar.absorptions,'width','m');
            % retinalPatchWidth = sceneGet(movingBar.scene,'width');
            
            % retinalPatchHeight = sensorGet(movingBar.absorptions,'height','m');
            retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;
            
            % % % coneSpacing = scene.wAngular*300
            % coneSpacing = sensorGet(sensor,'dimension','um');
            os = osSet(os, 'patchSize', retinalPatchWidth);
            
            timeStep = sensorGet(movingBar.absorptions,'time interval','sec');
            os = osSet(os, 'timeStep', timeStep);
            
            movingBar.sceneRGB = testmovieshort((iblock-1)*96+[1:96],(jblock-1)*96+[1:96],:);
            % movingBar.sceneRGB = (contrastElectrode)*(movingBar.sceneRGB - 0.5)+0.5;
            os = osSet(os, 'rgbData', movingBar.sceneRGB);
            
            sceneRGB_Healthy = (1)*(movingBar.sceneRGB - 0.5)+0.5;
            osHealthy = os;
            osHealthy = osSet(osHealthy, 'rgbData', sceneRGB_Healthy);                                  
            
            %% Build RGC array for healthy retina            
            
            %% Build RGC array for healthy retina
            
            clear paramsIR innerRetinaHealthy
            filenameRGC = ['C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_RGC.mat'];
            load(filenameRGC);
            innerRetinaHealthy3 = innerRetina;
    
            % filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OnMidget_RGC.mat'];
            filenameRGC = ['C:\Users\James\Documents\GitHub\onMidget1\WNstim_response_OnMidget_RGC.mat'];
            load(filenameRGC);
            innerRetinaHealthy4 = innerRetina;
            
            filenameRGC = ['C:\Users\James\Documents\GitHub\May26_offBig2\WNstim_response_OffParasol_RGC.mat'];
            load(filenameRGC);
            innerRetinaHealthy2 = innerRetina; % clear innerRetina;
            
            filenameRGC = ['C:\Users\James\Documents\GitHub\may26_onBig2\WNstim_response_OnParasol_RGC.mat'];
            load(filenameRGC);
            innerRetinaHealthy = innerRetina;

            innerRetinaHealthy = irCompute(innerRetinaHealthy,osHealthy,'coupling',false);
            innerRetinaHealthy2 = irCompute(innerRetinaHealthy2,osHealthy,'coupling',false);
            innerRetinaHealthy3 = irCompute(innerRetinaHealthy3,osHealthy,'coupling',false);
            innerRetinaHealthy4 = irCompute(innerRetinaHealthy4,osHealthy,'coupling',false);
            %% Do optimal reconstruction
            [movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalRecon(innerRetinaHealthy, innerRetinaHealthy2, innerRetinaHealthy3, innerRetinaHealthy4, percentDead);
            
            %% Save for tiling       
            
            if ismac || isunix
                save([phospheneRootPath '/dat/wnStatic_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
            else
                save([phospheneRootPath '\dat\wnStatic_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
            end
            toc
        end
    end
    
    %% Tile reconstructed movies 
        blockctr = 0; percentDead = 0; numbins = 8;

    movieRecon = zeros(rsFactor*96,rsFactor*96,size(movrecons_on_offHealthy{1},3));
%     movieRecon = zeros(rsFactor*96,rsFactor*96,599);
    tic
    for iblock = 1:rsFactor
        for jblock = 1:rsFactor
            blockctr = blockctr+1;
            if ismac || isunix
                load([phospheneRootPath '/dat/wnStatic_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
            else
                load([phospheneRootPath '\dat\wnStatic_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
            end
%             clear movrecons_on_offHealthy            
%             pOpt.innerRetina = innerRetina;
%             [movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalReconSingle(pOpt);
                        
            movieTmp =  movrecons_on_offHealthy{1};
            movieRecon((iblock-1)*96+[1:96],(jblock-1)*96+[1:96],1:size(movieTmp,3)) = movieTmp;% 255*ieScale(movieTmp - mean(movieTmp(:)));
            clear movrecons_on_offHealthy
        end
    end
     
    figure; ieMovie(movieRecon);
    
    %% Save as output with original movie
%     fig=figure;
%     
%     % set(fig,'position',[    624         437        1018         541]);
%     set(fig,'position',[624   704   537   274]);
%     if ismac || isunix
%         aviobj = avifile([phospheneRootPath '/dat/prosthesis_recon_' num2str(rsFactor) '_ns.avi'])
%     else
%         aviobj = avifile([phospheneRootPath '\dat\prosthesis_recon_' num2str(rsFactor) '_ns.avi'])
%     end
%     aviobj.Fps = 30;
%     shiftval = 4;
%     movieComb = 255*irradianceFraction*pulseDutyCycle*ieScale(movieRecon(:,:,1:567-shiftval+1));
%     movieComb(:,rsFactor*96+[1:rsFactor*96],:) = 255*irradianceFraction*pulseDutyCycle*ieScale(testmovieshort(:,:,shiftval+1:567+1));
%     
%     for k=1:size(movieRecon,3)-60
%         % imagesc(movieRecon(:,:,k)); colormap gray;
%         
%         image(movieComb(:,:,k)); colormap gray; axis image
%         caxis([0 255]);
%         F = getframe(fig);
%         aviobj = addframe(aviobj,F);
%     end
%     close(fig)
%     aviobj = close(aviobj);
    
    %%
    % clear movieComb movieRecon testmovieshort vidFrame
    toc
end