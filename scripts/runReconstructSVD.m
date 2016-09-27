


%%run Linear reconstruction for on, off, and joint on/off training 
windowsize = 30;
disp('Loading ON parasol spike responses...')
matfON = matfile('../dat/spikeResp_may26_on.mat');
disp('Loading OFF parasol spike responses')
matfOFF = matfile('../dat/spikeResp_may26_off.mat');
%load ON parasols
% disp('Loading stimulus movie...')
% load ../dat/movie_may26.mat

movielength = 360000;%size(stim,2);
% numOnCells = size(matfON.spikeResp,1);
numOffCells = size(matfOFF.spikeResp,1);


disp(['Total Movie Length in Frames: ' num2str(movielength)]);
% disp(['Number of ON Cells: ' num2str(numOnCells)])      
disp(['Number of OFF Cells: ' num2str(numOffCells)])



% disp('Training ON only filters...')

tic
% fileext = 'may26_on';
% [recons_test, filters]=linearReconstructSVD(stim,matfON.spikeResp,fileext, windowsize);

% fileext = 'may26_off';
% [recons_test, filters]=linearReconstructSVD(stim,matfOFF.spikeResp,fileext, windowsize);

% disp('Training ON/OFF filters...')
fileext = 'may26_on_off';
spikeResp1 = vertcat(matfON.spikeResp, matfOFF.spikeResp);
linearReconstructSVD0(1,spikeResp1,fileext, windowsize);
