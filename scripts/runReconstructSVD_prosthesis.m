


%%run Linear reconstruction for on, off, and joint on/off training 
windowsize = 30;
% disp('Loading ON parasol spike responses...')
% matfON = matfile('../dat/spikeResp_may26_on.mat');
% disp('Loading OFF parasol spike responses')
% matfOFF = matfile('../dat/spikeResp_may26_off.mat');
% %load ON parasols
% % disp('Loading stimulus movie...')
% % load ../dat/movie_may26.mat

matfON = matfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/spikeResp_prosthesis_stixds4');
% matfON = matfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/spikeResp_prosthesis_on_parasol_single2');
% matfOFF = matfile('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/spikeResp_prosthesis_off_parasol');

disp('Loading stimulus movie...')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/movie_prosthesis_on_parasol')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/movie_prosthesis_off_parasol')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/movie_prosthesis_on_parasol_single2')
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/movie_prosthesis_stixds4')


movielength = 240000;%size(stim,2);
% numOnCells = size(matfON.spikeResp,1);
% numOffCells = size(matfOFF.spikeResp,1);


disp(['Total Movie Length in Frames: ' num2str(movielength)]);
% disp(['Number of ON Cells: ' num2str(numOnCells)])      
% disp(['Number of OFF Cells: ' num2str(numOffCells)])



% disp('Training ON only filters...')

tic
% fileext = 'may26_on';
% [recons_test, filters]=linearReconstructSVD(stim,matfON.spikeResp,fileext, windowsize);

% fileext = 'may26_off';
% [recons_test, filters]=linearReconstructSVD(stim,matfOFF.spikeResp,fileext, windowsize);

% disp('Training ON/OFF filters...')
fileext = 'on_parasol_stixds4';
spikeResp1 = vertcat(matfON.spikeResp);%, matfOFF.spikeResp);
linearReconstructSVD0_pros(stim,spikeResp1,fileext, windowsize,1000,1);
