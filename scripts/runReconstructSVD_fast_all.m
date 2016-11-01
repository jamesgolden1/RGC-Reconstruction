


%%run Linear reconstruction for on, off, and joint on/off training 
windowsize = 8;
disp('Loading ON parasol spike responses...')
% matfON = matfile('../dat/spikeResp_onMidget.mat');
disp('Loading OFF parasol spike responses')
% matfOFF = matfile('../dat/spikeResp_onMidget_long300.mat');

% matfON = matfile('../dat/spikeResp_onMidget.mat');

% matfON = matfile('C:\Users\James\Documents\matlab\github\RGC-Reconstruction\dat\spikeResp_onParasol_fast');

matfON = matfile('C:\Users\James\Documents\matlab\github\RGC-Reconstruction\dat\spikeResp_all0');
% matfOFF = matfile('../dat/spikeResp_offMidget.mat');
% 
% matfOFFP = matfile('../dat/spikeResp_offParasol.mat');
% 
% matfONP = matfile('../dat/spikeResp_onParasol.mat');
%load ON parasols
disp('Loading stimulus movie...')
% load ../dat/movie_may26.mat

movielength = 1*240000;%size(stim,2);
% numOnCells = size(matfON.spikeResp,1);
% numOffCells = size(matfOFF.spikeResp,1);


disp(['Total Movie Length in Frames: ' num2str(movielength)]);
% disp(['Number of ON Cells: ' num2str(numOnCells)])      
% disp(['Number of OFF Cells: ' num2str(numOffCells)])



% disp('Training ON only filters...')

% tic
% fileext = 'may26_on';
% trainSizeArray = .6;%[.6/8:.6/8:.6]; 
% trainInd = 1;
% includedComponentsArray = [ 1200:400:3600];
% for incInd = 1:length(includedComponentsArray)
% % linearReconstructSVD_short_midgets(1,matfON.spikeResp,fileext, windowsize,includedComponentsArray(incInd),trainSizeArray(trainInd));
% end

% clear matfON
% % 
% matfOFF = matfile('../dat/spikeResp_may26_off.mat');
% numOffCells = size(matfOFF.spikeResp,1);
% fileext = 'may26_off_short';
% 
% % trainSizeArray = [.6/8:.6/8:.6];


% fileext = 'may26_on_midget_both_neg';
fileext = 'onParasol_fast';
trainSizeArray = 1;%[.6/8:.6/8:.6]; 
trainInd = 1;
includedComponentsArray = 3952;%3000;%[1000 4000 5000 6000];% [ 1200:400:3600];

srON = matfON.spikeResp;
% srOFF = matfOFF.spikeResp;
% clear matfON matfOFF
spikeResp1 = vertcat(srON(:,1:240000));%, srOFF(:,1:240000), matfONP.spikeResp, matfOFFP.spikeResp);

% spikeRespSingle = single(matfOFF.spikeResp);

% spikeResp1 = vertcat(matfON.spikeResp, matfOFF.spikeResp);

% spikeResp1 = matfOFFP.spikeResp;

% spikeResp1 = matfONP.spikeResp;

clear matfOFF

for incInd = 1%:length(includedComponentsArray)
% linearReconstructSVD_short_midgets(1,matfOFF.spikeResp,fileext, windowsize,includedComponentsArray(incInd),trainSizeArray(trainInd));
linearReconstructSVD_short_midgets_both(1,spikeResp1,fileext, windowsize,includedComponentsArray(incInd),trainSizeArray(trainInd));
end

% includedComponentsArray = [3000];
% linearReconstructSVD(1,matfOFF.spikeResp,fileext, windowsize,includedComponentsArray);

% % for trainInd = 1:length(trainSizeArray)
% % linearReconstructSVD_short(1,matfOFF.spikeResp,fileext, windowsize,includedComponentsArray,trainSizeArray(trainInd));
% % end
% % includedComponentsArray = [400];
% % for trainInd = 1:length(trainSizeArray)
% % linearReconstructSVD_short(1,matfOFF.spikeResp,fileext, windowsize,includedComponentsArray,trainSizeArray(trainInd));
% % end

% % disp('Training ON/OFF filters...')
% fileext = 'may26_on_off_short';
% spikeResp1 = vertcat(matfON.spikeResp, matfOFF.spikeResp);
% % 
% % includedComponentsArray = [1300:100:1500]
% % linearReconstructSVD(1,spikeResp1,fileext, windowsize,includedComponentsArray);
% 
% trainSizeArray = [.6/8:.6/8:.6];
% includedComponentsArray = [3000];
% % linearReconstructSVD(1,matfOFF.spikeResp,fileext, windowsize,includedComponentsArray);
% for trainInd = 1:length(trainSizeArray)
% linearReconstructSVD_short(1,spikeResp1,fileext, windowsize,includedComponentsArray,trainSizeArray(trainInd));
% end
% 
% includedComponentsArray = [400];
% for trainInd = 1:length(trainSizeArray)
% linearReconstructSVD_short(1,spikeResp1,fileext, windowsize,includedComponentsArray,trainSizeArray(trainInd));
% end
