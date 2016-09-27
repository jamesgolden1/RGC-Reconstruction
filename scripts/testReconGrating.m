

% clear
% matfOn = matfile('../dat/WNstim_response_OnParasol_36_grating_june10.mat');
% matfOff = matfile('../dat/WNstim_response_OffParasol_64_grating_june10.mat');
    
% filename1 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/Gratingresponse/prosthesis_on_parasol/Gratingstim_response_ProsONParasol4.mat'];

filename1 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/Gratingresponse/prosthesis_on_parasol/Gratingstimlong_response_ProsONParasol3.mat'];
matfOn = matfile(filename1);

filename2 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/Gratingresponse/prosthesis_off_parasol/Gratingstim_response_ProsOFFParasol4.mat'];
% filename2 = ['/Users/james/Documents/MATLAB/EJLExperimentalRGC/Gratingresponse/prosthesis_off_parasol/Barstim_response_ProsOFFParasol1.mat'];
matfOff = matfile(filename2);
% blocklength = 1000;

% In-sample
% frstart = 1;
% frend = frstart + blocklength-1;

% frstart = 336001;
% frend = frstart + blocklength-1;

% disp('Reconstructing with off parasols...')
%test with off parasols
% filtMat_off = matfile('../output/filters_may26_off.mat');
%

%%
%create response matrix for test stimulus

load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_off_parasol_prosthesis_svd_1021.mat')

numcells=34;

blocklength = 152;
frstart = 1;
frend = frstart + blocklength-1;
spikesout = double(matfOff.spikesoutsm);
stim = double(reshape(matfOff.whiteNoiseSmall,40*80,blocklength));
spikeRespOff = downSampRespPhys(spikesout, numcells, blocklength);
recons_stim_off = reconsFromFilt(filterMat, spikeRespOff);
  mov_off = reshape(recons_stim_off,[80 40 122]);
  figure; ieMovie(mov_off);
%%
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_prosthesis_svd_600.mat')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_prosthesis_single_transpose2_svd_1000.mat')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_short3_svd_1000.mat');
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_short4_svd_400.mat');
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_on_parasol_short6_svd_800.mat');
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_may26_on_svd_200.mat')

numcells= 48; 

% blocklength = 152;
blocklength = 352; 
spikesout = double(matfOn.spikesoutsm);
spikeRespOn = downSampRespPhys(spikesout, numcells, blocklength);
% spikeRespOn = downSampResp(spikesout, numcells, blocklength);
% recons_stim_on = reconsFromFilt(filtMat_on.filterMat, spikeRespOn);

recons_stim_on = reconsFromFilt(filterMat, spikeRespOn);
% cellNum= 8;
% recons_stim_on = reconsFromFilt(filterMat((cellNum-1)*20+[1:21],:), spikeRespOn((cellNum),:));
% 
% recons_stim_on = zeros(3200,332);
% for cellNum= [1:20];
% recons_stim_on = recons_stim_on + reconsFromFilt(filterMat((cellNum-1)*20+[1:21],:), spikeRespOn((cellNum),:));
% end

  mov_on = permute(reshape(recons_stim_on,[40 80 blocklength-20]),[1 2 3]);
  figure; 
  ieMovie(mov_on);