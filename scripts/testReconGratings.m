


% clear all;
% close all;
disp('Loading testing stimulus data...')
% name = 'bar';
name = 'grating';
%set which testing stimulus to use
if strcmp(name,'grating')
    matfOn = matfile('RGC-Reconstruction/gratings/gratingsON_midget.mat');
    matfOff = matfile('RGC-Reconstruction/gratings/gratingsOFF_midget.mat');
    
    matfOnP = matfile('RGC-Reconstruction/gratings/gratingsON.mat');
    matfOffP = matfile('RGC-Reconstruction/gratings/gratingsOFF.mat');
%     matfOnP = matfile('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OnParasol_36_grating_june10.mat');
%     matfOffP = matfile('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OffParasol_64_grating_june10.mat');

%     matfOn = matfile('C:\Users\James\Documents\GitHub\RGC-ReconAdd\output\gratings\gratingsON_midget.mat')    
%     matfOff = matfile('C:\Users\James\Documents\GitHub\RGC-ReconAdd\output\gratings\gratingsOFF_midget.mat')
%     matfOnP = matfile('../dat/WNstim_response_OnParasol_36_grating_june10.mat');
%     matfOffP = matfile('../dat/WNstim_response_OffParasol_64_grating_june10.mat');
    %load ../output/reconstructedSTA/onParasolGratingSTARecon.mat
    %load ../output/reconstructedSTA/offParasolGratingSTARecon.mat
elseif strcmp(name,'bar')
    matfOff = matfile('../dat/WNstim_response_OffParasol_offBig2_64_barGay_june10.mat');
    matfOn = matfile('../dat/WNstim_response_OnParasol_onBig2_36_barGray_june10.mat');
    %load ../output/reconstructedSTA/onParasolBarSTARecon.mat
    %load ../output/reconstructedSTA/offParasolBarSTARecon.mat
end


blocklength = 152;

numcells = 225+144+64+36;

onSR = matfOn.spikesoutsm;
offSR = matfOff.spikesoutsm;
onPSR = matfOnP.spikesoutsm;
offPSR = matfOffP.spikesoutsm;

spikesout = vertcat(onSR(:,1:15000), offSR(:,1:15000), onPSR(:,1:15000), offPSR(:,1:15000));

spikeRespOn = downSampResp(spikesout, numcells, blocklength);
spikeResp = spikeRespOn;

%%
% icArray = [400:400:4000];
icArray = 1000;%[1000:1000:6000];
for icind = 1:length(icArray)
    
    % load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/filters_may26_parasol_midget_combined_svd_3000_len_100.mat');
    % load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/oct8/filters_may26_fine_all2_svd_3000_len_100.mat')
    % load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/oct8/filters_may26_fine_all2_svd_1500_len_100.mat')
    
    
    rdt = RdtClient('isetbio');
    rdt.crp('/resources/data/reconstruction');
    data = rdt.readArtifact('filters_may26_parasol_midget_combined_svd_3000_len_100', 'type', 'mat');
    filterMat = data.filterMat; clear data;
    
    
    numbins = 8;
    recons_stim_on_off = reconsFromFiltLen(filterMat, spikeResp, numbins);
    
    
    %   % Threshold filter
    %   filterMatInd = find(abs(filterMat)<0.005); filterMat2 = filterMat; filterMat2(filterMatInd)=0;
    %     recons_stim_on_off = reconsFromFiltLen(filterMat2, spikeResp, numbins);
    
    
    % load('C:\Users\James\Documents\GitHub\RGC-ReconAdd\output\svd_reconstruct_shorttrain_midgets\filters_may26_on_long300_shorttime_svd_1000_len_100.mat')
    %     numbins = 12;
    %     recons_stim_on_off = reconsFromFiltLen(filterMat, spikeRespOn, numbins);
    
    %     mov = reshape(stim,96,96,size(stim,2));
    %     movrecons_on_off = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));
    movrecons_on_off_full = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));
    movrecons_on_off =(.25+.5*(movrecons_on_off_full)./median((movrecons_on_off_full(:))));
end
%%
% minmaxv = [min(mov(:)) max(mov(:))];
minmaxorig = [min(movrecons_on_off(:)),max(movrecons_on_off(:))];
minmaxv(1) = minmaxorig(1) + minmaxorig(1)*0.0;
minmaxv(2) = minmaxorig(2) - minmaxorig(1)*0.0;
% minmaxv = [-1 1]; 
figure; hold on;
for it = 1:121
    clf
    imagesc(movrecons_on_off(:,:,it));
%     imagesc(mov(:,:,it));
%     pause(0.05);
   caxis(minmaxv)
   colormap gray
    drawnow
end