
load('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/may22/spikeResp_hallway.mat');
stim = RGB2XWFormat(loadHallStimulus(560));
reconHealthy = recon();
mse = reconHealthy.testCV(stim,spikeResp(:,1:560));

mse = reconHealthy.testCV(stim,spikeResp(randperm(size(spikeResp,1)),1:560));
%    79.0022
%    79.0095
%    78.5099

mse = reconHealthy.testCV(stim,spikeResp(:,randperm(560)));

%    76.9084
%    76.9483
%    77.1655


%    47.7779 % 70 learning
%    38.7041 % 35
%    40.3098 % 9

onParasolInd = 27*31;
offParasolInd = onParasolInd+31*35;
onMidgetInd = offParasolInd+54*62;
%  load('filters_mosaic0_sv50_w1_sh4_tr80.mat')
spikeNoLearn = spikeResp;%(randperm(size(spikeResp,1)),randperm(size(spikeResp,2)));
spikeNoLearn(onParasolInd+1:offParasolInd,:) = 0;
spikeNoLearn(onMidgetInd+1:end,:) = 0;

% spikeNoLearn(1:onParasolInd,:) = 0;
% spikeNoLearn(offParasolInd+1:onMidgetInd,:) = 0;
mse = reconHealthy.testCV(stim,spikeNoLearn(:,1:560));

%    29.5275 %healthy w border
% 28.1300 % healthy wo border
% 37.4126 %only on
% 33.7504 %only on wo border

% healthy
% 
% 
% ans =
% 
%   34.6454
% 
% 
% ans =
% 
%   36.7515
% 
% 
% ans =
% 
%   43.3490
% 
% [12:11] 
load('/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/june16prima/spikeResp_hallway.mat');
stim = RGB2XWFormat(loadHallStimulus(560));
reconHealthy = recon();
mse = reconHealthy.testCV(stim,spikeResp(:,1:560));
% pros
% 42.825 % w border

%    38.4046 wo border
% ans =
% 
%   52.2066
% 
% 
% ans =
% 
%   51.1401
% 
% 
% ans =
% 
%   50.1939
% mean =    51.1803;
stdStim = std(stim(:));
mse = reconHealthy.testCV(128+stdStim*randn(size(stim)),spikeResp(:,1:560));
%  81.9407
%    82.0037
%    82.0179
meanSpike = mean(spikeResp(:));
stdSpike =  std(spikeResp(:));
% mse = reconHealthy.testCV(stim,2*meanSpike*(rand(size(spikeResp(:,1:560)))));

%    77.1653
%    76.6953
%    76.6321
mse = reconHealthy.testCV(stim,meanSpike+stdSpike*(randn(size(spikeResp(:,1:560)))));
%    77.8432
%    77.4910
%    77.4753

% w pitch decay spikes on pros decoding
%    78.1500
%    77.8925
%    77.8153
   
% w pitch decay spikes on healthy decoding

%    79.4777
%    79.3904
%    79.2331
% mean =   79.3671;
mse = reconHealthy.testCV(128+stdStim*randn(size(stim)),10*(rand(size(spikeResp(:,1:560)))));

%    82.1148
%    82.1477
%    82.1484
mse = reconHealthy.testCV(mean(stim(:))+stdStim*randn(size(stim)),meanSpike+stdSpike*(randn(size(spikeResp(:,1:550)))));
%    82.0435
%    82.0449
%    82.0550

% % shuffled, sv = .25% 
%          0         0   78.3440
%          0         0   78.4472
%          0         0   78.1523
%          0         0   77.7713

% % shuffled, sv = .5% 
%          0         0   78.8935
%          0         0   79.0824
%          0         0   78.7665
%          0         0   78.8289


%% testCV
% both
% 
%    75.3441   81.3350   79.5408   78.0975
%    74.1407   80.4186   79.9132   79.5959
%    75.5987   73.2060   72.8124   75.3700
%    75.6536   75.8101   76.5910   76.4740

% only off

% 
%    88.5405   87.0432   80.4448   79.2160
%    88.0536   81.8464   78.4393   78.7468
%    85.3061   75.6519   69.6691   72.4435
%    84.7291   84.0611   82.1206   81.7846

% only on
%    64.9192
%    63.9589
%    68.3832
%    68.9738
% % 
% ans =
% 
%          0   68.3093   76.1724   76.7640
%          0   75.7647   78.5251   78.2878
%          0   76.1740   81.9469   80.4041
%          0   69.8523   72.3767   72.5750

% on sv .25
% 
%    66.6303   68.8197   75.7738   76.2480
%    66.1709   75.2423   77.9555   77.4959
%    69.6325   76.5788   82.6378   80.9844
%    69.6717   70.4753   72.9874   73.1215


%%

% % no permute
% msePitchDecay = [...
%     79.0931   79.3456   77.5131   77.2088;...
%    79.0884   78.2584   77.9641   77.7447;...
%    79.8809   79.7761   79.3101   78.7424;...
%    79.6584   79.7212   79.1678   79.1792];

%  permute
% msePitchDecay = [...
%    75.3441   81.3350   79.5408   78.0975
%    74.1407   80.4186   79.9132   79.5959;...
%    75.5987   73.2060   72.8124   75.3700;...
%    75.6536   75.8101   76.5910   7629.5275 29.5275.4740];
 
% % w border
% msePitchDecay = [67.6513;...
%    64.9076;...
%    65.3225;...
%    70.1967];

% no border
  msePitchDecay = [ 63.3279;...
   60.1911;...
   61.7649;...
   67.1843];
   
mseShuffle = [77.88; 77.88; 77.88; 77.88];

% w border
% mseOnlyOn = [...
%     
%    51.5337   63.7182   81.1329   83.3590;...
%    44.2077   76.3429   89.2540   89.2598;...
%    38.1882   56.7302   89.0743   87.2737;...
%    40.6945   41.5071   45.0369   45.9518];

% w/o border
mseOnlyOn = [...
   48.4719;...
   41.5069;...
   35.9057;...
   35.5269];
%    
%          0   61.6673
%          0   74.4299
%          0   53.6097
%          0   36.4317
%          
%          0         0   78.8349
%          0         0   88.7946
%          0         0   89.4602
%          0         0   40.9375
% mseOnlyOn = [...    
%          0   68.3093   76.1724   76.7640;...
%          0   75.7647   78.5251   78.2878;...
%          0   76.1740   81.9469   80.4041;...
%          0   69.8523   72.3767   72.5750];
% mseOnlyOn(:,1) = [...
%    64.9192;...
%    63.9589;...
%    68.3832;...
%    68.9738];

pixelPitchArr = [70 70/2 70/4 70/8];% 70/16 70/32];
currentDecayArr = [2 4 8 16];

figure; hold on;
% plot(pixelPitchArr, msePitchDecay(:,1)/255,'-o','linewidth',4,'markersize',5,'color','k'); 

line([0 70],[   76.9454   76.9454 ]/255,'linestyle','--','linewidth',3,'color','r');
plot(pixelPitchArr, msePitchDecay/255,'-o','linewidth',4,'markersize',5,'color','k'); 
% plot(pixelPitchArr, mseShuffle/255,'-o','linewidth',4,'markersize',5,'color','c'); 
plot(pixelPitchArr, mseOnlyOn(:,1)/255,'-o','linewidth',4,'markersize',5,'color','m'); 
% line([0 70],[50.19 50.19 ]/255,'linestyle','--','linewidth',3,'color','g');
% % plot(70, 50.19/255,'linewidth',3,'linestyle','--','markersize',5,'color','g');
% %  scatter(70, 50.19/255,50,'og','filled')
% line([0 70],[34.6454 34.6454 ]/255,'linestyle','--','linewidth',3,'color','b');

% line([0 70],[ 42.825  42.825]/255,'linestyle','--','linewidth',3,'color','g');

line([0 70],[   38.4046    38.4046]/255,'linestyle','--','linewidth',3,'color','g');

% line([0 70],[37.4126 37.4126]/255,'linestyle','--','linewidth',3,'color','c');
line([0 70],[33.7504 33.7504]/255,'linestyle','--','linewidth',3,'color','c');

% line([0 70],[ 29.5275 29.5275]/255,'linestyle','--','linewidth',3,'color','b');
line([0 70],[ 28.1300 28.1300 ]/255,'linestyle','--','linewidth',3,'color','b');
axis([0 75 0 85/255]);
grid on;
legend('shuffled control','no learning','no learning, only ON', 'prosthesis learning','healthy only ON','healthy','location','se');
%  scatter(70, 50.19/255,50,'og','filled')
set(gca,'fontsize',16);
xlabel('Electrode Spacing (um)','fontsize',16);
 ylabel('Normalized RMS Error','fontsize',16);
 title('Reconstruction Error vs. Electrode Spacing','fontsize',16);
 
 %%
 
 load('mse1_new.mat')
 figure; scatter(squeeze(mse1(1,1,:))/255,squeeze(mse1(4,1,:))/255); 
 hold on; line([0 100]/255,[0 100]/255,'color','k','linewidth',2); 
axis equal; grid on;
axis([0 0.4 0 0.4]);
 
f1 = squeeze(mse1(1,1,:))/255\squeeze(mse1(4,1,:))/255
hold on; plot(.01:.01:.4,f1*[.01:.01:.4],'r','linewidth',3);

set(gca,'fontsize',16);
xlabel('RMSE, Electrode Spacing = 70 microns','fontsize',16);
 ylabel('RMSE, Electrode Spacing = 8.75 microns','fontsize',16);
 title('Reconstruction Error vs. Electrode Spacing','fontsize',16);
 
 %%
 
 
load('mse_35um_prosthesis.mat'); mse35 = mse;
load('mse_70um_prosthesis.mat'); mse70 = mse;
 
 figure; scatter(squeeze(mse70(1,1,:))/255,squeeze(mse35(1,1,:))/255); 
 hold on; line([0 100]/255,[0 100]/255,'color','k','linewidth',2); 
axis equal; grid on;
axis([0 0.4 0 0.4]);
 
f1 = squeeze(mse70(1,1,:))/255\squeeze(mse35(1,1,:))/255
hold on; plot(.01:.01:.4,f1*[.01:.01:.4],'r','linewidth',3);

xlabel('RMSE, Electrode Spacing = 70 microns','fontsize',16);
 ylabel('RMSE, Electrode Spacing = 35 microns','fontsize',16);
 title(sprintf('Recon Error vs. Electrode Spacing, Prosthesis Learning\nRegression Slope %1.2f',f1),'fontsize',16);
 
set(gca,'fontsize',16);

% title(sprintf('Recon Error vs. Electrode Spacing, No Learning, Only On\nRegression Slope %1.2f',.93),'fontsize',16);
 