% f_9_gratings
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 9: psychometric curve for gratings frequency detection.
%
% Plots the psychometric curves for gratings contrast detection. Fits
% datapoints to a psychometric function.
%
% 2017 JRG (c) isetbio team
% [formerly make_psychometric_curve_gratings_sep22_fit]

%% Load data and set up plot - healthy reconstruction

clear

% load('ws_gratings_healthy_sep20_4reps_filt05_comb.mat','contrastArr','Pbig');
% load('healthy_gratings_sep20_accuracy.mat');

% rd = RdtClient('isetbio');
% rd.crp('/resources/data/reconstruction');
% filterFile = 'healthy_gratings_sep20_accuracy.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% contrastArr = data.contrastArr;
% healthyAccuracy = data.Pbig;

load('ws_gratings_healthy_aug17_2018_4reps_out0.mat');
contrastArr = freqArr;
healthyAccuracy = Pbig;

% freqArr2(1,1,:) = contrastArr;

radVal=1;
[~,sortInd] = sort(contrastArr(1:end),'ascend');
contrastArrSorted = contrastArr((sortInd))/1.7;

healthyAccuracySorted = squeeze(mean(healthyAccuracy(radVal,:,:)));

%% Fit psychometric curve - healthy reconstruction

xData = contrastArrSorted(:);
yData = healthyAccuracySorted(sortInd);
% Set up fittype and options.
ft = fittype( '0.5*exp(-(x/a)^b + c) + .5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [2 0.976303691832645 0 ];
opts.lower = [-10 -4 0 ]; opts.upper = [100 50 0 ];

[fitresult, gof] = fit( xData(:), yData(:)/100, ft, opts );

%% Plot - healthy reconstruction

figure;
% hold on;
hold on;
set(gca,'xscale','log')
hhealthy = plot( fitresult, 'b', xData, yData/100 );
set(hhealthy(2),'linewidth',3,'linestyle','--');
set(hhealthy(1),'markersize',30);
hold on; scatter(xData(:),yData(:)/100,40,'b','filled');

hones = plot(xData(:),ones(size(xData(:))),'b','linewidth',3);
scatter(xData(:),ones(size(xData(:))),80,'b','filled');

grid on
thresh1h = fitresult.a;
set(gca,'fontsize',16)

xlabel('Contrast'); ylabel('Fraction Correct');
% title('
grid on;
set(gca,'xScale','log');
set(gca,'fontsize',14)
% axis([4e-4 .2 48 100]);
legend('Healthy','Prosthesis','location','nw');
title('Landolt C Orientation Discrimination');

%% Load data and set up plot - prosthesis reconstruction
% load('ws_gratings_prima_sep17_4reps_out0.mat','contrastArr','Pbig');
% load('prosthesis_gratings_sep20_accuracy.mat');

% rd = RdtClient('isetbio');
% rd.crp('/resources/data/reconstruction');
% filterFile = 'prosthesis_gratings_sep20_accuracy.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% contrastArr = data.contrastArr;
% prosthesisAccuracy = data.Pbig;


% load('ws_gratings_prima_aug17_2018_4reps_out0.mat');
load('ws_gratings_primaLearning_aug25_2018_8reps500_sh8_out0.mat');


freqArr = freqArr(1:end);
Pbig = Pbig(:,:,1:end);
contrastArr = freqArr;
prosthesisAccuracy = Pbig;

% freqArr2(1,1,:) = contrastArr;
% 
% contAdd = ones(1,1,6);
% contAdd(1,1,:) = [55.9375   63.4375   62.1875   69.6875   75.0000   89.6875];
% Pbig(1,:,17:22)=contAdd;
% contrastArr(17:22) = [.002:.001:.007];

radVal=1;
[~,sortInd] = sort(contrastArr(1:end),'ascend');
sortInd = sortInd([1:7 9 11]);
% sortInd = sortInd([1:7  11:12 14:16]);
% sortInd = sortInd([1:7 9 11:12 14:15]);
contrastArrSorted = contrastArr((sortInd))/1.7;

prosthesisAccuracySorted = squeeze(mean(prosthesisAccuracy(radVal,:,sortInd)));
% figure; plot(contrastArrSorted,prosthesisAccuracySorted)
%% Fit psychometric curve - prosthesis reconstruction
% [xData, yData] = prepareCurveData(freqArrSort(:),Pbigsort(sortInd));
xData = contrastArrSorted(:);
yData = prosthesisAccuracySorted(:);
% Set up fittype and options.
ft = fittype( '0.5*exp(-(x/a)^b + c) + .5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [2 0.976303691832645 0 ];
opts.lower = [-10 -4 0 ]; opts.upper = [100 50 0 ];

[fitresult, gof] = fit( xData(:), yData(:)/100, ft, opts );


%% Plot - prosthesis reconstruction
% figure;
hold on;

grid on;
set(gca,'xScale','log');
set(gca,'fontsize',14)
legend('Healthy','Prosthesis (learning)','location','ne');
% title('Landolt C Orientation Discrimination');
% axis([5e-6 .1 .45 1]);
% title(''); xlabel(''); ylabel('');
hold on;
hprosthesis = plot( fitresult,'r', xData, yData/100);
set(hprosthesis(2),'linewidth',3);

set(hprosthesis(1),'markersize',30,'color','r');
hold on; scatter(xData(:),yData(:)/100,40,'r','filled');

hprosthesisnl2 = plot( fitresult,'r', [xData; (3:3:30)'], fitresult([xData; (3:3:30)']));

set(hprosthesisnl2(2),'linewidth',3,'markersize',.5);

set(gca,'xscale','log')
% legend([hhealthy(2), hprosthesis(2)],{'Healthy', 'Prosthesis'} ,'Location', 'NorthWest')

grid on
thresh1 = fitresult.a;
% title(sprintf('Gratings Orientation Discrimination, thresholds %1.4f, %1.4f',thresh1h,thresh1));
set(gca,'fontsize',16)

xlabel('Contrast'); ylabel('Fraction Correct');


%% Load data and set up plot - prosthesis Only On reconstruction
% load('ws_gratings_prima_sep17_4reps_out0.mat','contrastArr','Pbig');
% load('prosthesis_gratings_sep20_accuracy.mat');

% rd = RdtClient('isetbio');
% rd.crp('/resources/data/reconstruction');
% filterFile = 'prosthesis_gratings_sep20_accuracy.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% contrastArr = data.contrastArr;
% prosthesisAccuracy = data.Pbig;


% load('ws_gratings_primaOnlyOn_aug17_2018_4reps_out0.mat');
% 
load('ws_gratings_primaOnlyOnCorrect_aug25_2018_8reps500_sh8_f05_out0.mat');
freqArr = freqArr(1:end);
Pbig = Pbig(:,:,1:end);
contrastArr = freqArr;
prosthesisAccuracy = Pbig;

% freqArr2(1,1,:) = contrastArr;
% 
% contAdd = ones(1,1,6);
% contAdd(1,1,:) = [55.9375   63.4375   62.1875   69.6875   75.0000   89.6875];
% Pbig(1,:,17:22)=contAdd;
% contrastArr(17:22) = [.002:.001:.007];

radVal=1;
% [~,sortInd] = sort(contrastArr(1:end),'ascend');
sortInd = sortInd([1:7  11]);
% sortInd = sortInd([1:7  11:12 14:16]);
% sortInd = sortInd([1:7 9 11:12 14:16]);
contrastArrSorted = contrastArr((sortInd))/1.7;

prosthesisAccuracySorted = squeeze(mean(prosthesisAccuracy(radVal,:,sortInd)));

contrastArrSorted = [.1 contrastArrSorted];
prosthesisAccuracySorted = [100; prosthesisAccuracySorted];
% figure; plot(contrastArrSorted,prosthesisAccuracySorted)
%% Fit psychometric curve - prosthesis Only On reconstruction
% [xData, yData] = prepareCurveData(freqArrSort(:),Pbigsort(sortInd));
xData = contrastArrSorted(:);
yData = prosthesisAccuracySorted(:);
% Set up fittype and options.
ft = fittype( '0.5*exp(-(x/a)^b + c) + .5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [2 0.976303691832645 0 ];
opts.lower = [-10 -4 0 ]; opts.upper = [100 50 0 ];

[fitresult, gof] = fit( xData(:), yData(:)/100, ft, opts );


%% Plot - prosthesis Only On reconstruction
% figure;
hold on;

grid on;
set(gca,'xScale','log');
set(gca,'fontsize',14)
legend('Healthy','Prosthesis (learning)','location','ne');
% title('Landolt C Orientation Discrimination');
% axis([5e-6 .1 .45 1]);
% title(''); xlabel(''); ylabel('');
hold on;
hprosthesisoo = plot( fitresult,'g', xData, yData/100);
set(hprosthesisoo(2),'linewidth',3);

set(hprosthesisoo(1),'markersize',30,'color','g');
hold on; scatter(xData(:),yData(:)/100,40,'g','filled');


hprosthesisoo2 = plot( fitresult,'g', [xData; (3:3:30)'], fitresult([xData; (3:3:30)']));

set(hprosthesisoo2(2),'linewidth',3,'markersize',.5);


set(gca,'xscale','log')
% legend([hhealthy(2), hprosthesis(2)],{'Healthy', 'Prosthesis'} ,'Location', 'NorthWest')

grid on
thresh1oo = fitresult.a;
title(sprintf('Gratings Orientation Discrimination, thresholds %1.4f, %1.4f',thresh1h,thresh1));
set(gca,'fontsize',16)

xlabel('Contrast'); ylabel('Fraction Correct');


%% Load data and set up plot - prosthesis No Learning reconstruction
% load('ws_gratings_prima_sep17_4reps_out0.mat','contrastArr','Pbig');
% load('prosthesis_gratings_sep20_accuracy.mat');

% rd = RdtClient('isetbio');
% rd.crp('/resources/data/reconstruction');
% filterFile = 'prosthesis_gratings_sep20_accuracy.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% contrastArr = data.contrastArr;
% prosthesisAccuracy = data.Pbig;



load('ws_gratings_primaNoLearn_aug17_2018_4reps_out0.mat');

freqArr = freqArr(2:end);
Pbig = Pbig(:,:,2:end);
contrastArr = freqArr;
prosthesisAccuracy = Pbig;

% freqArr2(1,1,:) = contrastArr;
% 
% contAdd = ones(1,1,6);
% contAdd(1,1,:) = [55.9375   63.4375   62.1875   69.6875   75.0000   89.6875];
% Pbig(1,:,17:22)=contAdd;
% contrastArr(17:22) = [.002:.001:.007];

radVal=1;
[~,sortInd] = sort(contrastArr(1:end),'ascend');
contrastArrSorted = contrastArr((sortInd))/1.7;

prosthesisAccuracySorted = squeeze(mean(prosthesisAccuracy(radVal,:,:)));
% figure; plot(contrastArrSorted,prosthesisAccuracySorted)
%% Fit psychometric curve - prosthesis No Learning reconstruction
% [xData, yData] = prepareCurveData(freqArrSort(:),Pbigsort(sortInd));
xData = contrastArrSorted(:);
yData = prosthesisAccuracySorted(sortInd);
% Set up fittype and options.
ft = fittype( '0.5*exp(-(x/a)^b + c) + .5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [2 0.976303691832645 0 ];
opts.lower = [-10 -4 0 ]; opts.upper = [100 50 0 ];

[fitresult, gof] = fit( xData(:), yData(:)/100, ft, opts );


%% Plot - prosthesis No Learning reconstruction
% figure;
hold on;

grid on;
set(gca,'xScale','log');
set(gca,'fontsize',14)
legend('Healthy','Prosthesis (learning)','location','ne');
% title('Landolt C Orientation Discrimination');
axis([.29 37 .45 1]);
% title(''); xlabel(''); ylabel('');
hold on;
hprosthesisnl = plot( fitresult,'k', xData, yData/100);
set(hprosthesisnl(2),'linewidth',3);

set(hprosthesisnl(1),'markersize',30,'color','k');
hold on; scatter(xData(:),yData(:)/100,40,'k','filled');
set(gca,'xscale','log')
legend([hhealthy(2), hones, hprosthesis(2),hprosthesisoo(2),hprosthesisnl(2)],{'healthy, low contrast','healthy', 'prosthesis, learned', 'prosthesis, only ON', 'prosthesis, ON+OFF'} ,'Location', 'ne','fontsize',12)

grid on
thresh1nl = fitresult.a;
% title(sprintf('Gratings Orientation Discrimination, \nthresholds %1.4f, %1.4f, %1.4f, %1.4f',thresh1h,thresh1,thresh1oo,thresh1nl));
title('');
set(gca,'fontsize',16)

% xlabel('Frequency (cpd)'); ylabel('Fraction Correct');
        xlabel('');
    ylabel('');
    
    %%
        hc = get(gca,'Children');
        set(gca,'Children',[hc(11:13); hc(6:10);  hc(1:5); hc(14:18)]);
%     set(gca,'Children',[hc(10:14) hc(7:9)  hc(4:6)   hc(1:3)])
%     set(gca,'Children',[hc(7:9) ; hc(4:6) ; hc(1:3)  ; hc(10:14)])
%%
% set(gcf,'position',[ 440   472   663   326]);
set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/gratings_freq_JNE2_826.pdf']);
