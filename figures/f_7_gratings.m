% f_7_gratings
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 7 psychometric curve for gratings contrast detection.
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

rd = RdtClient('isetbio');
rd.crp('/resources/data/reconstruction');
filterFile = 'healthy_gratings_sep20_accuracy.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
contrastArr = data.contrastArr;
healhyAccuracy = data.Pbig;

% freqArr2(1,1,:) = contrastArr;

radVal=1;
[~,sortInd] = sort(contrastArr(1:end),'ascend');
contrastArrSorted = contrastArr((sortInd));

healhyAccuracySorted = squeeze(mean(healhyAccuracy(radVal,:,:)));

%% Fit psychometric curve - healthy reconstruction

xData = contrastArrSorted(:);
yData = healhyAccuracySorted(sortInd);
% Set up fittype and options.
ft = fittype( '1 - 0.5*exp(-(x/a)^b + c)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.1 0.976303691832645 0];
opts.lower = [1e-6 -1 0]; opts.upper = [1 50 0];

[fitresult, gof] = fit( xData(:), yData(:)/100, ft, opts );

%% Plot - healthy reconstruction

figure;
% hold on;
hold on;
set(gca,'xscale','log')
hhealthy = plot( fitresult, 'b', xData, yData/100 );
set(hhealthy(2),'linewidth',3);
set(hhealthy(1),'markersize',30);
hold on; scatter(xData(:),yData(:)/100,40,'filled');

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

rd = RdtClient('isetbio');
rd.crp('/resources/data/reconstruction');
filterFile = 'prosthesis_gratings_sep20_accuracy.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
contrastArr = data.contrastArr;
prosthesisAccuracy = data.Pbig;

% freqArr2(1,1,:) = contrastArr;
% 
% contAdd = ones(1,1,6);
% contAdd(1,1,:) = [55.9375   63.4375   62.1875   69.6875   75.0000   89.6875];
% Pbig(1,:,17:22)=contAdd;
% contrastArr(17:22) = [.002:.001:.007];

radVal=1;
[~,sortInd] = sort(contrastArr(1:end),'ascend');
contrastArrSorted = contrastArr((sortInd));

prosthesisAccuracySorted = squeeze((prosthesisAccuracy(radVal,:,:)));

%% Fit psychometric curve - prosthesis reconstruction
% [xData, yData] = prepareCurveData(freqArrSort(:),Pbigsort(sortInd));
xData = contrastArrSorted(:);
yData = prosthesisAccuracySorted(sortInd);
% Set up fittype and options.
ft = fittype( '1 - 0.5*exp(-(x/a)^b + c)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.1 0.976303691832645 0];
opts.lower = [1e-6 -1 0]; opts.upper = [1 50 0];

[fitresult, gof] = fit( xData(:), yData(:)/100, ft, opts );


%% Plot - prosthesis reconstruction

hold on;

grid on;
set(gca,'xScale','log');
set(gca,'fontsize',14)
legend('Healthy','Prosthesis (learning)','location','nw');
title('Landolt C Orientation Discrimination');
axis([5e-6 .1 .45 1]);
% title(''); xlabel(''); ylabel('');
hold on;
hprosthesis = plot( fitresult,'r', xData, yData/100);
set(hprosthesis(2),'linewidth',3);

set(hprosthesis(1),'markersize',30,'color','r');
hold on; scatter(xData(:),yData(:)/100,40,'filled');
set(gca,'xscale','log')
legend([hhealthy(2), hprosthesis(2)],{'Healthy', 'Prosthesis'} ,'Location', 'NorthWest')

grid on
thresh1 = fitresult.a;
title(sprintf('Landolt C Gap Detection, thresholds %1.4f, %1.4f',thresh1h,thresh1));
set(gca,'fontsize',16)

xlabel('Contrast'); ylabel('Fraction Correct');
%%
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/landolt_psycometric_oct4.pdf']);
