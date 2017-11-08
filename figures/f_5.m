% f_5
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 5 comparison of reconstuction accuracy for prostheses with different
% parameters.
%
% Compares the reconstruction accuracy for the prostheses with different
% parameters: 35 vs. 70 um pixel pitch, on only activation with no
% learning, etc.
%
% 2017 JRG (c) isetbio team
% [formerlhy f4_dscatter_sep20]

%% Load data - no learning

cd('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/sep20_err')
load('errmean35_nolearn_sep20.mat')
errmean35 = errmeanout;
load('errmean70_nolearn_sep20.mat')
errmean70 = errmeanout;

zi = find(errmeanout~=0);

%% Plot data points - no learning
% dscatter shows density of scatter plot

figure; 
set(gcf,'position',[20         361        1408         413]);
subplot(131);
hold on;
dscatter(sqrt(errmean70(zi)')/255,sqrt(errmean35(zi)')/255);
dscatter(sqrt(errmean70(zi)')/255,sqrt(errmean35(zi)')/255,'plottype','contour',...
    'filled',1,'bins',4*[64 64],'msize',40,'smoothing',20);

hold on; line([0 0.35],[0 0.35],'color','k','linewidth',2);
grid on; axis square
axis([0 .3 0 .3]);

f1=sqrt(errmean70')/255\sqrt(errmean35')/255;
hold on; plot(.01:.01:.4,f1*[.01:.01:.4],'r','linewidth',1.5);

title(sprintf('Reconstruction Error, No Learning \nRegression Slope %1.2f',f1),'fontsize',16);

xlabel(sprintf('70 \\mum pixel spacing'));
ylabel(sprintf('35 \\mum pixel spacing'));
% xlabel(''); ylabel('');
set(gca,'fontsize',14);

set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/prosthesis_35vs70_scatter_nolearning_sep20.pdf']);

%% Load data - no learning, only on
load('errmean35_nolearn_onlyon_sep20.mat')
errmean35 = errmeanout;
load('errmean70_nolearn_onlyon_sep20.mat')
errmean70 = errmeanout;

zi = find(errmeanout~=0);

% figure;
subplot(132);
hold on;

dscatter(sqrt(errmean70(zi)')/255,sqrt(errmean35(zi)')/255);
dscatter(sqrt(errmean70(zi)')/255,sqrt(errmean35(zi)')/255,'plottype','contour','filled',1,'bins',4*[64 64]);

hold on; line([0 0.35],[0 0.35],'color','k','linewidth',2);
grid on; axis square
axis([0 .3 0 .3]);

f1=sqrt(errmean70')/255\sqrt(errmean35')/255;

hold on; plot(.01:.01:.4,f1*[.01:.01:.4],'r','linewidth',1.5);

title(sprintf('Reconstruction Error, No Learning, Only On \nRegression Slope %1.2f',f1),'fontsize',16);

xlabel(sprintf('70 \\mum pixel spacing'));
ylabel(sprintf('35 \\mum pixel spacing'));
set(gca,'fontsize',14);

set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/prosthesis_35vs70_scatter_nolearning_on_sep20.pdf']);
%% Load data - learning

% load('errmean35_2.mat')
% errmean35 = errmean;
% load('errmean70_2.mat')
% errmean70 = errmean;

load('errmean35_learn_sep20.mat')
errmean35 = errmeanout;
load('errmean70_learn_sep20.mat')

errmean70 = errmeanout;

zi = find(errmeanout~=0);

subplot(133);
hold on;

dscatter(sqrt(errmean70(zi)')/255,sqrt(errmean35(zi)')/255);
dscatter(sqrt(errmean70(zi)')/255,sqrt(errmean35(zi)')/255,'plottype','contour','bins',4*[64 64]);

grid on;
axis([0 .3 0 .3]); axis square
hold on; line([0 0.35],[0 0.35],'color','k','linewidth',2);

f1=sqrt(errmean70')/255\sqrt(errmean35')/255;

hold on; plot(.01:.01:.4,f1*[.01:.01:.4],'r','linewidth',1);

title(sprintf('Reconstruction Error, Learning, Only On \nRegression Slope %1.2f',f1),'fontsize',16);
xlabel(sprintf('70 \\mum pixel spacing'));
ylabel(sprintf('35 \\mum pixel spacing'));
set(gca,'fontsize',14);

set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/august24/figures827/prosthesis_35vs70_scatter_learning_sep20.pdf']);