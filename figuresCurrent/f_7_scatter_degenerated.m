
% f_7_scatter_degenerated
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 7 comparison of reconstuction accuracy for prostheses with different
% parameters for a degenerated retina.
%
% Compares the reconstruction accuracy for the prostheses with different
% parameters: 35 vs. 70 um pixel pitch, on only activation with no
% learning, etc.
%
% 2018 JRG (c) isetbio team
% [formerlhy resultsAll_...degen_aug26]

%%
clear;

%% Build example reconstructions

% pRecon.degenFlag = 1;
% reconHealthy = recon(pRecon);
% reconHealthy.buildPrimaLeaf();%'pixelWidth',35);
% % reconHealthy.buildPrimaWheel('pixelWidth',70);

%% 

datDir = [reconstructionRootPath '/figuresCurrent/dat/f_7_scatter_degenerated/'];

plotLin = 1;
degenFlag = '_degen';
dirName35 = 'prosthesis_35_testing_aug22_degen';
dirName70 = 'prosthesis_70_testing_aug22_degen';

% dirName35 = 'prosthesis_35_testing_aug26_degen';
% dirName70 = 'prosthesis_70_testing_aug26_degen';


%% Load data - RDT

% rd = RdtClient('isetbio');
% rd.crp('/resources/data/reconstruction');
% filterFile = 'errmeanall_sep20.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');

%% Load data - learning
ccB = [];

for ii = [1:3 5]; load([datDir dirName70 '/results/stats_' num2str(ii) '.mat']); ccB = [ccB mseAll]; end;
for ii = [4]; load([datDir dirName70 '/results/stats_' num2str(ii) '.mat']); ccB = [ccB mseAll]; end;

ccB2 = [];

for ii = [1:3 5]; load([datDir dirName35 '/results/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;
for ii = [4]; load([datDir dirName35 '/results/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;

%% Plot

if plotLin == 1
    %%
    % ccB(isnan(ccB))=0; %ccB2 = ccB2(1:end-1);
    ccB(isnan(ccB))=0; ccB2 = ccB2(1:end-1);
    figure; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)')/255,sqrt(ccB2(ccB~=0|ccB2~=0)')/255);
    % hold on; dscatter(ccB,ccB2,'plottype','contour');
    % hold on; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)'),sqrt(ccB2(ccB~=0|ccB2~=0)'),'plottype','contour');
    % axis([0 1 0 1]);
    
    hold on;
    scatter(sqrt(mean(ccB(ccB~=0|ccB2~=0)'))/255,sqrt(mean(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    % line([1 1],[100 100])
    grid on; colormap parula
    % colormap jet
    line([0 200],[0 200],'color','k')
    
    f1=sqrt((ccB(ccB~=0|ccB2~=0)'))/255\sqrt(ccB2(ccB~=0|ccB2~=0)')/255;
    hold on; plot(.01:.01:.6,f1*[.01:.01:.6],'r','linewidth',1.5);
    [sqrt(mean(ccB(ccB~=0)))/255, sqrt(mean(ccB2(ccB2~=0)))/255, f1]
    % xlabel(sprintf('reconstruction RMS error (70 \\mum)'),'fontsize',14)
    % ylabel(sprintf('reconstruction RMS error (35 \\mum)'),'fontsize',14)
    
    set(gca,'xticklabel','');
    
    set(gca,'yticklabel','');
    axis square; axis([0 0.6 0 0.6])
    % line([1e-5 1e-5],[0 0])
    % set(gca,'xscale','log');set(gca,'yscale','log')
    
    set(gcf,'PaperPositionMode','auto')
    % print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/figures_aug2018/prosthesis_35vs70_scatter_learning_aug24_2018' degenFlag '.pdf']);
    %
    
    %%
else
    figure; dscatter((ccB(ccB~=0|ccB2~=0)')/255^2,(ccB2(ccB~=0|ccB2~=0)')/255^2,'logy',1);
    % hold on; dscatter(ccB,ccB2,'plottype','contour');
    % hold on; dscatter((ccB(ccB~=0|ccB2~=0)')/255^2,(ccB2(ccB~=0|ccB2~=0)')/255^2,'plottype','contour','logy',1);
    % axis([0 1 0 1]);
    hold on;
    % scatter(mean((ccB(ccB~=0|ccB2~=0)')/255^2),mean((ccB2(ccB~=0|ccB2~=0)')/255^2), 32, 'r','filled');
    
    % scatter(sqrt(mean(ccB(ccB~=0|ccB2~=0)'))/255,sqrt(mean(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    line(10.^[-5 0],10.^[-5 0],'color','k','linewidth',2);
    
    grid on; colormap parula
    axis([ 0.0001    0.25    0.0001    0.25])
    
    xlabel(sprintf('reconstruction RMS error (70 \\mum)'),'fontsize',14)
    ylabel(sprintf('reconstruction RMS error (35 \\mum)'),'fontsize',14)
end
%%

%% Load data - Only On

ccB = [];
for ii = 1:4; load([datDir dirName70 '/resultsOnlyOnC/stats_' num2str(ii) '.mat']); ccB = [ccB mseAll]; end;

ccB2 = [];
for ii = 1:4; load([datDir dirName35 '/resultsOnlyOnC/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;

%% PLot

if plotLin == 1
    %%
    figure; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)')/255,sqrt(ccB2(ccB~=0|ccB2~=0)')/255);
    % hold on; dscatter(ccB,ccB2,'plottype','contour');
    % hold on; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)'),sqrt(ccB2(ccB~=0|ccB2~=0)'),'plottype','contour');
    % axis([0 1 0 1]);
    
    hold on;
    % scatter(mean(sqrt(ccB(ccB~=0|ccB2~=0)'))/255,mean(sqrt(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    
    scatter(sqrt(mean(ccB(ccB~=0|ccB2~=0)'))/255,sqrt(mean(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    % line([1 1],[100 100])
    grid on; colormap parula
    % colormap jet
    
    line([0 200],[0 200],'color','k')
    
    f1=sqrt((ccB(ccB~=0|ccB2~=0)'))/255\sqrt(ccB2(ccB~=0|ccB2~=0)')/255;
    [sqrt(mean(ccB(ccB~=0)))/255, sqrt(mean(ccB2(ccB2~=0)))/255, f1]
    hold on; plot(.01:.01:.6,f1*[.01:.01:.6],'r','linewidth',1.5);
    
    % xlabel(sprintf('reconstruction RMS error (70 \\mum)'),'fontsize',14)
    % ylabel(sprintf('reconstruction RMS error (35 \\mum)'),'fontsize',14)
    axis square; axis([0 0.6 0 0.6])
    
    % set(gca,'fontsize',14);
    % xlabel('RMS error 70 um'); ylabel('RMS error 35 um');
    set(gca,'xticklabel','');
    
    set(gca,'yticklabel','');
    % line([1 100],[1 100],'color','k')
    % line([1e-5 1e-5],[0 0])
    % set(gca,'xscale','log');set(gca,'yscale','log')
    % print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/figures_aug2018/prosthesis_35vs70_scatter_learning_aug13_2018.pdf']);
    
    set(gcf,'PaperPositionMode','auto')
    % print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/figures_aug2018/prosthesis_35vs70_scatter_onlyOn_aug26_2018' degenFlag '.pdf']);
    
    %%
else
    figure; dscatter((ccB(ccB~=0|ccB2~=0)')/255^2,(ccB2(ccB~=0|ccB2~=0)')/255^2,'logy',1);
    % hold on; dscatter(ccB,ccB2,'plottype','contour');
    % hold on; dscatter((ccB(ccB~=0|ccB2~=0)')/255^2,(ccB2(ccB~=0|ccB2~=0)')/255^2,'plottype','contour','logy',1);
    % axis([0 1 0 1]);
    hold on;
    % scatter(mean((ccB(ccB~=0|ccB2~=0)')/255^2),mean((ccB2(ccB~=0|ccB2~=0)')/255^2), 32, 'r','filled');
    
    % scatter(sqrt(mean(ccB(ccB~=0|ccB2~=0)'))/255,sqrt(mean(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    line(10.^[-5 0],10.^[-5 0],'color','k','linewidth',2);
    
    grid on; colormap parula
    axis([ 0.0001    0.25    0.0001    0.25])
    xlabel(sprintf('reconstruction RMS error (70 \\mum)'),'fontsize',14)
    ylabel(sprintf('reconstruction RMS error (35 \\mum)'),'fontsize',14)
end
%%


%% Load data - No Learning

ccB = [];
for ii = [1:3 5]; load([datDir dirName70 '/resultsNoLearn/stats_' num2str(ii) '.mat']); ccB = [ccB mseAll]; end;

ccB2 = [];
% for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/prima35TestLearning/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;
for ii = [1:3 5]; load([datDir dirName35 '/resultsNoLearn/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;

%% Plot

if plotLin == 1
    %%
    figure; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)')/255,sqrt(ccB2(ccB~=0|ccB2~=0)')/255);
    % hold on; dscatter(ccB,ccB2,'plottype','contour');
    % hold on; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)'),sqrt(ccB2(ccB~=0|ccB2~=0)'),'plottype','contour');
    % axis([0 1 0 1]);
    
    hold on;
    % scatter(mean(sqrt(ccB(ccB~=0|ccB2~=0)'))/255,mean(sqrt(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    
    scatter(sqrt(mean(ccB(ccB~=0|ccB2~=0)'))/255,sqrt(mean(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    % line([1 1],[100 100])
    grid on; colormap parula
    % colormap jet
    
    line([0 200],[0 200],'color','k')
    
    
    f1=sqrt((ccB(ccB~=0|ccB2~=0)'))/255\sqrt(ccB2(ccB~=0|ccB2~=0)')/255;
    hold on; plot(.01:.01:.6,f1*[.01:.01:.6],'r','linewidth',1.5);
    [sqrt(mean(ccB(ccB~=0)))/255, sqrt(mean(ccB2(ccB2~=0)))/255, f1]
    % xlabel(sprintf('reconstruction RMS error (70 \\mum)'),'fontsize',14)
    % ylabel(sprintf('reconstruction RMS error (35 \\mum)'),'fontsize',14)
    axis square; axis([0 0.6 0 0.6])
    
    set(gca,'xticklabel','');
    
    set(gca,'yticklabel','');
    % line([1 100],[1 100],'color','r')
    % line([1e-5 1e-5],[0 0])
    % set(gca,'xscale','log');set(gca,'yscale','log')
    set(gcf,'PaperPositionMode','auto')
    % print(gcf, '-dpdf',['/Users/james/Documents/matlab/EJLPhosphene/local/figures_aug2018/prosthesis_35vs70_scatter_noLearning_aug24_2018' degenFlag '.pdf']);
    
    %%
else
    figure; dscatter((ccB(ccB~=0|ccB2~=0)')/255^2,(ccB2(ccB~=0|ccB2~=0)')/255^2,'logy',1);
    % hold on; dscatter(ccB,ccB2,'plottype','contour');
    % hold on; dscatter((ccB(ccB~=0|ccB2~=0)')/255^2,(ccB2(ccB~=0|ccB2~=0)')/255^2,'plottype','contour','logy',1);
    % axis([0 1 0 1]);
    hold on;
    % scatter(mean((ccB(ccB~=0|ccB2~=0)')/255^2),mean((ccB2(ccB~=0|ccB2~=0)')/255^2), 32, 'r','filled');
    
    % scatter(sqrt(mean(ccB(ccB~=0|ccB2~=0)'))/255,sqrt(mean(ccB2(ccB~=0|ccB2~=0)'))/255, 32, 'r','filled');
    line(10.^[-5 0],10.^[-5 0],'color','k','linewidth',2);
    
    grid on; colormap parula
    axis([ 0.0001    0.25    0.0001    0.25])
    xlabel(sprintf('reconstruction RMS error (70 \\mum)'),'fontsize',14)
    ylabel(sprintf('reconstruction RMS error (35 \\mum)'),'fontsize',14)
    
end