


clear
% load('stats_6.mat')
ccB = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima35test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/prima35TestLearning/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;



ccB2 = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima70test/results/stats_' num2str(ii) '.mat']); ccB2 = [ccB2; ccAll]; end;
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/primaTestLearning/stats_' num2str(ii) '.mat']); ccB2 = [ccB2; ccAll]; end;
%%
figure; dscatter(ccB,ccB2); 

hold on; dscatter(ccB,ccB2,'plottype','contour'); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
axis([0 1 0 1]); grid on; colormap parula
% colormap jet
% set(gca,'xscale','log');set(gca,'yscale','log')
% figure; dscatter(ccB,ccB2); hold on; dscatter(ccB,ccB2,'plottype','contour'); axis([0 1 0 1]); grid on; colormap parula

%%
ccB = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/prima35TestLearning/stats_' num2str(ii) '.mat']); ccB = [ccB mseAll]; end;

ccB2 = [];
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/primaTestLearning/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;

%%
figure; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)'),sqrt(ccB2(ccB~=0|ccB2~=0)')); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
hold on; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)'),sqrt(ccB2(ccB~=0|ccB2~=0)'),'plottype','contour'); 
% axis([0 1 0 1]); 
grid on; colormap parula
% colormap jet

%%
figure; dscatter((ccB(ccB~=0|ccB2~=0)'),(ccB2(ccB~=0|ccB2~=0)')); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
hold on; dscatter((ccB(ccB~=0|ccB2~=0)'),(ccB2(ccB~=0|ccB2~=0)'),'plottype','contour'); 
% axis([0 1 0 1]); 
grid on; colormap parula
% colormap jet

% set(gca,'xscale','log');set(gca,'yscale','log')
%%
figure; dscatter(log(1+(ccB(ccB~=0|ccB2~=0)')),log(1+(ccB2(ccB~=0|ccB2~=0)'))); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
hold on; dscatter(log(1+(ccB(ccB~=0|ccB2~=0)')),log(1+(ccB2(ccB~=0|ccB2~=0)')),'plottype','contour'); 
% axis([0 1 0 1]); 
grid on; colormap parula
% set(gca,'xscale','log');set(gca,'yscale','log')
set(gca,'xscale','log');set(gca,'yscale','log')
% axis([4 9 4 9]);
% line([4 10],[4 10],'color','r','linewidth',4);
%%

ccB = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima35test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/healthyTest/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;



ccB2 = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima70test/results/stats_' num2str(ii) '.mat']); ccB2 = [ccB2; ccAll]; end;
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/primaTestLearning/stats_' num2str(ii) '.mat']); ccB2 = [ccB2; ccAll]; end;

figure; dscatter(ccB,ccB2); 
hold on; dscatter(ccB,ccB2,'plottype','contour'); 
axis([0 1 0 1]); grid on; colormap parula
% colormap jet
% set(gca,'xscale','log');set(gca,'yscale','log')
% figure; dscatter(ccB,ccB2); hold on; dscatter(ccB,ccB2,'plottype','contour'); axis([0 1 0 1]); grid on; colormap parula

%%

ccB = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima35test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/healthyTest/stats_' num2str(ii) '.mat']); ccB = [ccB mseAll]; end;



ccB2 = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima70test/results/stats_' num2str(ii) '.mat']); ccB2 = [ccB2; ccAll]; end;
for ii = 2:20; load(['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/nov_results/primaTestLearning/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;
%%
figure; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)')/255,sqrt(ccB2(ccB~=0|ccB2~=0)')/255); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
hold on; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)')/255,sqrt(ccB2(ccB~=0|ccB2~=0)')/255,'plottype','contour'); 
% axis([0 1 0 1]); 
grid on; colormap parula
% colormap jet
% axis square; 
% axis equal
% set(gca,'xscale','log');set(gca,'yscale','log')
%%
figure; dscatter(log10((ccB(ccB>0|ccB2>0)')/255),log10((ccB2(ccB~=0|ccB2~=0)')/255)); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
hold on; dscatter(log10((ccB(ccB>0|ccB2>0)')/255),log10((ccB2(ccB~=0|ccB2~=0)')/255),'plottype','contour'); 
% axis([0 1 0 1]); 
grid on; colormap parula
% colormap jet
% axis square; 
axis equal
% set(gca,'xscale','log');set(gca,'yscale','log')
%%
figure; dscatter(log(1+(ccB(ccB~=0|ccB2~=0)')),log(1+(ccB2(ccB~=0|ccB2~=0)'))); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
hold on; dscatter(log(1+(ccB(ccB~=0|ccB2~=0)')),log(1+(ccB2(ccB~=0|ccB2~=0)')),'plottype','contour'); 
% axis([0 1 0 1]); 
grid on; colormap parula
% set(gca,'xscale','log');set(gca,'yscale','log')
% set(gca,'xscale','log');set(gca,'yscale','log')
% axis([4 9 4 9]);
line([4 10],[.04 10],'color','r','linewidth',4);

set(gca,'xscale','log');set(gca,'yscale','log')