


clear
% load('stats_6.mat')
ccB = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima35test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;

ccB2 = [];
for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima70test/results/stats_' num2str(ii) '.mat']); ccB2 = [ccB2; ccAll]; end;
figure; dscatter(ccB,ccB2); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
axis([0 1 0 1]); grid on; colormap parula
colormap jet
% set(gca,'xscale','log');set(gca,'yscale','log')
% figure; dscatter(ccB,ccB2); hold on; dscatter(ccB,ccB2,'plottype','contour'); axis([0 1 0 1]); grid on; colormap parula

ccB = [];
% for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20test/results/stats_' num2str(ii) '.mat']); ccB = [ccB; ccAll]; end;
for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima35test/results/stats_' num2str(ii) '.mat']); ccB = [ccB mseAll]; end;

ccB2 = [];
for ii = 2:20; load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/sep20prima70test/results/stats_' num2str(ii) '.mat']); ccB2 = [ccB2 mseAll]; end;
figure; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)'),sqrt(ccB2(ccB~=0|ccB2~=0)')); 
% hold on; dscatter(ccB,ccB2,'plottype','contour'); 
hold on; dscatter(sqrt(ccB(ccB~=0|ccB2~=0)'),sqrt(ccB2(ccB~=0|ccB2~=0)'),'plottype','contour'); 
% axis([0 1 0 1]); 
grid on; colormap parula
colormap jet
