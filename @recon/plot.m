function plot(obj, varargin)
%PLOT - plot the properties of the recon object

load(obj.filterFile);
figure; for fr = 1:25; subplot(5,5,fr); imagesc(reshape(filterMat(000+fr,:),[100 100]));colormap parula;  end;