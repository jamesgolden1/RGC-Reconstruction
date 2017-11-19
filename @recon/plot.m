function plot(obj, plottype, varargin)
%PLOT - plot the properties of the recon object

%%
p = inputParser;
p.addRequired('obj');
p.addRequired('plottype');
p.addParameter('pixelWidth',70,@isnumeric);
p.KeepUnmatched = true;
p.parse(obj, plottype, varargin{:});
pixelWidth = p.Results.pixelWidth;
%%

switch plottype
    case 'filters'
        
        disp('Loading filter file from RDT...');
        
        % load(obj.filterFile);
        
        if isempty(pixelWidth)
            
            rd = RdtClient('isetbio');
            rd.crp('/resources/data/reconstruction');
            filterFile = 'filtersmosaic0_sv50_w1_sh15_dr0_aug27.mat';
            data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
            filterMat = data.filterMat; clear data;
        else
            %          load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv 5_w1_sh4_dr0_pitch_70_decay_2.mat')
            %
            rd = RdtClient('isetbio');
            rd.crp('/resources/data/reconstruction');
            filterFile = 'filtersmosaic0_sv05_w1_sh4_dr0_pitch_70_decay_2_aug29.mat';
            %         filterFile = 'filtersmosaic0_sv5_w1_sh4_dr0_pitch_70_decay_2_aug27.mat'
            data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
            filterMat = data.filterMat; clear data;
        end
        
        % size(filterMat)
        % figure; for fr = 1:25; subplot(5,5,fr); imagesc(reshape(filterMat(000+fr,:),[100 100]));colormap parula;  end;
        figure;
        for fr = 1:25;
            subplot(5,5,fr);
            filtIm = reshape(filterMat(000+fr,:),[100 100]);
            imagesc(filtIm);
            mabs = max(abs(filtIm(:)));
            caxis([-mabs mabs]); colormap parula;
        end;
        
        % figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterMat(000+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
        
end