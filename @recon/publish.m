function obj = publish(obj)
% Publish trained filters for recon object to RDT
% 
% 

%%
rd = RdtClient('isetbio');
rd.credentialsDialog;

% rd.crp('/resources/data/istim');
rd.crp('/resources/data/reconstruction/training');

%%
for ii = 101:360
% fname = [reconstructionRootPath '/dat/im_blocks/im_block_1.mat'];
ii
fname = [reconstructionRootPath '/dat/im_blocks/im_block_' num2str(ii) '.mat'];
% fname = [reconstructionRootPath '/dat/nov_results/healthyTest/stim500_1.mat'];
% fname = ('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/gratings_sep22/prosthesis_gratings_sep20_accuracy.mat');
% fname = ('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/sep20_err/errmeanall_sep20.mat');
% fname = ('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/landolt_sep20/prosthesis_landolt_sep20_accuracy.mat');
% fname = ('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug27prima70/filtersmosaic0_sv5_w1_sh4_dr0_pitch_70_decay_2_aug27.mat')
% fname = ('/Users/james//MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv05_w1_sh4_dr0_pitch_70_decay_2_aug29.mat');
% fname = fullfile(reconstructionRootPath, 'dat/', [obj.filterFile '.mat']);
% fname = fullfile(reconstructionRootPath, 'dat/current/', 'filters_mosaic0_sv75_w1_sh2_may26primaSmall.mat');
% save(fname,'iStim');

% Publish to RDT
rd.publishArtifact(fname); 
end