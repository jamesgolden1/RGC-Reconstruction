function obj = publish(obj)
% Publish trained filters for recon object to RDT
% 
% 

%%
rd = RdtClient('isetbio');
rd.credentialsDialog;
rd.crp('/resources/data/istim');

fname = fullfile(reconstructionRootPath, 'dat/', [obj.filterFile '.mat']);
% fname = fullfile(reconstructionRootPath, 'dat/current/', 'filters_mosaic0_sv75_w1_sh2_may26primaSmall.mat');
% save(fname,'iStim');

% Publish to RDT
rd.publishArtifact(fname); 