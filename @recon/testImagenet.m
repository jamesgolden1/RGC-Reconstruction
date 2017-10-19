function mse = testImagenet(obj,varargin)
% Computes MSe for movie.


p = inputParser;
p.addParameter('stimFile',[],@ischar);
p.addParameter('respFile',[],@ischar);
p.addParameter('filterFile',[],@ischar);
p.addParameter('windowSize',8,@isnumeric);
p.addParameter('percentSV',0.5,@isnumeric);
p.addParameter('mosaicFile',[],@ischar);
p.addParameter('trainFraction',1,@isnumeric);
p.addParameter('shiftTime',0,@isnumeric);
p.addParameter('dropout',0,@isnumeric);
p.addParameter('stimType','ns',@ischar);
p.addParameter('buildFile',[],@ischar);

p.addParameter('pixelWidth',[],@isnumeric);
p.addParameter('currentDecay',2,@isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
filterFile = p.Results.filterFile;
stimFileName = p.Results.stimFile;
respFileName = p.Results.respFile;
windowSize = p.Results.windowSize;
percentSV = p.Results.percentSV;
mosaicFile = p.Results.mosaicFile;
trainSizeArray = p.Results.trainFraction;
shiftTime = p.Results.shiftTime;
stimType = p.Results.stimType;
buildFile = p.Results.buildFile;
dropout = p.Results.dropout;

pixelWidth = p.Results.pixelWidth;
currentDecay = p.Results.currentDecay;

%
% p = inputParser;
% p.addRequired('obj');
% p.addParameter('permuteFlag',0,@isnumeric);
% p.addParameter('testFile',[],@ischar);
% p.addParameter('mosaicFile','mosaic0',@ischar);
% p.KeepUnmatched = true;
% p.parse(obj,varargin{:});
% permuteFlag = p.Results.permuteFlag;
% testFile = p.Results.testFile;
% mosaicFile = p.Results.mosaicFile;

%%

spikeFile = fullfile(reconstructionRootPath,'dat',stimFileName,['sp_hallway_' mosaicFile '.mat']);
load(spikeFile);

stimFile = fullfile(reconstructionRootPath,'dat',stimFileName,['mov_hallway_' mosaicFile '.mat']);
load(stimFile);


%%
szStim = size(stim);
% stimTest = stim(:,end-.05*szStim(2)+1:end);
% stimTest = stim(:,end-.25*szStim(2)+1:end);
stimTest = stim;
stimTest = single(stimTest);
% spikeTest = spikeResp(:,end-.05*szStim(2)+1:end);
spikeTest = spikeResp;
stimTestzm = ((stimTest)-(ones(size(stimTest,1),1)*mean(stimTest,1)));
stimTest = stimTestzm;

%%

% for percentSVind = 1:length(percentSVarr)
%     for trainSizeInd = 1:length(trainSizeArr)
%     filterFileFull = fullfile(reconstructionRootPath,'dat',[filterFile '.mat']);
    

rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
filterFile = 'filtersmosaic0_sv50_w1_sh17_dr0.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterMat = data.filterMat; clear data;
% % % % % % % % % % 
    if isempty(pixelWidth)
    filterFileFull  = fullfile(reconstructionRootPath,'/dat/', ...
    [filterFile '.mat']);
        
    else
    filterFileFull  = fullfile(reconstructionRootPath,'/dat/', ...
    [filterFile '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '.mat']);
    end
        load(filterFileFull);
% % % % % % % 


%     shiftval = 9;
%     shiftval = shiftTime+9;
    
    spikeAug(1,:) = ones(1,size(spikeTest,2));
    spikeAug(2:9716+1,:) = spikeTest;%spikeResp;
    % load('filters__mosaic0.mat')
    
%     movRecon = filterMat'*spikeAug;
    mse = 0;
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     nFramesPlay = 40;
%     figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));
    % save('hallwayReconMovie.mat','movRecon');
    %%
    lambda = .01;%.01;%.0075;
    filterMat2 = zeroFilter(filterMat,lambda);
    movRecon2 = filterMat2'*spikeAug;
    
%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon2.avi']
%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon_pros18_filt0025.avi']
%     p2.vname = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/aug20_hallway_recon_pros35_filt0005.avi']
%     p2.save = true;
%     p2.FrameRate = 25;
% 
%     figure; ieMovie(reshape(movRecon2,[100 100 size(movRecon2,2)]),p2);
    figure; ieMovie(reshape(movRecon2,[100 100 size(movRecon2,2)]));
% % % % % % % % % % % %     
%     % works for healthy
%     % figure; fr = 1;
%     fr = fr+10;
%     subplot(121);
%     imagesc(reshape(stim(:,fr),[100 100]));colormap gray;
%     
%     subplot(122);
%     imagesc(reshape(movRecon2(:,fr+16),[100 100])); colormap gray;
%     
%     % for prosthesis
%     figure;
% %     for frskip = 1:18
%     fr = 1;
%     frskip = 7;
% %         figure; fr = 1;
%     fr = fr+10;
%     subplot(121);
%     imagesc(reshape(stim(:,fr),[100 100]));colormap gray;
%     
%     subplot(122);
%     imagesc(reshape(movRecon2(:,fr+frskip),[100 100])); colormap gray;
% %     end
% % % % % % % % % % % 


% %     movReconPlay = reshape(movRecon2,[100 100 size(spikeResp,2)]);
% %     imSingle = testmovieshort(:,:,40); imMean = mean(imSingle(:));
% %     mtest = max(abs(imSingle(:)-imMean));
% %     imSingleRecon = movReconPlay(:,:,40); mtestrecon = max(abs(imSingleRecon(:)));
% %     figure; subplot(131); imagesc((testmovieshort(:,:,40)-imMean)/mtest); colorbar; colormap gray;
% %     
% %     subplot(132); imagesc((movReconPlay(:,:,40)/mtestrecon)); colorbar; colormap gray
% %     subplot(133); imagesc((movReconPlay(:,:,40)/mtestrecon)-(testmovieshort(:,:,40)-imMean)/mtest); colorbar; colormap gray
% 
% %     clear movRecon
%     %
% %     lambda = .01;
% %     filterMat2 = zeroFilter(filterMat,lambda);
% %     movRecon2 = filterMat2'*spikeAug;
% 
% %     shiftval = 5;
%     % shiftval = 3;
% %     mc1 = (movReconPlay(:,:,shiftval+1:szLen+1));
% %     mc2 = (testmovieshort(:,:,1:szLen-shiftval+1));
% 
% %     mc1 = ieScale(movReconPlay(:,:,shiftval+1:szLen+1));
% %     mc2 = ieScale(testmovieshort(:,:,1:szLen-shiftval+1));
% %     clear bigMovie
% %     bigMovie(1:100,101:200,:) = mc1;  bigMovie(1:100,1:100,:) = mc2;
% %     errmov =(mc1-mean(mc1(:)))-(mc2-mean(mc2(:)));
% %     errtot = ((errmov.^2));
% 
% %     mc1rs = RGB2XWFormat(mc1);
% %     mc2rs = RGB2XWFormat(mc2);
%     
% %     mc1rz = mc1rs;% - ones(size(mc1rs,1),1)*mean(mc1rs,1);    
% %     mc2rz = mc2rs;%%%% - ones(size(mc2rs,1),1)*mean(mc2rs,1);
% %     
% % %     std_recon = std(mc1rz(:));
% % %     std_test = std(mc2rz(:));
% % %     
% % %     mc1rz = mc1rz*std_test/std_recon;
%     
%     
% %     errmov =(mc1rz)-(mc2rz);
% %     movRecon = movRecon2;
% 
%     if permuteFlag
%     
%         movReconPlay = reshape(movRecon2,[100 100 size(spikeResp,2)]);
%         movRecon = RGB2XWFormat(permute(movReconPlay,[2 1 3]));
% %         movRecon = RGB2XWFormat((movReconPlay));
%         shiftArr=2;
%     else
%         movRecon = movRecon2;
%     end
%     
%     szLen = size(spikeTest,2)-1;
%     shiftval = shiftArr+1;
%     
%     movReconzm = movRecon - ones(size(movRecon,1),1)*mean(movRecon);
% %     movReconStd = std(movReconzm(:));
% %     stimTestStd = std(stimTest(:));
% %     movReconStd = std(movReconzm);
%     movReconedge = reshape(movReconzm,[100 100 size(movRecon,2)]);
%     movReconStd = std(RGB2XWFormat(movReconedge));
% %     movReconStd = std(RGB2XWFormat(movReconedge(10:90,10:90,:)));
%     stimTestStd = std(stimTest);
%     
% %     movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./std(movReconzm)));
%     
%     movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./movReconStd));
% %     errmov = (movReconStd/stimTestStd)*single(stimTest(:,1:szLen-shiftval+1))-movReconzm(:,shiftval+1:szLen+1);
% %     errmov = single(ones(size(stimTest,1),1)*mean(stimTest(:,1:szLen-shiftval+1))+stimTest(:,1:szLen-shiftval+1))-movRecon(:,shiftval+1:szLen+1);
% 
% %     errmov = single(stimTest(:,1:szLen-shiftval+1))-(stimTestStd/movReconStd)*movReconzm(:,shiftval+1:szLen+1);
%     errmov = single(stimTest(:,1:szLen-shiftval+1))-movReconNorm(:,shiftval+1:szLen+1);
%     
%     errtot = ((errmov.^2));
% 
%     mse(percentSVind,trainSizeInd,:) = sqrt(mean(errtot));
%     
%     
% %     rs1 = reshape((mean(errmov.^2,2)),[100 100]);
% %     sqrt(mean(rs1(:)));
% %     rs2 = rs1(6:end-5,6:end-5);
% %     mse(percentSVind,trainSizeInd,:) = sqrt(mean(rs2(:)))
% 
%    %% 
%        trainSizeMat(percentSVind,trainSizeInd) = trainSizeArr(trainSizeInd);
%     evMat(percentSVind,trainSizeInd) = percentSVarr(percentSVind);
% 
% %     end
% % % mse(percentSVind,:)
% % end
% % rs1 = reshape((mean(errmov.^2,2)),[100 100]);
% % sqrt(mean(rs1(:)))
% % rs2 = rs1(6:end-5,6:end-5);
% % sqrt(mean(rs2(:)))
% 
% %    47.7779 
% %    40.3098
% 
% % 
% %%
% clear movBig;
% shiftval = 3;
% movBig(:,101:200,:) = reshape(stim(:,1:szLen-shiftval+1),[100 100 size(movRecon(:,1:szLen-shiftval+1),2)]);
% % movBig(:,1:100,:) = reshape((stimTestStd/movReconStd)*movReconzm(:,shiftval+1:szLen+1),[100 100 size(movReconzm(:,shiftval+1:szLen+1),2)]);
% movBig(:,1:100,:) = reshape(movReconNorm(:,shiftval+1:szLen+1),[100 100 size(movReconzm(:,shiftval+1:szLen+1),2)]);
% figure; 
% ieMovie(movBig(:,:,1:40));
% % % 
% % shiftval = 4;
% % figure; subplot(121); imagesc(reshape(movReconNorm(:,shiftval+1+350),[100 100])); axis image; subplot(122); imagesc(reshape(stimTest(:,1+350),[100 100])); colormap gray; axis image
% %%    
% % for shiftval = 1:8
% %     movReconzm = movRecon - ones(size(movRecon,1),1)*mean(movRecon);
% % %     movReconStd = std(movReconzm(:));
% % %     stimTestStd = std(stimTest(:));
% %     movReconStd = std(movReconzm);
% %     stimTestStd = std(stimTest);
% %     
% %     movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./std(movReconzm)));
% %     
% % %     errmov = (movReconStd/stimTestStd)*single(stimTest(:,1:szLen-shiftval+1))-movReconzm(:,shiftval+1:szLen+1);
% % %     errmov = single(ones(size(stimTest,1),1)*mean(stimTest(:,1:szLen-shiftval+1))+stimTest(:,1:szLen-shiftval+1))-movRecon(:,shiftval+1:szLen+1);
% % 
% % %     errmov = single(stimTest(:,1:szLen-shiftval+1))-(stimTestStd/movReconStd)*movReconzm(:,shiftval+1:szLen+1);
% %     errmov = single(stimTest(:,1:szLen-shiftval+1))-movReconNorm(:,shiftval+1:szLen+1);
% %     
% %     errtot = ((errmov.^2));
% % 
% %     msesh(shiftval,1:size(errtot,2)) = (mean(errtot));
% %     end
% %     
% %        trainSizeMat(percentSVind,trainSizeInd) = trainSizeArr(trainSizeInd);
% %     evMat(percentSVind,trainSizeInd) = percentSVarr(percentSVind);
% % 
% % %     end
% % figure; plot(mean(msesh,2));
% %%
% 
% % figure; 
% % plot(1e6*trainSizeMat',mse'/255,'-x','linewidth',4)
% % hold on;
% % % plot(1e6*trainSizeMat',mseh','-x','linewidth',4)
% % grid on; % axis([0 1e6 .13 .19]);
% % xlabel('Training Set Size','fontsize',15); 
% % ylabel('Mean Square Error - Reconstruction','fontsize',15);
% % set(gca,'fontsize',15)
% % 
% % legend('25% SVs','50% SVs','75% SVs','location','sw')
% % % legend('10% SVs','20% SVs','30% SVs','40% SVs'
% 
% %%
% % 
% % for ri = 7561:9807
% %     filterSpace = filterMat(ri,:);
% %     filterIm = reshape(filterSpace,[100 100])';
% %     filterMatTrans(ri,:) = filterIm(:);
% % end