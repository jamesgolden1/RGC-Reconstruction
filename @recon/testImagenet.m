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

if isempty(regexpi(stimFileName,'hall'))
spikeFile = fullfile(reconstructionRootPath,'dat',stimFileName,['sp_' mosaicFile '.mat']);
load(spikeFile);

stimFile = fullfile(reconstructionRootPath,'dat',stimFileName,['mov_' mosaicFile '.mat']);
load(stimFile);
else
spikeFile = fullfile(reconstructionRootPath,'dat',stimFileName,['sp_hallway_' mosaicFile '.mat']);
load(spikeFile);

stimFile = fullfile(reconstructionRootPath,'dat',stimFileName,['mov_hallway_' mosaicFile '.mat']);
load(stimFile);
end

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
    
% old isetbio commit 2b8fe22
% rd = RdtClient('isetbio');
% rd.crp('/resources/data/istim');
% filterFile = 'filtersmosaic0_sv50_w1_sh17_dr0.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% filterMat = data.filterMat; clear data;


% up to date w master
% rd = RdtClient('isetbio');
% rd.crp('/resources/data/reconstruction');
% filterFile = 'filtersmosaic0_sv50_w1_sh15_dr0_aug27.mat';

% filterFile = 'filtersmosaic0_sv50_w1_sh17_dr0.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% filterMat = data.filterMat; clear data;
% % % % % % % % % % % 
    if isempty(pixelWidth)
    filterFileFull  = fullfile(reconstructionRootPath,'/dat/', ...
    [filterFile '.mat']);
        
    else
    filterFileFull  = fullfile(reconstructionRootPath,'/dat/', ...
    [filterFile '_pitch_' sprintf('%2.0f',pixelWidth) '_decay_' num2str(currentDecay) '.mat']);
    end
        load(filterFileFull);

%     shiftval = 9;
%     shiftval = shiftTime+9;
    
    spikeAug(1,:) = ones(1,size(spikeTest,2));
    spikeAug(2:9716+1,:) = spikeTest;%spikeResp;
    % load('filters__mosaic0.mat')
    
%     movRecon = filterMat'*spikeAug;%(:,randperm(size(spikeAug,2)));
    mse = 0;
%     movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
%     nFramesPlay = 40;
%     figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));
    % save('hallwayReconMovie.mat','movRecon');
    %%
    lambda = .001;%.0075;% .01;
%     lambda = .06;%.001;%1;%.01;%.0075;
    filterMat2 = zeroFilter(filterMat,lambda);
%     movRecon2 = filterMat2'*spikeAug;%(:,randperm(size(spikeAug,2)));;
    
   movRecon2 = filterMat2'*(spikeAug);
   
% %    movRecon2 = filterMat'*(spikeAug);
%     spikeOn = zeros(size(spikeAug));
%    spikeOn(1:1+28*32,:) = spikeAug(1:1+28*32,:);
%    spikeOn(28*32+31*35+1:28*32+31*35+55*63+1,:) = spikeAug(28*32+31*35+1:28*32+31*35+55*63+1,:);
%     filterMat2(1,:) = filterMat(1,:);
% %     spikeOff = spikeAug-spikeOn;
%    movRecon2 = filterMat2'*(spikeOn);
   
%     figure; ieMovie(reshape(movRecon2,[100 100 size(movRecon2,2)]));
%     mse = movRecon2;
% % % % % % % % % % % % %     
% %     % works for healthy
% %     % figure; fr = 1;
% %     fr = fr+25;
%     subplot(121);
%     imagesc(reshape(stim(:,fr),[100 100]));colormap gray;
%     
%     subplot(122);
%     %healthy
%     %imagesc(reshape(movRecon2(:,fr+18),[100 100])); colormap gray;
%     %pros
%      imagesc(reshape(movRecon2(:,fr+8),[100 100])); colormap gray;

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
    movRecon = movRecon2;
    szLen = size(spikeTest,2)-1;
%     shiftval = shiftArr+1;
%     
    movReconzm = movRecon - ones(size(movRecon,1),1)*mean(movRecon);
%     movReconStd = std(movReconzm(:));
%     stimTestStd = std(stimTest(:));
%     movReconStd = std(movReconzm);
    movReconedge = reshape(movReconzm,[100 100 size(movRecon,2)]);
    movReconStd = std(RGB2XWFormat(movReconedge));
%     movReconStd = std(RGB2XWFormat(movReconedge(10:90,10:90,:)));
    stimTestStd = std(stimTest);
    
%     movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./std(movReconzm)));
    
    movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./movReconStd));
%     errmov = (movReconStd/stimTestStd)*single(stimTest(:,1:szLen-shiftval+1))-movReconzm(:,shiftval+1:szLen+1);
%     errmov = single(ones(size(stimTest,1),1)*mean(stimTest(:,1:szLen-shiftval+1))+stimTest(:,1:szLen-shiftval+1))-movRecon(:,shiftval+1:szLen+1);

%     errmov = single(stimTest(:,1:szLen-shiftval+1))-(stimTestStd/movReconStd)*movReconzm(:,shiftval+1:szLen+1);
%     errmov = single(stimTest(:,1:szLen-shiftval+1))-movReconNorm(:,shiftval+1:szLen+1);
    
%     errtot = ((errmov.^2));

%     mse(percentSVind,trainSizeInd,:) = sqrt(mean(errtot));
    
%     errMov = single(stimTest(:,1:25:end-25)) - movReconNorm(:,18+[1:25:end-25]);
%     errMov = single(stimTest(:,1:25:end-25)) - movReconNorm(:,13+[1:25:end-25]);
%        errMov = single(stimTest(:,1:20:end-25)) - movReconNorm(:,13+[1:20:end-25]);
  errMov = single(stimTest(:,1:20:end-25)) - movReconNorm(:,18+[1:20:end-25]);
% errMov = single(stimTest(:,1:end-25)) - movReconNorm(:,15+[1:end-25]);
% errMov = single(stimTest(:,1:end-25)) - movReconNorm(:,3+[1:end-25]);

%        errMov = single(stimTest(:,1:100:end-25)) - movReconNorm(:,8+[1:100:end-25]);
%        errMov = single(stimTest(:,1:25:end-25)) - movReconNorm(:,8+[1:25:end-25]);
%          errMov = single(stimTest(:,1:end-25)) - movReconNorm(:,3+[1:end-25]);
 
         
    errmean = mean(errMov.^2);
    mse = sqrt(mean(errmean)); mse/255
    
    for tshift = 15%18%8:20
    
    ccrec = zeros(2700,1);
    recctr = 0;
%     for ii = 1:20:size(movReconNorm,2)
        
    for ii = 1:size(movReconNorm,2)-18
        recctr=recctr+1;
    ccrec(recctr) = (single(movReconNorm(:,ii+tshift))\stimTest(:,ii+0));
%     [rho,p1]=corr(single(movReconNorm(:,ii)),stimTest(:,ii),'type','spearman');
% mld2 = robustfit(single(movReconNorm(:,ii)),stimTest(:,ii),'welsch',.5);
% ccrec2(recctr) = rho;
% ccrec2(recctr) =  mld2(2);
    
    end
    medcc = median(ccrec(ccrec<1&ccrec~=0))
    
    end
    
    figure; hist(ccrec(ccrec<1&ccrec~=0),40)
    
%     shuff
% .114
% .1168
    % only on 
%     0.6368

% learning
% 0.6546
% 0.6512

% healthy
% 0.8325
    
stimTest2 = stimTest(:,0+[1:20:size(movReconNorm,2)]);
movReconNorm2 = single(movReconNorm(:,18+[1:20:size(movReconNorm,2)]));

% save('recon_prima70_learning_oct8.mat','movReconNorm2','stimTest2','ccrec','-v7.3');

    ii = 10;
figure; 
hold on;

hold on;
dscatter(single(movReconNorm(:,13+ii)),stimTest(:,ii));
dscatter(single(movReconNorm(:,13+ii)),stimTest(:,ii),'plottype','contour',...
    'filled',1,'bins',4*[64 64],'msize',40,'smoothing',20);
figure; imagesc(reshape(stimTest(:,1+289*20),[100 100])); colormap gray; axis image;
    
figure; imagesc(reshape(movReconNorm(:,1+289*20+13),[100 100])); colormap gray; axis image;
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
% for shiftval = 1:8
%     movReconzm = movRecon - ones(size(movRecon,1),1)*mean(movRecon);
% %     movReconStd = std(movReconzm(:));
% %     stimTestStd = std(stimTest(:));
%     movReconStd = std(movReconzm);
%     stimTestStd = std(stimTest);
%     
%     movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./std(movReconzm)));
%     
% %     errmov = (movReconStd/stimTestStd)*single(stimTest(:,1:szLen-shiftval+1))-movReconzm(:,shiftval+1:szLen+1);
% %     errmov = single(ones(size(stimTest,1),1)*mean(stimTest(:,1:szLen-shiftval+1))+stimTest(:,1:szLen-shiftval+1))-movRecon(:,shiftval+1:szLen+1);
%%
% save('errmean70_big.mat','errmean');
for shiftval = 1:20
    
    
%     movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./movReconStd(:,shiftval+1:20:szLen+1)));
    

%     stdratio = (std(stimTest(:,1:25:szLen-shiftval+1)))./movReconStd(:,shiftval+1:25:szLen-shiftval+1);
% %     stdratio(1)
%     stdrs = (ones(size(movRecon(:,1:25:szLen-shiftval+1),1),1))*stdratio;
%     movReconNorm =  movReconzm(:,shiftval+1:25:szLen-shiftval+1).*stdrs;
    
%     errmov = single(stimTest(:,1:szLen-shiftval+1))-(stimTestStd/movReconStd)*movReconzm(:,shiftval+1:szLen+1);

%     errmov = single(stimTest(:,1:25:szLen-shiftval+1))-movReconNorm(:,shiftval+1:25:szLen+1);
 errmov = single(stimTest(:,1:szLen-shiftval+1))-movReconNorm(:,shiftval+1:szLen+1);
%
%        errmov = single(stimTest(:,1:20:szLen-shiftval+1))-movReconNorm(:,shiftval+1:20:szLen+1);
%
    %         errmov = single(stimTest(:,1:20:szLen-shiftval+1))-movReconNorm(:,shiftval+1:20:szLen+1);

%         errmov = single(stimTest(:,1:szLen-shiftval+1))-movReconNorm(:,shiftval+1:szLen+1);
    
    errtot = ((errmov.^2));

    msesh(shiftval,1:size(errtot,2)) = (mean(errtot));
    mseshm(shiftval) = mean(errtot(:));
end
    figure; plot(mseshm);
    [mv,mi] = min(mseshm)
    errmeanout = msesh(mi,:);
    
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
