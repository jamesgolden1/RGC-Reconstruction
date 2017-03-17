function filterMat = linearReconstructSVD_short_midgets_both(stimName,resp,fileext, windowsize,includedComponentsArray, trainSize, shifttime, stimType)
%Linear reconstruction using least squares estimate for filters
%inputs: stim (movie stimulus reshaped to framewidth*frameheight x time bins)
%        resp (spike response of size numCells x time bins)

%split data into train, validate, and test
sizeStim = size(resp,2);
% sizeStim = 40*12000;%size(stim,2);
% trainSize = 0.6;
validSize = 0.2;
trainSize
%indices to block out train and validation data from the full stimulus
i1 = floor(trainSize*sizeStim);
i2 = floor((trainSize+validSize)*sizeStim);

%number of time bins used in reconstruction from single-timepoint
numbins = windowsize;
% shifttime = 9;
disp('Creating training response matrix...')
%Create the training response matrix (numtimepoints x numcells*numbins)
trainTimes = 1:i1-numbins-shifttime;
respTrain = uint8(zeros(length(trainTimes),size(resp,1)*numbins+1));
respTrain(:,1) = uint8(ones(length(trainTimes),1));
for t = 1:length(trainTimes)
    starttime = trainTimes(t)+1;
    endtime = trainTimes(t)+numbins;
%     size(resp)
%     size(respTrain)
%     numbins
    respTrain(t,2:end) = uint8(reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins));
end

% disp('Creating testing response matrix...')
% %Create the testing response matrix
% testTimes = i1+1:i2-numbins;
% respTest = zeros(length(testTimes),size(resp,1)*numbins+1);
% respTest(:,1) = ones(length(testTimes),1);
% for t = 1:length(testTimes)
%     starttime = testTimes(t);
%     endtime = testTimes(t)+numbins-1;
%     respTest(t,2:end) = reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins);
% end

disp('Calculating pseudoinverse...')
%get pseudoiverse of respose matrix
% corrTrain = pinv(respTrain);
if ~strcmpi(stimType,'wnZero')
    clear resp;
end

% [Utrain, Strain, Vtrain] = svd(respTrain, 'econ');
% wVector =  (V(:,includedComponents) * inv(S(includedComponents,includedComponents)) * (U(:,includedComponents))') * C;

% includedComponentsArray = [100:100:3000];
% includedComponentsArray = [575 584 612 637 663]%[500: 25: 700];
% includedComponentsArray = [212 225 237 250 262 275 288 312 325 337 350 363 ];
% includedComponentsArray = [2600:100:3000];

tic
% load('on_off_covar_svd_strain.mat');
% load('on_off_covar_svd_utrain.mat');
% load('on_off_covar_svd_vtrain.mat');

for icind = 1:length(includedComponentsArray)
    if includedComponentsArray(icind) < 1
       percentSV = includedComponentsArray(icind);
       includedComponentsArray(icind) = round(percentSV*size(respTrain,2));
    end
    includedComponents = [1:includedComponentsArray(icind)];
    
% load('on_off_covar_svd_utrain.mat');
size(respTrain);
tic
[Utrain, Strain, Vtrain] = svd(single(respTrain), 'econ');
toc
% save('svd_ns100_jan_1st3.mat','Utrain', 
% tic
% [Utrain, Strain, Vtrain] = svds(respTrain, includedComponentsArray);
% toc

% tic
% 
% % set options
% opts.tol = 1e-8;
% opts.maxit = 150;
% [Utrain, Strain, Vtrain, lmout] = lmsvd(respTrain, 2000,opts);
% toc
% disp('Calculating training matrix');
% corrTrainSVD =  (Vtrain(:,includedComponents) * inv(Strain(includedComponents,includedComponents)) * (Utrain(:,includedComponents))');


% clear respTrain
szStrain = size(Strain);
StrainInv = zeros(szStrain);
StrainInv(1:szStrain(1)+1:end) = 1./diag(Strain);
corrTrainSVD =  (Vtrain(:,includedComponents) * (StrainInv(includedComponents,includedComponents)) * (Utrain(:,includedComponents))');


clear Utrain
disp('loading stim movie');
% load ../dat/movie_may26.mat
% load ../dat/movie_onMidget.mat

% load('C:\Users\James\Documents\matlab\github\RGC-Reconstruction\dat\movie_spikeResp_onParasol_fast','stim')

% load('C:\Users\James\Documents\matlab\github\RGC-Reconstruction\dat\movie_spikeResp_all0')

if ismac || isunix
    load([reconstructionRootPath '/dat/' stimName]);
else
    load([reconstructionRootPath '\dat\' stimName]);
end

if strcmpi(stimType,'ns')
% Zero mean for NS
for blockNum = 1:floor(size(stim,2)/12000)
    stim(:,(blockNum-1)*12000+1:blockNum*12000) = ...
        uint8(128+127*(double(stim(:,(blockNum-1)*12000+1:blockNum*12000)) - ones(size(stim,1),1)*mean(stim(:,(blockNum-1)*12000+1:blockNum*12000),1)));
end
end

% load ../dat/movie_onMidget_long300.mat


%estimate reconstructed pixels. For memory efficiency, do pixels in batches. 

%%training and testing done in one loop while storing the filters
% recons_train = zeros(size(stim,1),length(trainTimes),'single');
% recons_test = zeros(size(stim,1), length(testTimes),'single');
batchsize = 200; %num pixels you want to do at once
% filterMat = zeros(size(respTrain,2),size(stim,1));
% filterMat = zeros(4321,size(stim,1));
% filterMat = zeros(6751,size(stim,1));
filterMat = zeros(szStrain(1),size(stim,1));
disp(['Reconstructing train and test stimuli with pixel batch size ' num2str(batchsize)])


if strcmpi(stimType,'ns')
    for pix = 1:batchsize:size(stim,1)-(batchsize-1)
        % pix
        % 		[pix size(stim,1)-(batchsize-1) ]
        pixelTC = (1/255)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)')-.5;
%         pixelTC = (1)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)');
        %     pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)')-.5;
        %     filter = corrTrain*pixelTC;
        filter = corrTrainSVD*pixelTC;
        %     recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
        %     recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
        filterMat(:,pix:pix+(batchsize-1)) = filter;
    end
    
    
elseif strcmpi(stimType,'wnZero')
    
    %%%%% Find STA
%     tstim = 2/119.5172;
    STA_length = 30;
%     movie_size = size(fitmovie{1});
    STA = zeros(size(stim,1),size(resp,1));
%     fitframes = movie_size(3);
    stimD = -.5+single(stim);%(:,9+[1:60000]))-ones(10000,1)*mean(stim(:,9+[1:60000]));
    for cellind = 1:size(resp,1)
%         sp_frame = floor(fitspikes{cellind}(:)/tstim);
        sp_rel = find(resp(cellind,:));
        for i = 25%:STA_length
            % STA(:,cellind) = sum(stimD(:,(sp_rel)-STA_length+1+i),2);%(stimD*single(resp(cellind,(sp_rel)-STA_length+1+i))); % (sp_rel)-STA_length+1+i
            STA(:,cellind) = sum(-.5+ single(stimD(:,(sp_rel)-STA_length+1+i)),2);
        end
    end
    % figure; for cellind = 1:64; subplot(8,8,cellind); imagesc(reshape(squeeze(STA(:,cellind)'),[80 40])); colormap parula;  end;

    clear stimD resp

    %%%%% Find max of STA
    clear mgr mgc fmaxr fmaxc
    [mgr,mgc] = meshgrid(1:80,1:40);
    
    [cmax,cind] = max(abs(STA));
    [fmaxr,fmaxc] = ind2sub([80 40],cind);
    
    mgrmat = mgr(:)*ones(1,size(fmaxr,2));
    fmaxrmat = ones(size(mgrmat,1),1)*fmaxr;
    mgrd = ((mgrmat - fmaxrmat)').^2;
    
    mgcmat = mgc(:)*ones(1,size(fmaxc,2));
    fmaxcmat = ones(size(mgcmat,1),1)*fmaxc;
    mgcd = ((mgcmat - fmaxcmat)').^2;
    
    [iv,iord] = sort(mgc(:),'ascend');
    
    dp = sqrt(mgrd(:,iord)+mgcd(:,iord));
    figure; imagesc(reshape(STA(:,40).*(dp(40,:)<12)',[80 40]));
%     figure; imagesc(reshape(STA(:,100),[80 40]).*(reshape(dp(100,:)',[40 80])'<6))
%     filterMat2 = filterMat;
%     filterMat2(dp>5) = 0;
    %%%%%%%
%     trainTimes = 1:i1-numbins-shifttime;
%     respTrain = uint8(zeros(length(trainTimes),size(resp,1)*numbins+1));
%     respTrain(:,1) = uint8(ones(length(trainTimes),1));
%     for t = 1:length(trainTimes)
%         starttime = trainTimes(t)+1;
%         endtime = trainTimes(t)+numbins;
%         %     size(resp)
%         %     size(respTrain)
%         %     numbins
%         respTrain(t,2:end) = uint8(reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins));
%     end
    %%%%%%%
    batchsize2 = 3200;
       pix = 1; pixelTC = (1/255)*double(stim(pix:pix+(batchsize2-1),shifttime+trainTimes)')-.5;
 
    batchsize = size(stim,1);
    zeroRadius = 10;
    filterMat = zeros(size(corrTrainSVD,1),size(stim,1));
%         pixelTC = (1/255)*double(stim')-.5;%(pix:pix+(batchsize-1),shifttime+trainTimes)')-.5;

% % % % % % % % % % % 
%     for cellind = 1:5%:size(dp,1)
%         cellind
%         filter = zeros(1861,3200); clear filterTemp
%     for pix = 1:batchsize:size(stim,1)-(batchsize-1)
%         % pix
% %         % 		[pix size(stim,1)-(batchsize-1) ]
% %         pixelTC = (1/255)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)')-.5;
% %         pixelTC = (1)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)');
% %             pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)')-.5;
%         %     filter = corrTrain*pixelTC;
% %         filter = corrTrainSVD*pixelTC;
% %         filterTemp = corrTrainSVD(2+(windowsize*(cellind-1)):2+(windowsize*(cellind))-1,:)*pixelTC(:,dp(cellind,:)<zeroRadius);
%         filterTemp = corrTrainSVD*pixelTC(:,dp(cellind,:)<zeroRadius);
%         filter(:,find(dp(cellind,:)<zeroRadius)) = filterTemp;
%         %     recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
%         %     recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
% %         filterMat(2+(windowsize*(cellind-1)):2+(windowsize*(cellind))-1,pix:pix+(batchsize-1)) = filter;
%         filterMat(:,pix:pix+(batchsize-1)) = filterMat+filter;
%     end
%     end
% % % % % % % % % % % % 
batchsize = 1;
dpThresh = dp<12;
filterMat = zeros(1861,size(stim,1));
%  for cellind = 1:5%:size(dp,1)
        cellind
        filter = zeros(1861,3200); clear filterTemp
    for pix = 1:batchsize:size(stim,1)-(batchsize-1)
        if ~mod(pix,500)
            pix
        end
        relcells = find(dpThresh(:,pix));
        % pix
%         % 		[pix size(stim,1)-(batchsize-1) ]
%         pixelTC = (1/255)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)')-.5;
%         pixelTC = (1)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)');
%             pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)')-.5;
        %     filter = corrTrain*pixelTC;
%         filter = corrTrainSVD*pixelTC;
%         filterTemp = corrTrainSVD(2+(windowsize*(cellind-1)):2+(windowsize*(cellind))-1,:)*pixelTC(:,dp(cellind,:)<zeroRadius);
        relcellsarr = []; for k = 1:length(relcells); relcellsarr = [relcellsarr (1+(windowsize*(k-1))+[1:10])]; end;
        filter = corrTrainSVD(relcellsarr,:)*pixelTC(:,pix);%(:,dp(cellind,:)<zeroRadius);
        
%         corrTemp(2+(windowsize*(k-1)):2+(windowsize*(k))-1,:) = corrTrainSVD(2+(windowsize*(k-1)):2+(windowsize*(k))-1,:);

%         filter = *pixelTC(:,pix);%(:,dp(cellind,:)<zeroRadius);
%         filter(:,find(dp(cellind,:)<zeroRadius)) = filterTemp;
        %     recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
        %     recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
%         filterMat(2+(windowsize*(cellind-1)):2+(windowsize*(cellind))-1,pix:pix+(batchsize-1)) = filter;
        filterMat(relcellsarr,pix:pix+(batchsize-1)) = filter;%...
%             filterMat(relcellsarr,pix:pix+(batchsize-1))+filter;
    end
%  end
    
    figure; for fr = 1:64; subplot(8,8,fr); imagesc(reshape(filterMat(100+fr,:),[80 40])'); colormap parula;  end;
    
else
    
    for pix = 1:batchsize:size(stim,1)-(batchsize-1)
%         pix
        % 		[pix size(stim,1)-(batchsize-1) ]
        %     pixelTC = (1/255)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)')-.5;
        pixelTC = (1)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)');
        %     pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)')-.5;
        %     filter = corrTrain*pixelTC;
        filter = corrTrainSVD*pixelTC;
        %     recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
        %     recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
        filterMat(:,pix:pix+(batchsize-1)) = filter;
    end

end
clear stim
fileext2 = [fileext '_svd_' num2str(includedComponentsArray(icind)) '_len_' num2str(round(100*trainSize))];

% save(strcat('../output/svd_reconstruct/recons_train_',fileext2),'recons_train','-v7.3');
% save(strcat('../output/svd_reconstruct/recons_test_', fileext2),'recons_test','-v7.3');
% % save(strcat('../output/svd_reconstruct/respTrain_',fileext2), 'respTrain', '-v7.3');
% save(strcat('../output/svd_reconstruct/respTest_',fileext2), 'respTest', '-v7.3');
% save(strcat([reconstructionRootPath '\dat\filters_'],fileext2), 'filterMat','-v7.3');

% save(strcat([reconstructionRootPath '/dat/filters_'],fileext2), 'filterMat','-v7.3');

clear recons_train recons_test stim
toc


end

end

