function [filterMat,keepIndicesSort] = linearReconstructSVD_short_midgets_both(stimName,resp,fileext, windowsize,includedComponentsArray, trainSize, shifttime, stimType, dropout)
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
trainTimes = shifttime+1:i1-numbins-shifttime;
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
clear resp
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
%     clear resp;
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
        includedComponentsArray(icind) = round((1-dropout)*percentSV*size(respTrain,2));
    end
    includedComponents = [1:includedComponentsArray(icind)];
    
% load('on_off_covar_svd_utrain.mat');
size(respTrain);

if dropout ~= 0
%     dropoutIndices = 1+round((size(respTrain,1)-1)*rand(ceil(dropout*(size(respTrain,1)-1)),1));
    fullIndices = [1:size(respTrain,2)];
    permIndices = fullIndices(randperm(size(respTrain,2)-1));
    dropoutIndices = permIndices(1:ceil(dropout*(size(respTrain,2)-1)));
    keepIndices = permIndices(1+ceil(dropout*(size(respTrain,2)-1)):end);
    dropoutIndicesSort = sort(dropoutIndices,'ascend'); keepIndicesSort = sort(keepIndices,'ascend');
    % [dropoutIndicesSort(1:10)' keepIndicesSort(1:10)'];
    
    respTrain(:,dropoutIndicesSort) = 0;
    disp('dropout executed');
else
    dropoutIndices = [];
    keepIndicesSort = [1:size(respTrain,2)];
end


% if exist([fileext  '_svd.mat'],'file') ~= 2
    
    tic
    [Utrain, Strain, Vtrain] = svd(single(respTrain(:,keepIndicesSort)), 'econ');
    toc
%     disp([fileext  '_svd_dropout80.mat'])
%     save(['spike_svd.mat'],'Utrain','Strain','Vtrain','-v7.3');
    
% else
%     
%     load([fileext  '_svd.mat']);
% end

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

% corrTrainSVD = zeros(length(keepIndicesSort)+length(dropoutIndicesSort),size(corrTrainSVDred,2));
% corrTrainSVD(1:keepIndicesSort,:) = corrTrainSVDred;
% corrTrainSVD(1:dropoutIndicesSort,:) = zeros(size(corrTrainSVD(1:dropoutIndicesSort)));

clear Utrain
disp('loading stim movie');
% load ../dat/movie_may26.mat
% load ../dat/movie_onMidget.mat

% load('C:\Users\James\Documents\matlab\github\RGC-Reconstruction\dat\movie_spikeResp_onParasol_fast','stim')

% load('C:\Users\James\Documents\matlab\github\RGC-Reconstruction\dat\movie_spikeResp_all0')


% load([reconstructionRootPath '/dat/' stimName]);
if ismac || isunix
    load([reconstructionRootPath '/dat/' stimName]);
    if abs(mean(stim(:,3))) > .01
        stimzm = (single(stim)-(ones(size(stim,2),1)*mean(stim,2)')');
        clear stim; stim = stimzm;
    end
%     load('/Volumes/Lab/Users/james/RGC-Reconstruction/dat/ns100_r2_10/ns100_jan1_mov3_mosaicAll_1246640.mat');
else
    load([reconstructionRootPath '\dat\' stimName]);
end



%estimate reconstructed pixels. For memory efficiency, do pixels in batches. 

%%training and testing done in one loop while storing the filters
% recons_train = zeros(size(stim,1),length(trainTimes),'single');
% recons_test = zeros(size(stim,1), length(testTimes),'single');
batchsize = 200; %num pixels you want to do at once
filterMat = zeros(szStrain(1),size(stim,1));
disp(['Reconstructing train and test stimuli with pixel batch size ' num2str(batchsize)])


for pix = 1:batchsize:size(stim,1)-(batchsize-1)
    % pix
    % 		[pix size(stim,1)-(batchsize-1) ]
    pixelTC = double(stim(pix:pix+(batchsize-1),-shifttime+trainTimes)')-0;
    %         pixelTC = (1)*double(stim(pix:pix+(batchsize-1),shifttime+trainTimes)');
    %     pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)')-.5;
    %     filter = corrTrain*pixelTC;
    filter = corrTrainSVD*pixelTC;
    %     recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
    %     recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
    filterMat(:,pix:pix+(batchsize-1)) = filter;
end

    
clear stim
% fileext2 = [fileext '_svd_' num2str(includedComponentsArray(icind)) '_len_' num2str(round(100*trainSize))];

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

