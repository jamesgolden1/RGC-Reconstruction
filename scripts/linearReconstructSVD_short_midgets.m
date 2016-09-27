function linearReconstructSVD_short_midgets(stimNull,resp,fileext, windowsize,includedComponentsArray, trainSize)
%Linear reconstruction using least squares estimate for filters
%inputs: stim (movie stimulus reshaped to framewidth*frameheight x time bins)
%        resp (spike response of size numCells x time bins)

%split data into train, validate, and test
sizeStim = 1.5*360000;%size(stim,2);
% trainSize = 0.6;
validSize = 0.2;

%indices to block out train and validation data from the full stimulus
i1 = floor(trainSize*sizeStim);
i2 = floor((trainSize+validSize)*sizeStim);

%number of time bins used in reconstruction from single-timepoint
numbins = windowsize;

disp('Creating training response matrix...')
%Create the training response matrix (numtimepoints x numcells*numbins)
trainTimes = 1:i1-numbins;
respTrain = single(zeros(length(trainTimes),size(resp,1)*numbins+1));
respTrain(:,1) = single(ones(length(trainTimes),1));
for t = 1:length(trainTimes)
    starttime = trainTimes(t)+1;
    endtime = trainTimes(t)+numbins;
    respTrain(t,2:end) = single(reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins));
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
clear resp;

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
    
    includedComponents = [1:includedComponentsArray(icind)];

% load('on_off_covar_svd_utrain.mat');
tic
[Utrain, Strain, Vtrain] = svd(respTrain, 'econ');
toc
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


clear respTrain
szStrain = size(Strain);
StrainInv = zeros(szStrain);
StrainInv(1:szStrain(1)+1:end) = 1./diag(Strain);
corrTrainSVD =  (Vtrain(:,includedComponents) * (StrainInv(includedComponents,includedComponents)) * (Utrain(:,includedComponents))');


clear Utrain
disp('loading stim movie');
% load ../dat/movie_may26.mat
% load ../dat/movie_onMidget.mat
load ../dat/movie_onMidget_long300.mat


%estimate reconstructed pixels. For memory efficiency, do pixels in batches. 

%%training and testing done in one loop while storing the filters
% recons_train = zeros(size(stim,1),length(trainTimes),'single');
% recons_test = zeros(size(stim,1), length(testTimes),'single');
batchsize = 96; %num pixels you want to do at once
% filterMat = zeros(size(respTrain,2),size(stim,1));
% filterMat = zeros(4321,size(stim,1));
% filterMat = zeros(6751,size(stim,1));
filterMat = zeros(szStrain(1),size(stim,1));
disp(['Reconstructing train and test stimuli with pixel batch size ' num2str(batchsize)])

for pix = 1:batchsize:size(stim,1)-(batchsize-1) 
		[pix size(stim,1)-(batchsize-1) ]
    pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)');
%     filter = corrTrain*pixelTC;
    filter = corrTrainSVD*pixelTC;
%     recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
%     recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
    filterMat(:,pix:pix+(batchsize-1)) = filter;
end
clear stim
fileext2 = [fileext '_svd_' num2str(includedComponentsArray(icind)) '_len_' num2str(round(100*trainSize))];

% save(strcat('../output/svd_reconstruct/recons_train_',fileext2),'recons_train','-v7.3');
% save(strcat('../output/svd_reconstruct/recons_test_', fileext2),'recons_test','-v7.3');
% % save(strcat('../output/svd_reconstruct/respTrain_',fileext2), 'respTrain', '-v7.3');
% save(strcat('../output/svd_reconstruct/respTest_',fileext2), 'respTest', '-v7.3');
save(strcat('../output/svd_reconstruct_shorttrain_midgets/filters_',fileext2), 'filterMat','-v7.3');

clear recons_train recons_test filterMat stim
toc


end

end

