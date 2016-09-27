function [filterMat] = linearReconstructSVD(stim,resp,fileext, windowsize)
%Linear reconstruction using least squares estimate for filters
%inputs: stim (movie stimulus reshaped to framewidth*frameheight x time bins)
%        resp (spike response of size numCells x time bins)

%split data into train, validate, and test
sizeStim = size(stim,2);
% sizeStim = 120000;
trainSize = 1;
validSize = 0.2;

%indices to block out train and validation data from the full stimulus
i1 = floor(trainSize*sizeStim);
i2 = floor((trainSize+validSize)*sizeStim);

%number of time bins used in reconstruction from single-timepoint
numbins = windowsize;

disp('Creating training response matrix...')
%Create the training response matrix (numtimepoints x numcells*numbins)
trainTimes = 1:i1-numbins;
respTrain = zeros(length(trainTimes),size(resp,1)*numbins+1);
respTrain(:,1) = ones(length(trainTimes),1);
for t = 1:length(trainTimes)
    starttime = trainTimes(t)+1;
    endtime = trainTimes(t)+numbins;
    respTrain(t,2:end) = reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins);
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
size(respTrain)
includedComponentsArray = [800];
icind=1;
    includedComponents = [1:includedComponentsArray(icind)];
    
[Utrain, Strain, Vtrain] = svd(respTrain, 'econ');
% wVector =  (V(:,includedComponents) * inv(S(includedComponents,includedComponents)) * (U(:,includedComponents))') * C;

% includedComponentsArray = [100:100:3000];
% includedComponentsArray = [575 584 612 637 663]%[500: 25: 700];
% includedComponentsArray = [212 225 237 250 262 275 288 312 325 337 350 363 ];


tic
% load('on_off_covar_svd_strain.mat');
% load('on_off_covar_svd_utrain.mat');
% load('on_off_covar_svd_vtrain.mat');

for icind = 1:length(includedComponentsArray)
    
    includedComponents = [1:includedComponentsArray(icind)];

% [Utrain, Strain, Vtrain] = svd(respTrain, 'econ');

corrTrainSVD =  (Vtrain(:,includedComponents) * inv(Strain(includedComponents,includedComponents)) * (Utrain(:,includedComponents))');

% corrTrainSVD =  (Vtrain * ((Strain) \ (Utrain)'));
% corrTrainSVD = inv(respTrain'*respTrain);

%estimate reconstructed pixels. For memory efficiency, do pixels in batches. 

%%training and testing done in one loop while storing the filters
recons_train = zeros(size(stim,1),length(trainTimes),'single');
% recons_test = zeros(size(stim,1), length(testTimes),'single');
batchsize = 80; %num pixels you want to do at once
filterMat = zeros(size(respTrain,2),size(stim,1));

disp(['Reconstructing train and test stimuli with pixel batch size ' num2str(batchsize)])

for pix = 1:batchsize:size(stim,1)-(batchsize-1) 
		pix
    pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)');
%     filter = corrTrain*pixelTC;
size(pixelTC)
size(corrTrainSVD)
    filter = corrTrainSVD*pixelTC;
%     recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
%     recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
    filterMat(:,pix:pix+(batchsize-1)) = filter;
end

fileext2 = [fileext '_svd_' num2str(includedComponentsArray(icind))];

% save(strcat('../output/svd_reconstruct/recons_train_',fileext2),'recons_train','-v7.3');
% save(strcat('../output/svd_reconstruct/recons_test_', fileext2),'recons_test','-v7.3');
% % save(strcat('../output/svd_reconstruct/respTrain_',fileext2), 'respTrain', '-v7.3');
% save(strcat('../output/svd_reconstruct/respTest_',fileext2), 'respTest', '-v7.3');
% save(strcat('../output/svd_reconstruct/filters_',fileext2), 'filterMat','-v7.3');
save(strcat('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_',fileext2), 'filterMat','-v7.3');
clear recons_train recons_test filterMat
toc
end

end
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/output/svd_reconstruct/filters_off_parasol_prosthesis_svd_400.mat')
