function [ recons_test, filterMat] = linearReconstruct(stim,resp,fileext, windowsize)
%Linear reconstruction using least squares estimate for filters
%inputs: stim (movie stimulus reshaped to framewidth*frameheight x time bins)
%        resp (spike response of size numCells x time bins)

%split data into train, validate, and test
sizeStim = size(stim,2);
trainSize = 0.6;
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

disp('Creating testing response matrix...')
%Create the testing response matrix
testTimes = i1+1:i2-numbins;
respTest = zeros(length(testTimes),size(resp,1)*numbins+1);
respTest(:,1) = ones(length(testTimes),1);
for t = 1:length(testTimes)
    starttime = testTimes(t);
    endtime = testTimes(t)+numbins-1;
    respTest(t,2:end) = reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins);
end

disp('Calculating pseudoinverse...')
%get pseudoiverse of respose matrix
corrTrain = pinv(respTrain);
clear resp;


%estimate reconstructed pixels. For memory efficiency, do pixels in batches. 

%%training and testing done in one loop while storing the filters
recons_train = zeros(size(stim,1),length(trainTimes),'single');
recons_test = zeros(size(stim,1), length(testTimes),'single');
batchsize = 96; %num pixels you want to do at once
filterMat = zeros(size(respTrain,2),size(stim,1));

disp(['Reconstructing train and test stimuli with pixel batch size ' num2str(batchsize)])

for pix = 1:batchsize:size(stim,1)-(batchsize-1) 
		pix
    pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)');
    filter = corrTrain*pixelTC;
    recons_train(pix:pix+(batchsize-1),:) = (respTrain*filter)';
    recons_test(pix:pix+(batchsize-1),:) = (respTest*filter)';
    filterMat(:,pix:pix+(batchsize-1)) = filter;
end
% save(strcat('../output/recons_train_',fileext),'recons_train','-v7.3');
% save(strcat('../output/recons_test_', fileext),'recons_test','-v7.3');
% save(strcat('../output/respTrain_',fileext), 'respTrain', '-v7.3');
% save(strcat('../output/respTest_',fileext), 'respTest', '-v7.3');
% save(strcat('../output/filters_',fileext), 'filterMat','-v7.3');


end

