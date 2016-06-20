function [ recons_test, filterMat] = linearReconstruct2(stim,resp, pixCells,fileext)
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
numbins = 30;

%Create the training response matrix (numtimepoints x numcells*numbins)
trainTimes = 1:i1-numbins;
respTrain = zeros(length(trainTimes),size(resp,1)*numbins+1);
respTrain(:,1) = ones(length(trainTimes),1);
for t = 1:length(trainTimes)
    starttime = trainTimes(t);
    endtime = trainTimes(t)+numbins-1;
    respTrain(t,2:end) = reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins);
end

%Create the testing response matrix
testTimes = i1+1:i2-numbins;
respTest = zeros(length(testTimes),size(resp,1)*numbins+1);
respTest(:,1) = ones(length(testTimes),1);
for t = 1:length(testTimes)
    starttime = testTimes(t);
    endtime = testTimes(t)+numbins-1;
    respTest(t,2:end) = reshape(resp(:,starttime:endtime)',1,size(resp,1)*numbins);
end
%get pseudoiverse of respose matrix
%corrTrain = inv(respTrain'*respTrain)*respTrain';
clear resp;
%respTrain = sparse(respTrain);
%respTest = sparse(respTest);
%estimate reconstructed pixels. For efficiency, do pixels in batches. 



%%training and testing done in one loop while storing the filters
recons_train = zeros(size(stim,1),length(trainTimes),'single');
recons_test = zeros(size(stim,1), length(testTimes),'single');
batchsize = 1; %num pixels you want to do at once
filterMat = cell(size(stim,1));
for pix = 1:batchsize:size(stim,1)-(batchsize-1)
    pix
    pixelTC = double(stim(pix:pix+(batchsize-1),trainTimes)');
    startid = (pixCells{pix}-1)*numbins+2;
    endid = pixCells{pix}*numbins+1;
    values = arrayfun(@colon, startid, endid, 'Uniform', false);
    coltake = cell2mat(values);
    respT = respTrain(:,coltake);
    corrTrain = inv(respT'*respT)*respT';
    filter = corrTrain*pixelTC;
    recons_train(pix:pix+(batchsize-1),:) = (respT*filter)';%(full(respTrain*filter))';
    recons_test(pix:pix+(batchsize-1),:) = (respTest(:,coltake)*filter)';%(full(respTest*filter))';
    filterMat{pix} = filter;
end
save(strcat('../output/recons_trainRF_',fileext),'recons_train','-v7.3');
save(strcat('../output/recons_testRF_', fileext),'recons_test','-v7.3');
save(strcat('../output/respTrainRF_',fileext), 'respTrain', '-v7.3');
save(strcat('../output/respTestRF_',fileext), 'respTest', '-v7.3');
save(strcat('../output/filtersRF_',fileext), 'filterMat','-v7.3');


end

