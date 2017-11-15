function [imRefNorm, imTestNorm, matchShift, mse1, medcc, errmean, ccrec] = normalizeImages(imRef, imTest)
% Takes two sets of images in XW format and sets mean to zero and STD equal
% to the STD of imRef.


%% Set reference images to zero mean

stimTest = single(imRef);
stimTestzm = ((stimTest)-(ones(size(stimTest,1),1)*mean(stimTest,1)));
stimTest = stimTestzm;
movRefStd = std(stimTest);
imRefNorm = stimTest;

%% Set test images to zero mean and scale each image to same STDEV as reference

szLen = size(imRef,2)-1;

% Take mean of each image and subtract it
movTestZeroMean = imTest - ones(size(imTest,1),1)*mean(imTest);

movTestZMXYT = reshape(movTestZeroMean,[100 100 size(movTestZeroMean,2)]);
movTestStd = std(RGB2XWFormat(movTestZMXYT));

imTestNorm =  movTestZeroMean.*(ones(size(movTestZeroMean,1),1)*(movRefStd./movTestStd));
%    movReconNorm =  movReconzm.*(ones(size(movRecon,1),1)*(std(stimTest(:,1:size(movRecon,2)))./movReconStd));
%

%% Find optimal shift value

for shiftval = 0:30
    
    errmov = single(stimTest(:,1:20:szLen-shiftval+0))-imTestNorm(:,shiftval+1:20:szLen+0);
    errtot = ((errmov.^2));
    
    msesh(shiftval+1,1:size(errtot,2)) = (mean(errtot));
    mseshm(shiftval+1) = mean(errtot(:));
end


[mv,mi] = min(mseshm);
matchShift = mi-1;

%% Find mse
errMov = imRefNorm - imTestNorm;

errmean = mean(errMov.^2);
mse1 = sqrt(mean(errmean))/255;

%% Find correlation
tshift = matchShift;
ccrec = zeros(size(imTestNorm,2),1);
recctr = 0;
for ii = 1:size(imTestNorm,2)-tshift
    
    %     for ii = 1:size(movReconNorm,2)-18
    recctr=recctr+1;
    ccrec(recctr) = (single(imTestNorm(:,ii+tshift))\stimTest(:,ii+0));
    
end
medcc = median(ccrec(ccrec<1&ccrec~=0));

end

%     figure; hist(ccrec(ccrec<1&ccrec~=0),40)