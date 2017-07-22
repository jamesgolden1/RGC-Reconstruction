
clear

contrastArr =  [.002 .008 .02 .04 .06 .08];
pool = parpool(length(contrastArr));
for gap = [-2]
    parfor contrastInd = 1:length(contrastArr)

        reconLandolt = recon();
        reconLandolt.buildLandolt('nTrials',200,'contrast',contrastArr(contrastInd),'gap',gap);
    end
end

% for contrast = [.02 .04 .06 .08 .11 .15]
%     for gap = [0 ]
%         reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaLandolt('nTrials',200,'contrast',contrast,'gap',gap);
%     end
% end

% for contrast = [.02 .04 .06 .08 .11 .15]
%     for gap = [2]
%         reconPrimaLandolt = recon();
%         reconPrimaLandolt.buildPrimaLandolt('nTrials',200,'contrast',contrast,'gap',gap);
%     end
% end
%         