





%%

mosaicFile = '_mosaic0';
movieFile  = 'may22/mov'; 
spikesFile = 'may22/sp';
buildFile  = 'may22/raw/build';

windowSize = 1;

shiftArr = [4 ];
percentSVarr =[ .75 .50 .25 ]%[.05 .075 .1 .125]% [.2 .4 .6 .8 1];
% percentSVarr = [.625 .375];
% percentSVarr = [.375];
trainSizeArr = [.2 .4 .6 .8];
% trainSizeArr =.1+ [0 .2 .4 .6 .8];
% trainSizeArr = .1:.1:.9
    for shiftind = 1:length(shiftArr);
for percentSVind = 1:length(percentSVarr)
    for trainSizeInd = 1:length(trainSizeArr)
        
        percentSV = percentSVarr(percentSVind);
        shifttime = shiftArr(shiftind);
        trainSize = trainSizeArr(trainSizeInd);
        filterFile = ['may22/filters_wmean/filters'  mosaicFile sprintf('_sv%2.0f',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime) sprintf('_tr%2.0f',100*trainSize)];
        
        clear pRecon
        pRecon.buildFile = buildFile;
        pRecon.stimFile = movieFile;
        pRecon.respFile = spikesFile;
        pRecon.filterFile = filterFile;
        pRecon.mosaicFile = mosaicFile;
        pRecon.windowSize = windowSize;
        pRecon.percentSV = percentSV;
        pRecon.trainFraction = trainSize;
        % pRecon.stimType = 'ns';
        
        reconHealthy = recon(pRecon);
        
        reconHealthy.train(pRecon,'shifttime',shifttime);
    end
    end
end
