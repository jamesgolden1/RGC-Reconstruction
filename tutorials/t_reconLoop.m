





%%

mosaicFile = '_mosaic0';
movieFile  = 'may26primaSmall/mov'; 
spikesFile = 'may26primaSmall/sp';
buildFile  = 'may26primaSmall/build';

windowSize = 1;

shiftArr = [2 3 ];
percentSVarr =[.50 .25 ]%[.05 .075 .1 .125]% [.2 .4 .6 .8 1];

    for shiftind = 1:length(shiftArr);
for percentSVind = 1:length(percentSVarr)
    
        
        percentSV = percentSVarr(percentSVind);
        shifttime = shiftArr(shiftind);
        filterFile = ['may26primaSmall/filters'  mosaicFile sprintf('_sv%2d',100*percentSV) sprintf('_w%d',windowSize) sprintf('_sh%d',shifttime)];
        
        clear pRecon
        pRecon.buildFile = buildFile;
        pRecon.stimFile = movieFile;
        pRecon.respFile = spikesFile;
        pRecon.filterFile = filterFile;
        pRecon.mosaicFile = mosaicFile;
        pRecon.windowSize = windowSize;
        pRecon.percentSV = percentSV;
        % pRecon.stimType = 'ns';
        
        reconHealthy = recon(pRecon);
        
        reconHealthy.train(pRecon,'shifttime',shifttime);
        
    end
end
