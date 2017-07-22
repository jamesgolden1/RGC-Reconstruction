clear

% contrastArr =  .05*.5*[0 .00002 .0002 .0003:.00002:.0004 .001 .002 .005 .008 .02 .04 .06 .08];
% 
% % contrastArr =  .05*.5*[.0003:.00002:.0004];
% % 
% % contrastArr = .5*[ .0025 .003 .004 ];
% % contrastArr =  .5*[0  .008 .04 .08];
% % contrastArr = .5* [.002  .02 .06 .1];
% pool = parpool(length(contrastArr));

freqArr = [.05 .1 .2 .5 1 2 4 5 8 10 16];
% pool = parpool(length(freqArr));
    
for freqInd = 7%:length(freqArr)

% parfor contrastInd = 1:length(contrastArr)    
    reconGratings = recon();
%     reconGratings.buildGratings('nTrials',100,'contrast',contrastArr(contrastInd),'gratingSpFreq',10);

     reconGratings.buildGratings('nTrials',1,'contrast',1,'gratingSpFreq',freqArr(freqInd),'horizontalFlag',true);
end
% delete(gcp);
% contrastArr =  [0 .1]
