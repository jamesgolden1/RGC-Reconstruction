function natScenes = generateStimulus(stimTypeBuild, blockNum, nSteps, testFlag)

    if strcmpi(stimTypeBuild,'ns')
        % Training/testing with natural scenes data on server
        if blockNum <= 288
            movsm = parload(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(mod(blockNum-1,12)+1) '.mat']);            
            natScenes = movsm(1:100,1:100,nSteps*floor((blockNum-1)/12)+randperm(nSteps));
            
%             rd = RdtClient('isetbio');
%             rd.crp('/resources/data/reconstruction/training');            
%             fname = ['im_block_' num2str(blockNum)];
%             data  = rd.readArtifact(fname, 'type', 'mat');
%             whiteNoiseSmall = data.whiteNoiseSmall; clear data;
                  
%             natScenes = whiteNoiseSmall;%(1:100,1:100,nSteps*floor((blockNum-1)/12)+randperm(nSteps));
        else
            movsm = parload(['/Volumes/Lab/Users/james/RGC-Reconstruction/dat/imagenetBlocks/movsm_' num2str(12+mod(blockNum-1,12)+1) '.mat']);            
            natScenes = movsm(1:100,1:100,nSteps*floor((blockNum-1)/12)+randperm(nSteps));
            
%             rd = RdtClient('isetbio');
%             rd.crp('/resources/data/reconstruction/training');            
%             fname = ['im_block_' num2str(blockNum)];
%             data  = rd.readArtifact(fname, 'type', 'mat');
%             whiteNoiseSmall = data.whiteNoiseSmall; clear data;
%                   
%             % natScenes = whiteNoiseSmall(1:100,1:100,nSteps*floor((blockNum-1)/12)+randperm(nSteps));
%             
%             % natScenes = movsm(1:100,1:100,nSteps*(floor((-288+blockNum-1)/12))+randperm(nSteps));
%             natScenes = whiteNoiseSmall;%(1:100,1:100,nSteps*(floor((-288+blockNum-1)/12))+randperm(nSteps));
        end
    elseif strcmpi(stimTypeBuild,'wn')
        % Training with white noise
        natScenesRaw = (rand(100,100,nSteps));
        natScenes = 255*round(natScenesRaw); clear natScenesRaw;
    elseif strcmpi(stimTypeBuild,'ns500')
        % Testing with natural scenes data on RDT for demo purposes
        rd = RdtClient('isetbio');
        rd.crp('/resources/data/reconstruction');
        filterFile = 'stim500_1.mat';
        data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
        stim = data.stim; clear data;
        
        %         load([reconstructionRootPath '/dat/nov_results/healthyTest/stim500_1.mat']);
        natScenes = reshape(stim(:,1:500),[100 100 500]);
    end
    
    if testFlag
        
        testInds = (startInd-1)+[1:20:500-20];
        natScenesAll = natScenes;
        natScenes = zeros(size(natScenesAll));
        for ti = 0:19
            natScenes(:,:,testInds+ti) = natScenesAll(:,:,testInds);
        end
        
    end
    
    %     testInds = [1:25:500-20];
    %     natScenesAll = natScenes;
    %     natScenes = zeros(size(natScenesAll));
    %     for ti = 0:24%19s
    %     natScenes(:,:,testInds+ti) = natScenesAll(:,:,testInds);
    %     end
end

