function mse1 = t_reconPitchTest()


pixelPitchArr = [70 70/2 70/4 70/8];% 70/16 70/32];
currentDecayArr = [2 4 8 16];

% for pitchInd = [1:length(pixelPitchArr)]
%     for currentDecayInd = 1:length(currentDecayArr)
%     
%     buildFile = 'july6prima/raw/build';
%     respFile = ['july6prima/sp_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_' num2str(currentDecayArr(currentDecayInd))];
%     stimFile = ['july6prima/mov_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_' num2str(currentDecayArr(currentDecayInd))];
%     mosaicFile = '_mosaic0';
%     
%     loadSpikesPitch(buildFile,stimFile,respFile,mosaicFile,pixelPitchArr(pitchInd),currentDecayArr(currentDecayInd))
%     
%     end
% end

%%
if 0
movieIn = loadHallStimulus(55);
for pitchInd = 4%[1:length(pixelPitchArr)]
    for currentDecayInd = 1%1:length(currentDecayArr)
%         if ~((pitchInd == 4) && currentDecayInd ==4)
        [pitchInd currentDecayInd]
        primaParams.pixelWidth = pixelPitchArr(pitchInd)*1e-6; % meters
        primaParams.currentDecay = currentDecayArr(currentDecayInd);
        primaParams.ecc = 1.8;       % deg
        primaParams.fov = 1.7/1;     % deg
        
        primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
        primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
        primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance
        
        primaRecon = primaArray(movieIn,primaParams);
        
        primaRecon.compute(movieIn)
        
        spikeResp = mosaicSpikes(primaRecon.innerRetina);
%         respFile = [reconstructionRootPath '/dat/july12hallway/sp_hallway_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_' num2str(currentDecayArr(currentDecayInd)) '.mat'];
        respFile = [reconstructionRootPath '/dat/july12hallway/sp_hallway_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_p25.mat'];
%         save(respFile,'spikeResp');
%         end
    end
end
end
%%
if 1
stim = RGB2XWFormat(loadHallStimulus(600));
reconHealthy = recon();
for pitchInd = 2%1:2%1:length(pixelPitchArr)
    
    for currentDecayInd = 1%1:length(currentDecayArr)
        
%         pixelWidth = pixelPitchArr(pitchInd);
%         currentDecay = currentDecayArr(currentDecayInd);
%     mosaicFile = '_mosaic0';
%     filename1 = [reconstructionRootPath '/dat/' buildFile '_block_' num2str(blockNum) '_pitch_' sprintf('%2.0f',pixelPitch) '_decay_' num2str(currentDecay) '' mosaicFile '.mat'];
%     filename1 = ['spikeResp_hallway_pitch' sprintf('%d',pitchInd) '_decay' num2str(currentDecayInd) '.mat'];
%     load(filename1);
       respFile = [reconstructionRootPath '/dat/july12hallway/sp_hallway_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_' num2str(currentDecayArr(currentDecayInd)) '.mat'];
%           respFile = [reconstructionRootPath '/dat/july12hallway/sp_hallway_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_0.mat'];
     
    load(respFile);
%     respFile = [reconstructionRootPath '/dat/july6prima/sp_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_' num2str(currentDecayArr(currentDecayInd)) mosaicFile '.mat'];
%     stimFile = [reconstructionRootPath '/dat/july6prima/mov_pitch' sprintf('%2.0f',pixelPitchArr(pitchInd)) 'decay_' num2str(currentDecayArr(currentDecayInd)) mosaicFile '.mat'];
%     load(respFile);
%     load(stimFile);


%     mse(pitchInd,currentDecayInd) = reconHealthy.testCV(stim,spikeResp);
    
onParasolInd = 27*31;
% offParasolInd = onParasolInd+31*35;
% onMidgetInd = offParasolInd+54*62;
%  load('filters_mosaic0_sv50_w1_sh4_tr80.mat')
spikeNoLearn = spikeResp;%(randperm(size(spikeResp,1)),randperm(size(spikeResp,2)));
% spikeNoLearn(onParasolInd+1:offParasolInd,:) = 0;
% spikeNoLearn(onMidgetInd+1:end,:) = 0;

% spikeNoLearn(1:onParasolInd,:) = 0;
% spikeNoLearn(offParasolInd+1:onMidgetInd,:) = 0;


    mse1(pitchInd,currentDecayInd,:) = reconHealthy.testCV(stim,spikeNoLearn,'permuteFlag',1);
    
%     spikeRespSh1 = spikeResp( randperm(size(spikeResp,1)),:);
%     spikeRespSh2 = spikeRespSh1(:,randperm(size(spikeRespSh1,2)));
%     mse(pitchInd,currentDecayInd) = reconHealthy.testCV(stim,spikeRespSh2);
   
    end
end

ph=1;
end
%%
% figure; scatter(squeeze(mse1(2,2,1:250)),squeeze(mse1(3,2,1:250))); hold on; line([0 100],[0 100]); axis equal; grid on;
% both
% 
%    75.3441   81.3350   79.5408   78.0975
%    74.1407   80.4186   79.9132   79.5959
%    75.5987   73.2060   72.8124   75.3700
%    75.6536   75.8101   76.5910   76.4740

% only off

% 
%    88.5405   87.0432   80.4448   79.2160
%    88.0536   81.8464   78.4393   78.7468
%    85.3061   75.6519   69.6691   72.4435
%    84.7291   84.0611   82.1206   81.7846

% only on
%    64.9192
%    63.9589
%    68.3832
%    68.9738
% % 
% ans =
% 
%          0   68.3093   76.1724   76.7640
%          0   75.7647   78.5251   78.2878
%          0   76.1740   81.9469   80.4041
%          0   69.8523   72.3767   72.5750

% on sv .25
% 
%    66.6303   68.8197   75.7738   76.2480
%    66.1709   75.2423   77.9555   77.4959
%    69.6325   76.5788   82.6378   80.9844
%    69.6717   70.4753   72.9874   73.1215



function loadSpikesPitch(buildFile,stimFile,respFile,mosaicFile,pixelPitch,currentDecay)
loadFile = buildFile;

dNames = (dir([reconstructionRootPath '/dat/' loadFile '*block_*' '_pitch_' sprintf('%2.0f',pixelPitch) '_decay_' num2str(currentDecay) mosaicFile '.mat']));

numReps = length(dNames);
filename1 = [reconstructionRootPath '/dat/' loadFile '_block_' num2str(451) '_pitch_' sprintf('%2.0f',pixelPitch) '_decay_' num2str(currentDecay) mosaicFile '.mat'];

matf = matfile(filename1);
szMov = size(matf.whiteNoiseSmall);
blocklength = szMov(3);
stim = zeros(szMov(1)*szMov(2),blocklength*numReps,'uint8');
clear matf
blockNum = 0;

for blockNumInd =[1:length(dNames) ]
    % for blockNumInd =[1:12 21:50]
    blockNum = blockNum+1
    
    matf = matfile(filename1);
    spikesoutsm = matf.spikesoutsm;
    % Spikes in this variable for each block
    spikesout = double(spikesoutsm);
    pointer = (blockNum-1)*blocklength;
    
    for i = 1:blocklength
        blocksize = 10;
        endval = i*blocksize;
        if endval > size(spikesout,2)
            endval = size(spikesout,2);
        end
        startval = (i-1)*blocksize + 1;
        spikeResp(:,pointer+i) = sum(spikesout(:,startval:endval),2);
    end
    clear spikesout
    clear spikesoutsm
    szMov = size(matf.whiteNoiseSmall);
    stimtmp = reshape(matf.whiteNoiseSmall,szMov(1)*szMov(2),blocklength);
    %     stimtmp = uint8(128+double(stimtmp) - ones(size(stimtmp,1),1)*mean(stimtmp,1));
    stim(:,(blockNum-1)*blocklength +1 : blockNum*blocklength) = stimtmp;
    
    % Stimulus here
    % whiteNoiseSmall;
end



save([reconstructionRootPath '/dat/' respFile mosaicFile ],'spikeResp','-v7.3');
save([reconstructionRootPath '/dat/' stimFile mosaicFile],'stim','-v7.3')
