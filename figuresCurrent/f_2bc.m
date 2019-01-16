% f_2bc
%
% Golden, J. R., Erickson-Davis, C., Cottaris, N. P., Parthasarathy, N.,
% Rieke, F., Brainard, D. H., Wandell, B. A., & Chichilnisky, E. J. (2017). Simulation
% of visual perception and learning with a retinal prosthesis. bioRxiv,
% 206409.
%
% Fig 2bc RGC spatial RF, generator function
%
% Generates a small RGC mosaic with four types of cells and illustrates
% their temporal impulse responses.
%
% 2017 JRG (c) isetbio team

%% Generate impulse response to an arbitrary stimulus
% In order to find the impulse response at the output of the RGCs, it is
% necessary to run an impulse stimulus through the entire pipeline.

% Compute optical image (oi) response for arbitrary image with a number of time points
sampleRate = .008;
oi = oiCreate;
oi = oiCompute(oi,sceneCreate('rings rays'));
cm = coneMosaic(oi);
cm.integrationTime = sampleRate;
cm.emGenSequence(250);
cm.compute(oi);

% Set impulse and scale up magnitude
cm.absorptions(:,:,1) = 1000*cm.absorptions(:,:,10);
% cm.absorptions(:,:,11:end) = 0;
cm.absorptions(:,:,2:end) = 0;
cm.computeCurrent();

% Choose an RGC cell to observe impulse response
spLocR = 25; spLocC = 25;
bpL = bipolarLayer(cm);

%% Temporal impulse response, 2b
figure;
typeString{1} = 'on diffuse';
typeString{2} = 'off diffuse';
typeString{3} = 'on midget';
typeString{4} = 'off midget';

signMult = [-1 1 -1 1];

for mosaicNumber = 1:4
    
    bpL.mosaic{1} = bipolarMosaic(cm,typeString{mosaicNumber});
    % bpL.mosaic{1} = bipolarMosaic(cm,'off midget');
    bpL.mosaic{1}.compute();
    rgcGain = (bpL.mosaic{1}.responseCenter(spLocR,spLocC,:)-bpL.mosaic{1}.responseSurround(spLocR,spLocC,:));
    rgcGain = rgcGain-0*rgcGain(1,1,end);
    
    subplot(2,2,mosaicNumber);
    plot(squeeze(signMult(mosaicNumber)*rgcGain))
    axis([0 50 -80 max(signMult(mosaicNumber)*rgcGain)]);
    grid on;
    
end

cellType = {'on parasol','off parasol','on midget','off midget'};%,'onsbc'};
innerRetina = rgcLayer(bpL);
innerRetina.mosaic{1} = rgcGLM(innerRetina,bpL.mosaic{1},cellType{1});

%%

% line([100 timePts],[0 0],'color','k'); line([0 0],[-.2 1.2],'color','k');
% % axis off;
% hold on;
% xtick = -[-4000:500:0]; 
% for ti = 1:length(xtick)
% line([xtick(ti) xtick(ti)],[-.1 .1]/2,'color','k');
% end
% 
% ytick = [0:.2:1.2]; 
% for ti = 1:length(ytick)
% line([-200 200]/2,[ytick(ti) ytick(ti)],'color','k');
% end

%% Generator function, 2c

mosaicNumber = 1;
endVal = 4; dt = .01;

nonlinearFunction = innerRetina.mosaic{1}.generatorFunction;

figure; plot(nonlinearFunction(0:dt:endVal),'linewidth',3);
line([0 endVal/dt],[0 0],'color','k','linewidth',2); line([0 0],[0 exp(endVal)],'color','k','linewidth',2);
axis off;
hold on;
xtick = [0:50:400]; 
for ti = 1:length(xtick)
line([xtick(ti) xtick(ti)],[0 10]/4,'color','k');
end

ytick = [0:10:50]; 
for ti = 1:length(ytick)
line([0 50]/4,[ytick(ti) ytick(ti)],'color','k');
end
