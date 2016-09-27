

% function ieMovie(movieMatrix)
% Show a movie of an (x,y,t) matrix with minimal fuss.

%% Show test movie

movieMatrix = movrecons_on;
movieMatrix2 = 0;%movrecons_off;
movieMatrix3 = 0;%movrecons_on_off;

szMovie = size(movieMatrix);

h1=figure; % set(gcf,'position',[160 60 1070 740]);
% vcNewGraphWin([],'upperleftbig'); 

% set(gcf,'position',[463        1078        1097         260]);
cmin = min(movieMatrix(:));
cmax = max(movieMatrix(:));


cmin2 = min(movieMatrix2(:));
cmax2 = max(movieMatrix2(:));


cmin3 = min(movieMatrix3(:));
cmax3 = max(movieMatrix3(:));

% name_str = ['gratingH_20Hz_width_' num2str(params.barWidth) '_freq_' num2str(freqL) '_onM_8_hz_ON_IMS1_' num2str(cputime*100) '.mp4'];
% name_str = 'moving_bar_all_svd.mp4';
name_str = 'gratings_on_midget_all_SVs.mp4';
path_str = ['/Users/james/Documents/MATLAB/RGC-Reconstruction/output/'];
    
vObj = VideoWriter([path_str name_str],'MPEG-4');
vObj.FrameRate = 30;
vObj.Quality = 100;
open(vObj);

hold on;
for frame1 = 1:szMovie(3)
% subplot(1,3,1);
        imagesc(uint8(255*movieMatrix(3:end-2,:,frame1)./cmax));
        
title('on only','fontsize',12);
    colormap gray; 
% subplot(1,3,2);

%         imagesc(uint8(255*movieMatrix2(3:end-2,:,frame1)./cmax2));
%         title('off only','fontsize',12);
%     colormap gray; 
%     subplot(1,3,3);
% 
%         imagesc(uint8(255*movieMatrix3(3:end-2,:,frame1)./cmax3));
%         title('on and off','fontsize',12);
%     colormap gray; 

    drawnow;
    
    F = getframe(h1);
    writeVideo(vObj,F);
end
close(vObj);
% close;