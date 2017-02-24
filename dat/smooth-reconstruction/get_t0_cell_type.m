clear;
addpath(genpath('/Volumes/Lab/Users/james/matlab'));
javaaddpath /home/vision/Nishal/Java/Java/vision7/bin/
% datarun = load_data('2016-02-17-6/data025');
% datarun = load_data('2016-02-17-6/data025-cf/edited/data025-cf/data025-cf');
datarun = load_data('2016-02-17-6/data026');
datarun=load_params(datarun);
% 


% RGB-16-2-0.48-22222
cell_str = {'on parasol','off parasol','on midget','off midget','on smooth','off smooth','all','all5','all4'};

for type_number = 1% 1:5

cell_type = cell_str{type_number}

indices = [];
if strcmp(cell_type,'all')
    for cell_type_number = 1:6
        indices = [indices get_cell_indices(datarun, cell_str
        {cell_type_number})];
    end
elseif strcmp(cell_type,'all5')
    for cell_type_number = 1:5
        indices = [indices get_cell_indices(datarun, cell_str{cell_type_number})];
    end
elseif strcmp(cell_type,'all4')
    for cell_type_number = 1:4
        indices = [indices get_cell_indices(datarun, cell_str{cell_type_number})];
    end
else
    indices = get_cell_indices(datarun,cell_type);
end
cell_ids = datarun.cell_ids(indices)
% indices = get_cell_indices(datarun,'on midget');
% cell_ids = datarun.cell_ids(indices);
% d_save_str = '/Volumes/Lab/Users/james/fits-2016-04-21-10/on_midget';
% d_save_str = '/Volumes/Lab/Users/james/fits-2015-11-09-3/on_parasol';
d_save_str = '/Volumes/Lab/Users/james/smooth-reconstruction/';
% fittedGLM = glm_fit_from_WN([7], '2009-04-3/data008', 'BW-1-4-11111')
cellind= cell_ids;
% for indind =59:length(cellind)

% fitspikesmat = zeros(length(cellind),60000);
% 
% for indind =[1:4 6:length(cellind)]
% [indind length(cellind)]
% % fittedGLM = glm_fit_from_WN(cellind(indind), '2016-02-17-6/data025', 'RGB-16-2-0.48-22222','d_save',d_save_str)
% [fitspikes,fitmovie] = get_spikes_stim(cellind(indind), '2016-02-17-6/data025', 'RGB-16-2-0.48-22222','d_save',d_save_str);
% % fitspikesmat(indind,round(60*fitspikes{indind}))=1+fitspikesmat(indind,round(60*fitspikes{indind}));
% end

indind =[1:length(cellind)];
[fitspikes,fitmovie] = get_spikes_stim(cellind(indind), '2016-02-17-6/data025', 'RGB-16-2-0.48-22222','d_save',d_save_str);

% IS FITMOVIE SAME FOR EACH CELL
% 
% fitmoviemat = reshape(fitmovie,20*40,size(fitmovie,3));
% figure;
% subplot(121);
% imagesc(reshape(fitmoviemat(:,1),40,20))
% subplot(122);
% % figure; 
% imagesc(squeeze(fitmovie(:,:,1)))

% fitmoviemean = ones(800,1)*mean(fitmoviemat);

% figure;
% for ii = 1%:48
%     clear evIdx
%     evIdx = round(60*fitspikes);
%     rf0 = zeros(20*40,59);
%     mnstim =  mean(fitmoviemat(:));
%     for evind = 20:length(evIdx)-20
%         rf0 = rf0 + ((fitmoviemat(:,evIdx(evind)-29:evIdx(evind)+29)));%-((fitmoviemean(:,evIdx(evind)-29:evIdx(evind)+29)));
%     end
%     
%     %     rf2 = RGB2XWFormat(rf);
%     hold on;
%     %     figure;
%     plot(mean(rf0,1))
% end

%%

tstim = 2/119.5172;
STA_length = 30;
movie_size = size(fitmovie{1});
STA = zeros(movie_size(1),movie_size(2),STA_length,length(fitspikes));
fitframes = movie_size(3);

% for cellind = 1:length(fitspikes)
%     cellind
%     tic
% for i = 1:length(fitspikes{cellind})
%     sp_frame = floor(fitspikes{cellind}(i)/tstim);
%     if sp_frame > STA_length && sp_frame<fitframes
%         STA(:,:,:,cellind) = STA(:,:,:,cellind)+double(fitmovie{1}(:,:,(sp_frame-STA_length+1):sp_frame));
%     end
% end
% toc
% end

for cellind = 1:length(fitspikes)
%     cellind
%     tic
% for i = 1:length(fitspikes{cellind})
    sp_frame = floor(fitspikes{cellind}(:)/tstim);
%     if sp_frame > STA_length && sp_frame<fitframes

sp_rel = find((sp_frame>STA_length)&(sp_frame<fitframes));
    for i = 1:STA_length
        STA(:,:,i,cellind) = sum((fitmovie{1}(:,:,(sp_frame(sp_rel)-STA_length+1)+i)),3);
        
        
    end
%     end
% end
% toc
end

t1 = squeeze(sum(sum(sum(STA(21:40,6:15,:,1:12),1),2),4));
figure; plot(t1);
% save([cell_type '_sta.mat'],'STA');
end

numberCells = [103 112 117 182 12];
% save('t02.mat','t0','numberCells');

% figure; plot((t0{2}-mean(t0{2})));./sum(t0{2}-mean(t0{2})));
% figure; plot(horzcat(t0{:}));

%%

tstim = 2/119.5172;
STC_length = 1;
movie_size = size(fitmovie{1});
STC= zeros(200*STC_length,length(fitspikes{cellind}),1);%length(fitspikes));
fitframes = movie_size(3);

for cellind = 2%:length(fitspikes)
    cellind
    tic
for i = 1:length(fitspikes{cellind})
    sp_frame = floor(fitspikes{cellind}(i)/tstim);
    if sp_frame > STA_length && sp_frame<fitframes
        frtemp = (fitmovie{1}(11:30,6:15,sp_frame-3));
        STC(:,i,1) = frtemp(:);
    end
end
toc
end

STAC = mean(STC(:,1,2),2); figure; imagesc(reshape(STAC,[20 10]));
figure; imagesc((STC-mean(STC(:)))*(STC-mean(STC(:)))' ); colormap parula
[evec,eval] = eig((STC-mean(STC(:)))*(STC-mean(STC(:)))');
figure; plot(diag(eval));
figure; for i = 1:25; subplot(5,5,i); imagesc(reshape(evec(:,end-(i-1)),[20 10]));  end;

% 
% for cellind = 1:length(fitspikes)
%     sp_frame = floor(fitspikes{cellind}(:)/tstim);
%     
%     sp_rel = find((sp_frame>STA_length)&(sp_frame<fitframes));
%     for i = 1:STA_length
%         
%         STC(:,:,i,cellind) = sum((fitmovie{1}(11:30,6:15,(sp_frame(sp_rel)-STA_length+1)+i)),3);
%     end
% end
% 
% t0{type_number} = squeeze(sum(sum(sum(STA,1),2),4));
