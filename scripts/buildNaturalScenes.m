% 
% 
% Save natural scene images

images1 = dir([reconstructionRootPath '\dat\FoliageBig\*.tif']);

im0 = imread([reconstructionRootPath '\dat\FoliageBig\' images1(1).name]);
figure; imagesc(im0)

% mov = uint8(zeros(size(im0,1),size(im0,2),1));%length(images1));
% movsm = uint8(zeros(100,100,475*length(images1)));
movsm = uint8(zeros(100,100,12000));

vsm = zeros(12000,1);

pctr = 0;
bctr = 1;
for i =1:length(images1)
    
%     while pctr <= 12350
    im1 = imread([reconstructionRootPath '\dat\FoliageBig\' images1(i).name]);
%     imrs = reshape(rgb2gray(im1),size(im0,1),size(im0,2));
    imrs = rgb2gray(im1);
%     mov(:,:,i) = imrs;
    % movsm(:,:,35*(i-1)+1:35*i) = permute(reshape(imrs(1:500,1:700)', 35,100,100),[3 2 1]);
    
    for k = 1:floor(size(im1,1)/100)
        for j = 1:floor(size(im1,2)/100)
            pctr = pctr+1;
            % movsm(:,:,35*(i-1)+5*(k-1) + j) = imrs(100*(k-1)+1:100*k, 100*(j-1)+1:100*j);
            imrstmp = imrs(100*(k-1)+1:100*k, 100*(j-1)+1:100*j);
            movsm(:,:,pctr) = imrstmp;
            
            
            vsm(pctr) = var(double(imrstmp(:)));
        end
    end
    
%     end
    if pctr >= 12350;
        save([ reconstructionRootPath  '\dat\movsm_' num2str(bctr) '.mat'],'movsm');
        pctr = 0;
        bctr = bctr+1;
    end
    
    
end
pctr
% figure; imagesc(imrs); colormap gray
% figure; imagesc(movsm(:,:,i)); colormap gray;

%%
% figure;
% spctr=0;
% for k = 0+[1:5]%19
%     for j =0+[1:5]%25
%         spctr=spctr+1;
% %         movsm(:,:,35*(i-1)+5*(k-1) + j) = imrs(100*(k-1)+1:100*k, 100*(j-1)+1:100*j);
% %         subplot(19,25,spctr);
%         subplot(5,5,spctr);
%         imagesc(movsm(:,:,19*(k-1)+1 + j)); colormap gray
%     end
% end