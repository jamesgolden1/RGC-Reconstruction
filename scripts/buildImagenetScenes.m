% buildImagenetScenes
% 
% Save natural scene images from the imagenet directory on Bertha in
% /Volumes/Lab/Users/nikhilp/imagenet

% Get dir names
imagenetDir = dir(['/Volumes/Lab/Users/nikhilp/imagenet']);
dirFlags = [imagenetDir.isdir];
imagenetFolders = imagenetDir(dirFlags);

% Test on one image
% images1 = dir(['/Volumes/Lab/Users/nikhilp/imagenet/' imagenetFolders(10).name '/*.JPEG']);
% im0 = imread(['/Volumes/Lab/Users/nikhilp/imagenet/' imagenetFolders(10).name '/' images1(1).name ]);
% figure; imagesc(im0)

movsm = uint8(zeros(100,100,12000));

vsm = zeros(12000,1);

pctr = 0;
bctr = 22;

% Iterate folders
for folderInd = [200:280]%3:100%length(length(imagenetFolders))
    

images1 = dir(['/Volumes/Lab/Users/nikhilp/imagenet/' imagenetFolders(folderInd).name '/*.JPEG']);

% Iterate images, save by 100x100 blocks
for i =1:length(images1)
    imgFileName = (['/Volumes/Lab/Users/nikhilp/imagenet/' imagenetFolders(folderInd).name '/' images1(i).name ]);
    imgInfo = imfinfo(imgFileName);
    if strcmpi(imgInfo.ColorType, 'truecolor')
        im1 = imread(['/Volumes/Lab/Users/nikhilp/imagenet/' imagenetFolders(folderInd).name '/' images1(i).name ]);
   
    % imrs = reshape(rgb2gray(im1),size(im0,1),size(im0,2));
    imrs = rgb2gray(im1);
    % mov(:,:,i) = imrs;
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
    
    if pctr >= 12350;
        save([ reconstructionRootPath  '/dat/imagenetBlocks/movsm_' num2str(bctr) '.mat'],'movsm','vsm');
        pctr = 0;
        bctr = bctr+1;
    end
    
    
    end
end
[folderInd i]

end

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