

% % % % % % % % % Feb 6
% evArr = [.2 .1 .3 .4];
% 
% % evArr = [.01 .02 .03 .04];
% trainFraction = [.16 .5 .66 .83 1];
% 
% for evInd = 1:4
% for trainFractionInd = 1:5
% [evInd trainFractionInd]
% filterFile = ['ns100_r2_10/filters3_ns100_feb6_sh9_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];

% % % % % % % % % % % % Feb 22
evArr = [.01 .02 .03 .04 .1 .2];
trainFraction = [.16 .5 .66 .83 1];
% already did trainFraction =1
for evInd = 4%1:6
for trainFractionInd = 5%1:5
    [evInd trainFractionInd]
% if ~(evInd==1 && trainFractionInd==1)

% filterFile = ['ns100_r2_10/filters3_ns100_feb6_sh9_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];
% filterFile = ['pixium15_100/filters3_pix1_feb6_sh0_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];
filterFile = ['pixium15_sm/filters_pix_feb21_sh0_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];



%     shiftval = 9;
    shiftval = 0;
    
% mosaicFile = '_mosaicAll_35336498'; %pOpt.mosaicFile = mosaicFile;
% filterFile = ['pixium15_100/filters3_pix1_feb6_sh0_sv' sprintf('%2d',100*evArr(evInd)) '_tr' sprintf('%2d',100*trainFraction(trainFractionInd))  mosaicFile];
%     shiftval = 0;
        

            pOpt.filterFile = filterFile;
            [movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalReconSingle(pOpt);
            
    movieRecon = movrecons_on_offHealthy;

%     clear movieComb
    szLen = 596;
%     movieComb = 255*irradianceFraction*pulseDutyCycle*ieScale(movieRecon(:,:,1:szLen-shiftval+1));
%     movieComb(:,rsFactor*stimSize+[1:rsFactor*stimSize],:) = 255*irradianceFraction*pulseDutyCycle*ieScale(testmovieshort(:,:,shiftval+1:szLen+1));
%     maxc = max(movieComb(:)); minc = min(movieComb(:));


%     mc01 = (movieRecon(:,:,1:szLen-shiftval+1));
%     mc02 = (testmovieshort(:,:,shiftval+1:szLen+1));
%     figure; subplot(121); hist(mc01(:),40); subplot(122); hist(mc02(:),40);
%     
%     fr0=25;
%     mc011 = mc01(:,:,fr0); mc022 = mc02(:,:,fr0);
%     
%     figure; subplot(121); hist(mc011(:),40); subplot(122); hist(mc022(:),40);
%     mc011rs = ieScale(mc011); mc022rs = ieScale(mc022);
%     
%     figure; 
%     subplot(131); hist(mc011rs(:)-mean(mc011rs(:)),40); 
%     subplot(132); hist(mc022rs(:)-mean(mc022rs(:)),40);
%     subplot(133); hist(sqrt((mc011rs(:)-mean(mc011rs(:))-(mc022rs(:)-mean(mc022rs(:)))).^2),40);

    mc1 = ieScale(movieRecon(:,:,1:szLen-shiftval+1));
    mc2 = ieScale(testmovieshort(:,:,shiftval+1:szLen+1));
    
%     errmov =(mc1-mean(mc1(:)))-(mc2-mean(mc2(:)));
%     errtot = ((errmov.^2));

    mc1rs = RGB2XWFormat(mc1);
    mc2rs = RGB2XWFormat(mc2);
    
    mc1rz = mc1rs - ones(size(mc1rs,1),1)*mean(mc1rs,1);    
    mc2rz = mc2rs - ones(size(mc2rs,1),1)*mean(mc2rs,1);
    
    errmov =(mc1rz)-(mc2rz);
    errtot = ((errmov.^2));

    mse(evInd,trainFractionInd) = sqrt(mean(errtot(:)));
    mss(evInd,trainFractionInd) = sqrt(var(errtot(:)));
    
    trainSizeMat(evInd,trainFractionInd) = trainFraction(trainFractionInd);
    evMat(evInd,trainFractionInd) = evArr(evInd);
    
    clear mc1 mc2 errmov errtot movieRecon
end
mse(evInd,:)
end

% mse =

%     0.1502    0.1474    0.1474    0.1470
%     0.1523    0.1526    0.1533    0.1530
%     0.1505    0.1477    0.1481    0.1472
%     0.1509    0.1486    0.1483    0.1476
% 
% mse =
% 
%     0.1369    0.1348    0.1346    0.1344
%     0.1384    0.1383    0.1385    0.1384
%     0.1370    0.1350    0.1349    0.1345
%     0.1375    0.1355    0.1353    0.1347
% 
% trainSizeMat =...
%     [0.1600    0.5000    0.6600    0.8300    1.0000;...
%     0.1600    0.5000    0.6600    0.8300    1.0000;...
%     0.1600    0.5000    0.6600    0.8300    1.0000;...
%     0.1600    0.5000    0.6600    0.8300    1.0000];
%     
% evMat=...
%     [0.2000    0.2000    0.2000    0.2000    0.2000;...
%     0.1000    0.1000    0.1000    0.1000    0.1000;...
%     0.3000    0.3000    0.3000    0.3000    0.3000;...
%     0.4000    0.4000    0.4000    0.4000    0.4000];
%     
% mseh =...
%     [0.1423    0.1370    0.1359    0.1354    0.1352;...
%     0.1409    0.1369    0.1363    0.1360    0.1359;...
%     0.1432    0.1369    0.1367    0.1356    0.1354;...
%     0.1455    0.1383    0.1382    0.1361    0.1361];
% % prosthesis
% msep =...
%     [0.1687    0.1592    0.1564    0.1549    0.1536;...
%     0.1740    0.1619    0.1586    0.1549    0.1534;...
%     0.1805    0.1738    0.1671    0.1637    0.1610;...
%     0.1816    0.1753    0.1729    0.1661    0.1635];
% 
% 
% 
figure; 
plot(1e6*trainSizeMat',mse','-x','linewidth',4)
hold on;
% plot(1e6*trainSizeMat',mseh','-x','linewidth',4)
grid on; axis([0 1e6 .13 .19]);
xlabel('Training Set Size','fontsize',15); 
ylabel('Mean Square Error - Reconstruction','fontsize',15);
set(gca,'fontsize',15)
legend('10% SVs','20% SVs','30% SVs','40% SVs');