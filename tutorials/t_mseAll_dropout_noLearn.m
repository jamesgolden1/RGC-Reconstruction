
clear
resultsName = 'resultsOnlyOnC';
% resultsName = 'resultsNoLearn';
ccB =[];
pwctr = 0;
for pwid =[ 35 70]%70];
    pwctr = pwctr+1;
drctr = 0;
for dr = 100*[.1 .5 .7 .9]
    drctr = drctr+1;
    
ccB =[];
for ii = [1:3 ];
    

 load(['/Volumes/Lab/Users/james/current/RGC-Reconstruction/dat/prosthesis_' num2str(pwid) '_testing_aug13_dr' num2str(dr) '/' resultsName '/stats_' num2str(ii) '.mat']); 
 ccB = [ccB mseAll]; 
%   ccB = [ccB ccAll]; 

end;
% ccB = [ccB mseAll];
% end


msebig(pwctr, drctr) = sqrt(mean(ccB(ccB>0)))/255;

% msebig(drctr) = (mean(ccB(ccB>0)))/1;
end

msebig
% msebig =[]
end
%%
figure; hist(sqrt(ccB(ccB>0))/255,80)
%%
% only on 100*[.1 .5 .7 .9]
%      0.1757    0.1993    0.2086    0.2250
%     0.1978    0.2083    0.2182    0.2310
    
% msebig35 = fliplr([0.1266    0.1274    0.1471    0.1549    0.1518    0.1498    0.1953    0.2311    0.2289]);
% 
% msebig70 = fliplr([0.1159    0.1169    0.1356    0.1475    0.1456    0.1421    0.2179    0.2279    0.2267]);

% msebig35 = fliplr([.12  0.1260    0.1269    0.1497    0.1584    0.1546    0.1520    0.1975    0.2144    0.2259]);
% 
% msebig70 = fliplr([.11 0.1138    0.1156    0.1369    0.1495    0.1471    0.1430    0.1883    0.2091    0.2238]);

% msebig70 = fliplr([.12  0.1260    0.1269    0.1347    0.1584    0.1546    0.1520    0.1975    0.2144    0.2259]);
% 
% msebig35 = fliplr([.11 0.1138    0.1156    0.1249    0.1495    0.1471    0.1430    0.1883    0.2091    0.2238]);
% 
% 
% figure; plot((1:10)/10,(msebig70)/1,'-x','linewidth',4,'markersize',12)
% hold on; plot((1:10)/10,(msebig35)/1,'-xr','linewidth',4,'markersize',12)


figure; plot((0:9)/10,fliplr(msebig(1,:))/1,'-x','linewidth',4,'markersize',12)
hold on; plot((0:9)/10,fliplr(msebig(2,:))/1,'-xr','linewidth',4,'markersize',12)
% figure; plot((1:9)/10,sqrt(mse70)/255,'-x')
% hold on; plot((1:9)/10,sqrt(mse35)/255,'-xr')
% figure; plot(([.5 .6 .8 1:9])/10,(mse70)/1,'-x')
% hold on; plot(([.5 .6 .8  1:9])/10,(mse35)/1,'-xr')
grid on
xlabel('Fraction Remaining RGCs','fontsize',22); ylabel('RMSE','fontsize',22)
set(gca,'fontsize',32);

legend(sprintf('70 um'),sprintf('35 um'));
axis([0.1 1 0.1 0.18])


