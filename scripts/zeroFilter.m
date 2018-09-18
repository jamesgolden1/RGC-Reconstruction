function filterZero = zeroFilter(filterMat,lambda)
% Takes a set of decoding filters in the columns and zeros out weights far
% away from max by weighting with an exponential.
% 
% The smaller lambda is, the closer to the original filter the output is.
% Typical values are 0.001 - 0.05.
%%
[mgr,mgc] = meshgrid(1:100,1:100);

[cmax,cind] = max(abs(filterMat),[],2);
[fmaxc,fmaxr] = ind2sub([100 100],cind);

mgrmat = mgr(:)*ones(1,size(fmaxr,1));
fmaxrmat = ones(size(mgrmat,1),1)*fmaxr';
mgrd = ((mgrmat - fmaxrmat)').^2;

mgcmat = mgc(:)*ones(1,size(fmaxc,1));
fmaxcmat = ones(size(mgcmat,1),1)*fmaxc';
mgcd = ((mgcmat - fmaxcmat)').^2;


% imagesc(signValWN(mind)*staim.*(.1+.9*reshape(exp(-.05*dp(1+paraIndPlus,:).^2),[100 100])));

dp = sqrt(mgrd+mgcd);
clear fmaxrmat fmaxcmat mgrmat mgrd mgcmat mgcd 
% filterMat2 = filterMat;
% filterMat2(dp>5) = 0;
expFilter = 1.2*(exp(-lambda*dp.^2));
expFilter(expFilter>1) = 1;
filterZero = filterMat.*expFilter;
filterZero(1,:) = filterMat(1,:);