function [mtx,mtxi] = interpMtx(M,N, start, stop)
x = linspace(start,stop,M);
x = x.^2;
val = eye(M,M);
xq = linspace(start^2,stop^2,M*N);
mtx = interp1(x,val,xq);

mtx = mtx ./ repmat(sum(mtx),size(mtx,1),1);
mtx(isnan(mtx)) = 0;
mtxi = mtx';
mtxi = mtxi ./ repmat(sum(mtxi),size(mtxi,1),1);
mtxi(isnan(mtxi)) = 0;

mtx = full(mtx);
mtxi = full(mtxi);