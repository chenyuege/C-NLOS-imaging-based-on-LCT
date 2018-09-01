function [mtx,mtxi] = resamplingOperator(M)
% Local function that defines resampling operators

mtx = sparse([],[],[],M.^2,M,M.^2);

x = 1:M.^2;
mtx(sub2ind(size(mtx),x,ceil(sqrt(x)))) = 1;
mtx  = spdiags(1./sqrt(x)',0,M.^2,M.^2)*mtx;
mtxi = mtx';

K = log(M)./log(2);
for k = 1:round(K)
    mtx  = 0.5.*(mtx(1:2:end,:)  + mtx(2:2:end,:));
    mtxi = 0.5.*(mtxi(:,1:2:end) + mtxi(:,2:2:end));
end
end