function psf = nlosPsf(U,V,slope,dx,dy)

x = linspace(-1,1,2*U);
y = linspace(-1,1,2*U);
z = linspace(0,2,2*V);
[grid_z,grid_y,grid_x] = ndgrid(z,y,x);

% Define PSF
psf = abs((4*slope.^2).*((dx.*grid_x).^2 + (dy.*grid_y).^2) - grid_z);
psf = double(psf == repmat(min(psf,[],1),[2*V 1 1]));
psf = psf ./ sum(psf(:));

end