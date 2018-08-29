function ladmm(scene_args) 
    % display plots
    plots = 0;

    % retrieve arguments
    scene_name = scene_args.scene_name;
    max_iters = scene_args.max_iters;
    slope = scene_args.psf_slope;
    dx = scene_args.dx;
    dy = scene_args.dy;
    snr = scene_args.snr;
    width = scene_args.width;
    d = scene_args.dark_count;
    lambda = scene_args.lambda_tv;
    lambda_l1 = scene_args.lambda_l1;
    rect_data = scene_args.rect_data;
    N = size(rect_data, 1);
    KK = size(rect_data, 3) / size(rect_data, 1);
    S = 2;      % Upsample ratio for transformation operator
    z_start = 0; 
    z_stop = 1./slope;
    rho = 1e-5;
    z_offset = 200;

    % load data
    disp('Closed-form ADMM Processing');
    b = permute(rect_data,[3 2 1]);

    % compression/stretching matrices
    [T, Z] = interpMtx(KK*N, S, z_start, z_stop);

    % Spatial volume/area for objects (in meters)
    x = linspace(-0.5,0.5,N);
    y = linspace(-0.5,0.5,N);
    z = linspace(0,z_stop,KK*N);
    [grid_z,~,~] = ndgrid(z,y,x);

    % utility functions
    trim_array = @(x) x(S*KK*N/2+1:end-S*KK*N/2, N/2+1:end-N/2, N/2+1:end-N/2);
    pad_array = @(x) padarray(x, [S*KK*N/2, N/2, N/2]);
    square2cube = @(x) reshape(x, [],N,N);
    cube2square = @(x) x(:,:);

    % calculate nlosPsf
    psf = nlosPsf(N,S.*KK.*N,slope,dy,dx);
    psf = circshift(psf,[0 N N]);
    psfFT = fftn(psf);
    data = grid_z.^2.*b;

    % perform LCT reconstruction to get initialization volume
    invpsf = conj(psfFT) ./ (abs(psfFT).^2 + 1./snr);
    tdata  = reshape(T*data(:,:),[S.*KK.*N N N]);
    tvol = trim_array(ifftn(fftn(pad_array(tdata)).*invpsf));
    vol  = real(reshape(Z*tvol(:,:),[KK.*N N N]));
    x  = vol;
    vol(1:z_offset,:,:) = 0;

    % radiometric fall off
    r = 1./grid_z.^2;
    r(1:z_offset,:) = 0;

    % initialize ADMM variables
    z1 = (zeros(KK*N, N, N));
    z2 = (zeros(KK*N, N, N));
    z3 = (zeros(KK*N, N, N));
    z4 = (zeros(KK*N, N, N));
    z5 = (zeros(KK*N, N, N));
    z6 = (zeros(KK*N, N, N));

    u1 = (zeros(KK*N, N, N));
    u2 = (zeros(KK*N, N, N));
    u3 = (zeros(KK*N, N, N));
    u4 = (zeros(KK*N, N, N));
    u5 = (zeros(KK*N, N, N));
    u6 = (zeros(KK*N, N, N));

    % gradient kernels
    d2 = [0 0 0; 0 1 -1; 0 0 0];
    d2 = padarray(d2, [0,0,1]);
    d1 = [0 0 0; 0 1 0; 0 -1 0];
    d1 = padarray(d1, [0,0,1]);
    d3 = zeros(3,3,3);
    d3(2,2,2) = 1; d3(2,2,3) = -1;

    % operator functions
    p2o = @(x) psf2otf(x, [KK*N,N,N]);
    d1FT = p2o(d1);
    d2FT = p2o(d2);
    d3FT = p2o(d3);

    K3 = @(x) real(ifftn(d1FT .*  fftn(x)));
    K4 = @(x) real(ifftn(d2FT .* fftn(x)));
    K5 = @(x) real(ifftn(d3FT .*  fftn(x)));
    K6 = @(x) x;

    K3T = @(x) real(ifftn(conj(d1FT) .*  fftn(x)));
    K4T = @(x) real(ifftn(conj(d2FT) .* fftn(x)));
    K5T = @(x) real(ifftn(conj(d3FT) .*  fftn(x)));
    K6T = @(x) x;

    vec = @(x) x(:);
    Ffun = @(x)  fftn(x);
    Ftfun = @(x) ifftn(x);

    % the full operator and transpose
    K = @(x) [square2cube(Z*cube2square(trim_array(Ftfun(psfFT .* Ffun(pad_array(square2cube(T*x(:,:)))))))); ...
        x; ...
        K3(x); ...
        K4(x);...
        K5(x);...
        K6(x)];
    KT = @(x) square2cube(T' * cube2square(trim_array(Ftfun(conj(psfFT) .* Ffun(pad_array(square2cube(Z'*cube2square(x(1:end/6,:,:)))))))))...
        + x(end/6+1:2*end/6,:,:) ...
        + K3T(x(2*end/6+1:3*end/6,:,:)) ...
        + K4T(x(3*end/6+1:4*end/6,:,:))...
        + K5T(x(4*end/6+1:5*end/6,:,:))...
        + K6T(x(5*end/6+1:6*end/6,:,:));

    % quadratic operator and transpose
    K1 = @(x) square2cube(Z*cube2square(trim_array(Ftfun(psfFT .* Ffun(pad_array(square2cube(T*x(:,:))))))));
    K1T = @(x) square2cube(T' * cube2square(trim_array(Ftfun(conj(psfFT) .* Ffun(pad_array(square2cube(Z'*x(:,:))))))));
    shrinkage = @(a,kappa) repmat(max(0, 1-kappa./sqrt(a(:,:,:,1).^2 + a(:,:,:,2).^2 + a(:,:,:,3).^2)),1,1,1,3).*a;

    % compute the eigenvalues of the linear operator using power iteration
    if exist('data/ladmm_combined_knorm.mat','file')
        load('data/ladmm_combined_knorm.mat');
    else
        [Knorm, Ksmall]= compute_operator_norm(K, KT, [KK*N N N]);
        save('data/ladmm_combined_knorm.mat', 'Knorm', 'Ksmall');
    end
    fprintf('Operator eigenvalue ratio: %.02f\n', Knorm/Ksmall);

    % run linearized ADMM
    mu = rho * (Knorm).^2;
    for ii = 1:max_iters   
       tic

       % keep track of previous iteration
       xold = x;
       z1old = z1;
       z2old = z2;
       z3old = z3;
       z4old = z4;
       z5old = z5;
       z6old = z6;

       % z1 update -- quadratic term
       v = K1(x) + u1;
       k = (r.^2 + rho*d - rho.*v.*r) ./ (2.*rho.*r);
       z1 = -k + sqrt(k.^2 + (r.*b + rho.*d.*v - r.*d)./(rho.*r)); 

       % z2 update -- non-negativity
       v = x + u2;
       z2 = max(0,v);

       % z3/z4/z5 update -- TV regularizer
       v = cat(4, K3(x), K4(x), K5(x)) + cat(4, u3, u4, u5);
       kappa = lambda/rho;
       shrunk = shrinkage(v, kappa);
       z3 = squeeze(shrunk(:,:,:,1));
       z4 = squeeze(shrunk(:,:,:,2));
       z5 = squeeze(shrunk(:,:,:,3));

       % z6 update -- Sparsity penalty
       v = x + u6;
       kappa = lambda_l1/rho;
       z6 = max(v - kappa,0) - max(-v - kappa,0);   

       % don't reconstruct up to z_offset
       z1(1:z_offset,:,:) = 0;
       z2(1:z_offset,:,:) = 0;
       z3(1:z_offset,:,:) = 0;
       z4(1:z_offset,:,:) = 0;
       z5(1:z_offset,:,:) = 0;
       z6(1:z_offset,:,:) = 0;

       % u update
       u1 = u1 + K1(x) - z1;
       u2 = u2 + x - z2;
       u3 = u3 + K3(x) - z3;
       u4 = u4 + K4(x) - z4;
       u5 = u5 + K5(x) - z5;
       u6 = u6 + K6(x) - z6;

       % x update
       v1 = z1 - u1;
       v2 = z2 - u2;
       v3 = z3 - u3;
       v4 = z4 - u4;
       v5 = z5 - u5;
       v6 = z6 - u6;
       x = x - (rho/mu) * (KT(K(x)) - KT([v1; v2; v3; v4; v5; v6]));

       times(ii) = toc;   
       fprintf('Iteration: %d, Elapsed: %.01f\n', ii, times(ii));

       if plots 
           % compute loss/primal dual residuals
           K1x = K1(x);
           K1x = r.*K1x + d;
           mask = K1x>0;
           loss(ii) = -sum(log(K1x(mask)).*b(mask),'omitnan') + sum(K1x(:)) + sum(logfactorial(b(:)),'omitnan');
           rpri(ii) = norm(vec(K(xold) - [z1;z2;z3;z4;z5;z6]));
           rdual(ii) = norm(vec(rho.*(KT([z1;z2;z3;z4;z5;z6] - [z1old;z2old;z3old;z4old;z5old;z6old]))));

           set(0,'defaultaxesfontsize',14); 
           subplot(3,3,1);
           plot(real(loss)); title('Objective'); grid on;
           subplot(3,3,2);
           plot(real(rpri)); title('Primal residual'); grid on;
           subplot(3,3,3);
           plot(real(rdual)); title('Dual residual'); grid on;
           subplot(3,3,4);
           imagesc(squeeze(max(abs(vol),[],1))); colorbar; set(gca,'xtick',[],'ytick',[]);
           xlabel('x'); ylabel('y');
           subplot(3,3,5);
           imagesc(squeeze(max(abs(vol),[],2))); colorbar; set(gca,'xtick',[],'ytick',[]);
           xlabel('x'); ylabel('z'); title('LCT');
           subplot(3,3,6);
           imagesc(squeeze(max(abs(vol),[],3))); colorbar; set(gca,'xtick',[],'ytick',[]);
           xlabel('y'); ylabel('z');
           subplot(3,3,7);
           imagesc(squeeze(max(abs(x),[],1))); colorbar; set(gca,'xtick',[],'ytick',[]);
           xlabel('x'); ylabel('y');
           subplot(3,3,8);
           imagesc(squeeze(max(abs(x),[],2))); colorbar; set(gca,'xtick',[],'ytick',[]);
           xlabel('x'); ylabel('z'); title('Reconstruction');
           subplot(3,3,9);
           imagesc(squeeze(max(abs(x),[],3))); colorbar; set(gca,'xtick',[],'ytick',[]);
           xlabel('y'); ylabel('z');          
           set(gcf,'color','white');
           drawnow; 
       end
       if mod(ii,50) == 0
            save(sprintf('results/%s_%d_%0.02e_%0.02e.mat',scene_name,ii,lambda,lambda_l1),'x','lambda','rho','ii','d','times','lambda_l1','snr','vol');   
       end
    end
end