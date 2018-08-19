function [scene_args] = getSceneArgs(scene)    
    fprintf('Scene: %d\n', scene);

    switch scene
        case {1}
            scene_name = 'data_mannequin'
            load(['data/' scene_name '.mat'], 'rect_data', 'width');
            dx = 1;
            dy = 1;
            snr = 1e4;
            dark_count = 5;
            lambda_tv = 1e-4;
            lambda_l1 = 1e-3;
            max_iters = 400;
        case {2}
            scene_name = 'data_exit_sign'
            load(['data/' scene_name '.mat'], 'rect_data', 'width');
            dx = 1;
            dy = 1;
            snr = 1e4;
            dark_count = 5;
            lambda_tv = 3e-5;
            lambda_l1 = 1e-4;
            max_iters = 400;
        case {3}
            scene_name = 'data_s_u'
            load(['data/' scene_name '.mat'], 'rect_data', 'width');
            dx = 1;
            dy = 1;
            snr = 1e5;
            dark_count = 3;
            lambda_tv = 5e-4;
            lambda_l1 = 1e-8;
            max_iters = 400;
        case {4}
            scene_name = 'data_outdoor_s'
            load(['data/' scene_name '.mat'], 'rect_data', 'width');
            snr = 1e3;       
            dx = 1; 
            dy = 1;  
            dark_count = 4;
            lambda_tv = 4.2e-4;
            lambda_l1 = 1.3e-2;
            max_iters = 800;
        case {5}
            scene_name = 'data_diffuse_s'            
            load(['data/' scene_name '.mat'], 'rect_data', 'width');
            dx = 1; 
            dy = 1;
            snr = 1e4;
            dark_count = 3;
            lambda_tv = 3.25e-4;
            lambda_l1 = 2.75e-3;
            max_iters = 800;
    end
    
    % Downsample data to 16 picoseconds
    bin_resolution = 4e-12;
    for k = 1:2        
        rect_data = rect_data(:,:,1:2:end) + rect_data(:,:,2:2:end);
        bin_resolution = 2 * bin_resolution;
    end
    N = size(rect_data,1);        % Spatial resolution of data
    M = size(rect_data,3);        % Temporal resolution of data
    c = 3e8;
    range = M.*c/2.*bin_resolution; % Maximum range for histogram
    psf_slope = width / range;
    
    scene_args.scene_name = scene_name;
    scene_args.psf_slope = psf_slope;
    scene_args.dx = dx;
    scene_args.dy = dy;
    scene_args.snr = snr;
    scene_args.width = 2 * width;
    scene_args.dark_count = dark_count;
    scene_args.lambda_tv = lambda_tv;
    scene_args.lambda_l1 = lambda_l1;
    scene_args.rect_data = rect_data;
    scene_args.max_iters = max_iters;
    
end