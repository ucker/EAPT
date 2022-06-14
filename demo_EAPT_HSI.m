function [Error, PSNR, SSIM] = demo_EAPT_HSI(dataname)

    switch nargin
    case 1
        pic_name = strcat('./data/', dataname);
    case 0
        dataname = 'toys.mat';
        pic_name = strcat('./data/', dataname);
    end
    addpath('tensor_toolbox_2.6')
    addpath('toolbox')
    load(pic_name, 'orig');
    X = double(orig) / 255;
    r = 5;
    
    [n1, n2, n3] = size(X);

    Xn = permute(X, [3, 1, 2]);

    rhos = 0.3;
    ind = find(rand(n1*n2*n3,1)<rhos);
    Xn(ind) = rand(length(ind), 1);


    maxP = max(abs(X(:)));

    para.tol       = 1e-4;
    para.gamma     = 0.85;
    para.X         = X;
    para.maxP = maxP;
    para.mu = 500;
    para.dataname = dataname;
    UU = dctmtx(n2);
    tic
    r_list = r * ones(n2, 1);
    
    for i=1:n2
        if i > 100
            r_list(i) = 3;
        end
        if i > 300
            r_list(i) = 1;
        end
        if i > 500
            r_list(i) = 0;
        end
    end
    if strcmp(dataname, 'sponges_ms.mat')
        for i = 1:n2
            if i > 50
                r_list(i)=3;
            end
            if i > 200
                r_list(i)=2;
            end
            if i > 300
                r_list(i)=1;
            end
            if i > 400
                r_list(i)=0;
            end
        end
    end
    [Xhat, ~] = EAPT_mr(UU, Xn, r_list, r, para);
    cost_time = toc;
    fprintf('cost time is %.2f\n', cost_time);
    Xhat = permute(Xhat, [2, 3 ,1]);

    Xhat = max(Xhat,0);
    Xhat = min(Xhat,maxP);
    PSNR = psnr(X, Xhat, maxP);
    SSIM = ssim(X, Xhat);
    Error = norm(X(:) - Xhat(:), 2) / norm(X(:), 2);
    fprintf('Error=%f, PSNR, SSIM=%2.2f, %.4f\n', Error, PSNR, SSIM);

    end