function demo_EAPT_video_ShoppingMall(dataname)

    switch nargin
    case 1
        pic_name = strcat('./data/', dataname);
    case 0
        dataname = 'ShoppingMall.mat';
        pic_name = strcat('./data/', dataname);
    end
    addpath('tensor_toolbox_2.6')
    load(pic_name, 'x');
    X = double(x) / 255;
    
    [~, n2, n3] = size(X);
    Xn = permute(X, [3, 1, 2]);

    UU = dctmtx(n2);

    maxP = max(abs(X(:)));

    para.tol       = 1e-4;
    para.gamma     = 0.7;
    para.max_iter  = 100;
    para.X         = Xn;
    para.mu = 5;
    tic
    r = 20;
    r_list = r * ones(n2, 1);
    for i=1:n2
        if i > 1
            r_list(i) = 5;
        end
        if i > 3
            r_list(i) = 1;
        end
        if i > 100
            r_list(i) = 0;
        end
        if i > 120
            r_list(i) = 0;
        end
    end
    [Xhat, E] = EAPT_mr_video(UU, Xn, r_list, r, para);
    cost_time = toc;
    fprintf('cost time is %.2f\n', cost_time);
    Xhat = permute(Xhat, [2, 3, 1]);
    E = permute(E, [2, 3, 1]);

    Xhat = max(Xhat,0);
    Xhat = min(Xhat,maxP);

    FrameRate = 5;
    v = VideoWriter(['./res_' dataname]);
    v.FrameRate = FrameRate;
    open(v);
    for i = 1 : n3
        writeVideo(v, Xhat(:,:,i));
    end
    close(v);

    E(E<0) = 0;
    E(E>1) = 1;
    v = VideoWriter(['./foreground_' dataname]);
    v.FrameRate = FrameRate;
    open(v);
    for i = 1 : n3
        writeVideo(v, E(:,:,i));
    end
    close(v);

    v = VideoWriter(['./orig_' dataname]);
    v.FrameRate = FrameRate;
    open(v);
    for i = 1 : n3
        writeVideo(v, X(:,:,i));
    end
    close(v);
end