%%
for i=1:1
    
    nband = size(img{1},3);
    
    if(nband == 3)
        img_bf{i} = bfilter2(img{i}, 5, [3 0.05]);
    else
        for j=1:nband
            img_bf{i}(:,:,j) = bfilter2(img{i}(:,:,j), 5, [3 0.05]);
        end
    end
    
    model_pca_bf{i} = imgpca_2(uint8(mat2gray(img_bf{i})*255),'auto', nband);
    
    fset_bf{i}(:,:,1) = mat2gray(model_pca_bf{i}.scores(:,:,1))*2 - 1;
    fset_bf{i}(:,:,2) = mat2gray(model_pca_bf{i}.scores(:,:,2))*2 - 1;
    
end
%%
for i=1:1
    fset{i}(:,:,1) = mat2gray(model_pca{i}.scores(:,:,1))*2 - 1;
    fset{i}(:,:,2) = mat2gray(model_pca{i}.scores(:,:,2))*2 - 1;
    
end
%%
for i=1:40
    img_ace{i} = ace( img_gray{i}, 5, 25, 10, 15, 1, img_name{i});
end

%%
for i=1:40
    fset{i} = img{i}*2 - 1;
end



%%
for i=1:40
    
    for j=1:1
        nu = [0.05 0.1 0.2];
        range = [-1 1];

        for k=1:size(nu,2)
            if(sum(img_ace{i}(:)) ~= 0)
                [R] = imclassify2(fset{i}, fset_bf{i}, img_ace{i}, zeros(size(img_ace{i})), 'oc-svm', nu(k), range, img_name{i}, sprintf('Results_2013_02_03_2_ocsvm/%d',j));
                close all;
            end
        end
    end
    
end

%%
for i=1:1
    
    for j=1:3
        nu = [0.01 0.05 0.1];
        range = [-1 1];

        for k=1:size(nu,2)
            if(sum(img_ace{i}(:)) ~= 0)
                [R] = imclassify(fset{i}, img_ace{i}, zeros(size(img_ace{i})), 'oc-svm', nu(k), range, img_name{i}, sprintf('Results_2013_02_14_gozdehoca_ace/%d',j));
                close all;
            end
        end
    end
    
end


%% MRF Relaxation / Simulated Annealing with Gibbs Sampler %%

    nLevel = 7; % # of temp levels
    step = 50;  % # of step in each level
    c = 0.5;    % cooling schedule
    Ti = 2;     % initial temperature

    mrf = [];
    mrf.sorm = 1;
    
    mrf.params.Ke = 0.3;   
    mrf.params.Kl = 0.01; % Since lengths of segments are similar, Kl parameter is kept so small.  
    mrf.params.Kc = 0.3; 
    mrf.params.Ki = 0.4; 
    
    mrf.params.Vo = 0.7;
    
    mrf.params.eratio = (1/25);
    mrf.params.allowjunc = 0;

%   <<Eski parametre degerleri>>
%
%     mrf.params.Ke = 0.17; 
%     mrf.params.Kl = 0.13; 
%     mrf.params.Kc = 0.3; 
%     mrf.params.Ki = 0.2; 
% 
%     mrf.params.eratio = (1/20);
%
    
  %%  
%%for i=1:15
    
    % detect all possible connections
    [linelist{i}, Cliques{i}] = SORMNetworkConstruct(CC{i}, [], clres{i}, 50, 0.1);
%%
for i=1:3    
    
    % Draw all possible connections
    % DrawLabelledLineList(clres{i}, linelist{i}, ones(size(linelist{i},2))); 

    % save connection info for later use
    save_linelist_info(linelist{i}, sprintf('linelist_info_%s.txt', img_name{i}));   
    
    % initial labeling / not important for stochastic relaxation
    labels = ones(1,size(linelist{i},2))';

    % optimization
    [labels_GlbMin{i}] = SAwithGibbsSampler(linelist{i}, labels, step, nLevel, Ti, mrf, Cliques{i}, 0);

   
    % draw optimum network
    DrawLabelledLineList(clres{i}, linelist{i}, labels_GlbMin{i});  
    
end

%%
for i=1:15
    if(~isempty(gt{i}))
        % draw optimum network
        DrawLabelledLineList(gt{i}, linelist{i}, opt_labels{i});  
    end
end
%%
    nLevel = 7; % # of temp levels
    step = 6;  % # of step in each level
    c = 0.5;    % cooling schedule
    Ti = 2;     % initial temperature

    mrf = [];
    mrf.sorm = 1;
    
    mrf.params.Ke = 0.3;   
    mrf.params.Kl = 0.01; % Since lengths of segments are similar, Kl parameter is kept so small.  
    mrf.params.Kc = 0.3; 
    mrf.params.Ki = 0.4; 
    
    mrf.params.Vo = 0.7;
    
    mrf.params.eratio = (1/25);
    mrf.params.allowjunc = 0;

for i=1:15
    if(~isempty(gt{i}))
        %[linelist{i}, Cliques{i}] = SORMNetworkConstruct(CC{i}, [], clres{i}, 50, 0.1);

        iLabels = ([linelist{i}.prob] > 0.5)';

        icm_labels{i} = MRFRelaxWithICM(linelist{i}, iLabels, step, mrf, Cliques{i});
    end
end

%%

for i=1:15
    
    if(~isempty(gt{i}))
        
        ext{i} = ones(size(img{i},1),size(img{i},2));
   
        for j=1:size(linelist{i},2)
            if(opt_labels{i}(j) == 1)
                [~,~,ext{i}] = bresenham(ext{i},[ linelist{i}(j).s(2), ...
                                                      linelist{i}(j).s(1); ...
                                                      linelist{i}(j).e(2),...
                                                      linelist{i}(j).e(1)],0);          
            end
        end
        
        ext{i} = imcomplement(ext{i});   
        figure;imshow(ext{i},[]);
    end
end

%%
for i=1:15 
    if(~isempty(gt{i}))
        ref{i} = bwmorph(gt{i},'skel',Inf);
        ref_buf{i} = bwmorph(ref{i},'dilate',5);
        ext_buf{i} = bwmorph(ext{i},'dilate',5);
    end
end

%%
for i=1:15
    if(~isempty(gt{i}))
        match_ref{i} = ext{i} & ref_buf{i};

        match_ext{i} = ref{i} & ext_buf{i};
        
        completeness{i} = sum(match_ref{i}) / sum(ref{i});
        correctness{i}  = sum(match_ext{i}) / sum(ext{i});
    end
end

%%

for i=1:15
    if(~isempty(gt{i}))
        eval{i} = zeros(size(img{i},1),size(img{i},1),3);
        R = bwmorph(ref{i} & imcomplement(ext_buf{i}), 'dilate');
        G = bwmorph(match_ext{i},'dilate');
        B = bwmorph(ext{i} & imcomplement(ref_buf{i}),'dilate');
        
        W = imcomplement(R) & imcomplement(G) & imcomplement(B);
        
        eval{i}(:,:,1) = R + 0.5 * W;
        eval{i}(:,:,2) = G + 0.5 * W;
        eval{i}(:,:,3) = B + 0.5 * W;
        
        figure;imshow(eval{i},[])
        
        imwrite(eval{i}, sprintf('%s_eval_result.png',img_name{i}));
    end
end
%%


