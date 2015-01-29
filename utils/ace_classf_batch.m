%% ACE

for i=1:4
    
    img_ace{i} = ace( rgb2gray(img{i}), 6, 15, 15, 30, 1, img_name{i});
    
    figure;imshow(img_ace{i},[]);
end



%% PCA
for i=1:4
    pca_model{i} = imgpca_2(uint8(img{i}*255), 'auto', 3);
end

%% Bilateral filtering
for i=3:3
    img_bf{i} = bfilter2(img{i});
end

%% feature set 
fset = [];

for i=3:3
    fset{i,1} = img{i} * 2 - 1;
    
    fset{i,2}(:,:,1) = mat2gray(pca_model{i}.scores(:,:,1),[0 255]) * 2 -1;
    fset{i,2}(:,:,2) = mat2gray(pca_model{i}.scores(:,:,2),[0 255]) * 2 -1;
    
    fset{i,3} = img_bf{i} * 2 - 1;
end


%% Classification (OCSVM)

ocsvm_qmetrics = cell(4,2);

for i=3:3

    for j=1:1
        nu = [0.01 0.05 0.1];
        range = [-1 1];

        for k=1:size(nu,2)
            [ocsvm_qmetrics{i,j}] = imclassify(fset{i,j}, bwmorph(img_ace{i},'thicken',1), ref{i}, 'oc-svm', nu(k), range, img_name{i}, sprintf('Results_2013_01_08_ace/%d',j));
            close all;
        end
    end
end

clear i j k range nu

%% 

for i=1:3
    [Rmin] = gmm_classification_ex((fset{i,2}*0.5 + 0.5), mask{i}, ref{i}, img_name{i}, 'Results_2012_12_28/gmm');
    close all
end

%% Template Matching Filter

% 5 - ocsvm_clres_img1_(nu_0.01)_(g =   8).png
% 7 - ocsvm_clres_img2_(nu_0.10)_(g =  32).png
%10 - ocsvm_clres_img3_(nu_0.10)_(g = 256).png
% 7 - ocsvm_clres_img4_(nu_0.10)_(g =  32).png

ace_clres_flt    = cell(1,4);
ace_clres_flt_th = cell(1,4);

for i=1:4
    
    ace_clres_flt{i} = RoadTemplateMatchingFilterEx2(ace_clres{i}, 5, 15, [5 15], 0, 2);
    ace_clres_flt_th{i} = cc_threshold(ace_clres_flt{i} > 0.11, 100);
    
    figure;
    subplot(2,2,1); imshow(img_ace{i}, []);
    subplot(2,2,2); imshow(ace_clres{i}, []);
    subplot(2,2,3); imshow(ace_clres_flt{i}, []);
    subplot(2,2,4); imshow(ace_clres_flt_th{i}, []);
end


%% SORM
CC_ace = cell(1,4);

for i=3:3
    [CC_ace{i}, ~] = sorm(ace_clres_flt_th{i}, 20, 13, 5);
end

%% MRF Relaxation / Simulated Annealing with Gibbs Sampler %%

% !!!! mrf_batch_run.m dosyasindan buraya kopyalandi. !!!!!!!!

    nLevel = 6; % # of temp levels
    step = 50;  % # of step in each level
    c = 0.5;    % cooling schedule
    Ti = 1;     % initial temperature

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
    
    
for i=3:3
    
    % detect all possible connections
    [linelist_ace{i}, Cliques_ace{i}] = SORMNetworkConstruct(CC_ace{i}, [], ace_clres{i}, 45, 0.0);
    
    % Draw all possible connections
    DrawLabelledLineList(ace_clres{i}, linelist_ace{i}, ones(size(linelist_ace{i},2))); 

    % save connection info for later use
    save_linelist_info(linelist_ace{i}, sprintf('linelist_ace_info_%s.txt', img_name{i}));   
    
    % initial labeling / not important for stochastic relaxation
    labels = ones(1,size(linelist_ace{i},2))';

    % optimization
    [labels_GlbMin{i}] = SAwithGibbsSampler(linelist_ace{i}, labels, step, nLevel, Ti, mrf, Cliques_ace{i});

    %opt_label = labels_GlbMin{i};
    
    %save(sprintf('2-opt_labels_%s.mat',img_name{i}), 'opt_label');
    
    % draw optimum network
    DrawLabelledLineList(ace_clres{i}, linelist_ace{i}, labels_GlbMin{i});  
    
end

clear labels
%%

for i=2:2
    
    tmp = [cursor_info_ace{i}.Position];
    tmp = [tmp(2:2:end);tmp(1:2:end)]';

    line_s = tmp(1:2:end,:);
    line_e = tmp(2:2:end,:);

    new_lines = [];
    for j=1:size(line_s,1)
        new_lines{i}(j).s = line_s(j,:);
        new_lines{i}(j).e = line_e(j,:);
    end

    [linelist_ace_ex{i}, labels_final_ace_ex{i}] = AddLinesToLineList(new_lines{i}, linelist_ace{i}, opt_labels_ace{i}, CC_ace{i}, img{i});

    DrawLabelledLineList(img{i}, linelist_ace_ex{i}, labels_final_ace_ex{i}); 

    % labels_final_ace_ex{1}([44 152]) = [0 1];
    % labels_final_ace_ex{2}([112 330]) = [0 0];
    % labels_final_ace_ex{3}([147 305 346]) = [0 0 0];
end

clear tmp new_lines line_s line_e 
%% Quality measures before post-processing

for i=2:2
    [completeness{i}, correctness{i}, rmse{i}] = NetworkQuality(lineSet{i}, linelist_ace{i}, opt_labels_ace{i}, img{i});

    quality{i} = completeness{i}*correctness{i} / (completeness{i} - completeness{i}*correctness{i} + correctness{i});
    
    fprintf('Quality metrics(%s): completeness=%6.4f, correctness=%6.4f, quality=%6.4f, rmse=%6.4f\n', ...
        img_name{i}, completeness{i},correctness{i},quality{i},rmse{i});
end

% Quality metrics(img1): completeness=0.6374, correctness=0.9878, quality=0.6325, rmse=1.0455
% Quality metrics(img2): completeness=0.8345, correctness=0.9201, quality=0.7781, rmse=1.1858
% Quality metrics(img3): completeness=0.8629, correctness=0.9870, quality=0.8532, rmse=1.0458
% Quality metrics(img3): completeness=0.8909, correctness=0.9380, quality=0.8414, rmse=1.1028
% Quality metrics(img4): completeness=0.7649, correctness=0.9933, quality=0.7610, rmse=1.1992


%% Quality measures after post-processing

for i=2:2
    [completeness{i}, correctness{i}, rmse{i}] = NetworkQuality(lineSet{i}, linelist_ace_ex{i}, labels_final_ace_ex{i}, img{i});

    quality{i} = completeness{i}*correctness{i} / (completeness{i} - completeness{i}*correctness{i} + correctness{i});
    
    fprintf('Quality metrics(%s): completeness=%6.4f, correctness=%6.4f, quality=%6.4f, rmse=%6.4f\n', ...
        img_name{i}, completeness{i},correctness{i},quality{i},rmse{i});
end


% Quality metrics(img1): completeness=0.6816, correctness=0.9894, quality=0.6767, rmse=1.0596
% Quality metrics(img2): completeness=0.8979, correctness=0.9182, quality=0.8313, rmse=1.1889
% Quality metrics(img3): completeness=0.9588, correctness=0.9805, quality=0.9408, rmse=1.0475
% Quality metrics(img3): completeness=0.9728, correctness=0.9430, quality=0.9188, rmse=1.0895
% Quality metrics(img4): completeness=0.8270, correctness=0.9899, quality=0.8201, rmse=1.2479


%% plot resulting road network
for i=1:3
    DrawResultingNetwork(lineSet{i}, linelist_ace_ex{i}, labels_final_ace_ex{i}, img_plot{i});
end
%%

for i=3:3
    DrawLabelledLineList(img{i}, linelist_ace{i}, opt_labels_ace{i}); 
    DrawLabelledLineList(img{i}, linelist_ace_ex{i}, labels_final_ace_ex{i}); 
end

%% Junctions

n_ref_junctions{1} = 16;
n_missed_juctions{1} = 11;
n_false_junctions{1} = 1;
n_detected_junctions{1} = 6;

n_ref_junctions{2} = 22;
n_missed_juctions{2} = 5;
n_false_junctions{2} = 0;
n_detected_junctions{2} = 17;

n_ref_junctions{3} = 20;
n_missed_juctions{3} = 0;
n_false_junctions{3} = 1;
n_detected_junctions{3} = 20 + 1;

for i=1:3
    tmp = [detected_junctions_ace{i}.Position];
    tmp = [tmp(2:2:end);tmp(1:2:end)]';

    det_loc = tmp(1:2:end,:);
    ref_loc = tmp(2:2:end,:);

    rmse_junc{i} = sqrt(sum(((det_loc(:,1)-ref_loc(:,1)).^2 + (det_loc(:,2)-ref_loc(:,2)).^2))/size(det_loc,1));
    
    completeness_junc{i} = (n_ref_junctions{i} - n_missed_juctions{i}) / n_ref_junctions{i};
    correctness_junc{i} = (n_detected_junctions{i} - n_false_junctions{i}) / n_detected_junctions{i};
    
    quality_junc{i} = completeness_junc{i}*correctness_junc{i} /...
        (completeness_junc{i} - completeness_junc{i}*correctness_junc{i} + correctness_junc{i});
    
    fprintf('Junction Quality metrics(%s): completeness=%6.4f, correctness=%6.4f, quality=%6.4f, rmse=%6.4f\n', ...
        img_name{i}, completeness_junc{i},correctness_junc{i},quality_junc{i},rmse_junc{i});
    
end

% Junction Quality metrics(img1): completeness=0.3125, correctness=0.8333, quality=0.2941, rmse=2.4979
% Junction Quality metrics(img2): completeness=0.7727, correctness=1.0000, quality=0.7727, rmse=2.9842
% Junction Quality metrics(img3): completeness=1.0000, correctness=0.9524, quality=0.9524, rmse=2.6946

clear det_loc ref_loc

