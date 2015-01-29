
CC{1} = CC_im17; 
CC{2} = CC_im27;
CC{3} = CC_img1;
CC{4} = CC_img2;
CC{5} = CC_img3;
CC{6} = CC_img4;

CM{1} = CM_im17; 
CM{2} = CM_im27;
CM{3} = CM_img1;
CM{4} = CM_img2;
CM{5} = CM_img3;
CM{6} = CM_img4;

clres{1} = clres_17;
clres{2} = clres_27;
clres{3} = clres_img1;
clres{4} = clres_img2;
clres{5} = clres_img3;
clres{6} = clres_img4;
 
ref{ 1} = im_17_mask;
ref{ 2} = im_27_mask;
ref{ 3} = img1_mask;
ref{ 4} = img2_mask;
ref{ 5} = img3_mask;
ref{ 6} = img4_mask;

img_name{ 1} = 'im17';
img_name{ 2} = 'im27';
img_name{ 3} = 'img1';
img_name{ 4} = 'img2';
img_name{ 5} = 'img3';
img_name{ 6} = 'img4';


%% MRF Relaxation / Simulated Annealing with Gibbs Sampler %%

    nLevel = 6; % # of temp levels
    step = 70;  % # of step in each level
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
    
    
for i=8:8
    
    % detect all possible connections
    [linelist{i}, Cliques{i}] = SORMNetworkConstruct(CC{i}, [], clres{i}, 45, 0.1);
    
    % Draw all possible connections
    % DrawLabelledLineList(clres{i}, linelist{i}, ones(size(linelist{i},2))); 

    % save connection info for later use
    save_linelist_info(linelist{i}, sprintf('linelist_info_%s.txt', img_name{i}));   
    
    % initial labeling / not important for stochastic relaxation
    labels = ones(1,size(linelist{i},2))';

    % optimization
    [labels_GlbMin{i}] = SAwithGibbsSampler(linelist{i}, labels, step, nLevel, Ti, mrf, Cliques{i}, 0);

    %opt_label = labels_GlbMin{i};
    
    %save(sprintf('2-opt_labels_%s.mat',img_name{i}), 'opt_label');
    
    % draw optimum network
    DrawLabelledLineList(clres{i}, linelist{i}, labels_GlbMin{i});  
    
end

%%
for i=5:8
    DrawLabelledLineList(img{i}, linelist{i}, opt_labels{i});      
end

%% Junction detection (Automatic) %% 

% 
for i=1:4
    [labelsExt{i}, LineSetEx{i}, CliqSetEx{i}] = ConnectTerminationPoints(linelist{i}, Cliques{i}, opt_labels{i}, mrf, CC{i}, clres{i});
    DrawLabelledLineList(clres{i}, LineSetEx{i}, labelsExt{i}); 
end

%% Junction detection (Manual) %%

remove_lines{1} = [];
remove_lines{2} = [313 117];

%%
for i=8:8
    
    tmp = [cursor_info{i}.Position];
    tmp = [tmp(2:2:end);tmp(1:2:end)]';
    line_s = tmp(1:2:end,:);
    line_e = tmp(2:2:end,:);

    if(i == 2)
        % img2 icin optimum noktaya manuel olarak set ettik. 
        labels_final{i} = labels_new; 
    else
        labels_final{i} = opt_labels{i};
    end

    if(~isempty(remove_lines{i}))
        labels_final{i}(remove_lines{i}) = zeros(size(remove_lines));   
    end
    
    for j=1:size(line_s,1)
        new_lines{i}(j).s = line_s(j,:);
        new_lines{i}(j).e = line_e(j,:);
    end

    [linelist_ex{i}, labels_final_ex{i}] = AddLinesToLineList(new_lines{i}, linelist{i}, labels_final{i}, CC{i}, img{i});

    DrawLabelledLineList(img{i}, linelist_ex{i}, labels_final_ex{i}); 

end

%%
for i=3:3
    
    tmp = [cursor_info{i}.Position];
    tmp = [tmp(2:2:end);tmp(1:2:end)]';
    line_s = tmp(1:2:end,:);
    line_e = tmp(2:2:end,:);
    
    for j=1:size(line_s,1)
        new_lines{i}(j).s = line_s(j,:);
        new_lines{i}(j).e = line_e(j,:);
    end

    [linelist_ex{i}, labels_final_ex{i}] = AddLinesToLineList(new_lines{i}, linelist{i}, opt_labels{i}, CC{i}, img{i});

    DrawLabelledLineList(ref{i}, linelist_ex{i}, labels_final_ex{i}); 

end

%%
for i=3:3

    DrawLabelledLineList(img_plot{i}, linelist_ex{i}, labels_final_ex{i}); 

end

%% Quality measures before post-processing
for i=5:8
    [completeness{i}, correctness{i}, rmse{i}] = NetworkQuality(lineSet{i}, linelist{i}, opt_labels{i}, img{i});

    quality{i} = completeness{i}*correctness{i} / (completeness{i} - completeness{i}*correctness{i} + correctness{i});
    
    fprintf('Quality metrics(%s): completeness=%6.4f, correctness=%6.4f, quality=%6.4f, rmse=%6.4f\n', ...
        img_name{i}, completeness{i},correctness{i},quality{i},rmse{i});
end

% Quality metrics(img1): completeness=0.8825, correctness=0.9231, quality=0.8221, rmse=1.0012
% Quality metrics(img2): completeness=0.8624, correctness=0.9718, quality=0.8413, rmse=1.0800
% Quality metrics(img3): completeness=0.9438, correctness=0.9932, quality=0.9378, rmse=1.0016
% Quality metrics(img4): completeness=0.8864, correctness=0.9857, quality=0.8751, rmse=1.2080


%% Quality measures after post-processing

for i=6:8
    [completeness{i}, correctness{i}, rmse{i}] = NetworkQuality(lineSet{i}, linelist_ex{i}, labels_final_ex{i}, img{i});

    quality{i} = completeness{i}*correctness{i} / (completeness{i} - completeness{i}*correctness{i} + correctness{i});
    
    fprintf('Quality metrics(%s): completeness=%6.4f, correctness=%6.4f, quality=%6.4f, rmse=%6.4f\n', ...
        img_name{i}, completeness{i},correctness{i},quality{i},rmse{i});
end


% Quality metrics(img1): completeness=0.9394, correctness=0.9091, quality=0.8587, rmse=1.0365
% Quality metrics(img2): completeness=0.9232, correctness=0.9720, quality=0.8992, rmse=1.0793
% Quality metrics(img3): completeness=0.9810, correctness=0.9904, quality=0.9717, rmse=1.0153
% Quality metrics(img4): completeness=0.9557, correctness=0.9858, quality=0.9427, rmse=1.2293

% Quality metrics(im14): completeness=0.9853, correctness=0.8586, quality=0.8477, rmse=1.2982
% Quality metrics(im17): completeness=0.8298, correctness=0.9079, quality=0.7654, rmse=1.3329
% Quality metrics(im21): completeness=0.9702, correctness=0.8955, quality=0.8716, rmse=1.2104
% Quality metrics(im27): completeness=0.6059, correctness=0.9491, quality=0.5868, rmse=1.3783

% missed junctions(img1) = 2/16 , false junctions(img1) = 3/16 , rmse = 


%% plot resulting road network (before post-processing)
for i=1:4
    DrawResultingNetwork(lineSet{i}, linelist{i}, opt_labels{i}, img{i});
end

%% plot resulting road network
for i=5:8
    DrawResultingNetwork(lineSet{i}, linelist_ex{i}, labels_final_ex{i}, zeros(1000,1000));
end

%% Junctions

n_ref_junctions{1} = 16;
n_missed_juctions{1} = 3;
n_false_junctions{1} = 4;
n_detected_junctions{1} = 18;

n_ref_junctions{2} = 22;
n_missed_juctions{2} = 5;
n_false_junctions{2} = 0;
n_detected_junctions{2} = 17;

n_ref_junctions{3} = 20;
n_missed_juctions{3} = 0;
n_false_junctions{3} = 0;
n_detected_junctions{3} = 20;

for i=3:3
    tmp = [detected_junctions{i}.Position];
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

% Junction Quality metrics(img1): completeness=0.8125, correctness=0.7778, quality=0.6594, rmse=3.5993
% Junction Quality metrics(img2): completeness=0.7727, correctness=1.0000, quality=0.7727, rmse=1.9251
% Junction Quality metrics(img3): completeness=1.0000, correctness=1.0000, quality=1.0000, rmse=2.6110

clear det_loc ref_loc

%% plot reference road map %%

for i=3:3
    
    figure;imshow(img{i})
    hold on

    for j=1:size(lineSet{i},2)

        line([lineSet{i}{j}.s(1) lineSet{i}{j}.e(1)], [lineSet{i}{j}.s(2) lineSet{i}{j}.e(2)],'Color',[0 1 0],'LineWidth',3);

        if(lineSet{i}{j}.conn_status(1) == 0)
            scatter(lineSet{i}{j}.s(1), lineSet{i}{j}.s(2), 'r','.','LineWidth', 3);
        else
            scatter(lineSet{i}{j}.s(1), lineSet{i}{j}.s(2), 'g','.','LineWidth', 3);
        end

        if(lineSet{i}{j}.conn_status(2) == 0)
            scatter(lineSet{i}{j}.e(1), lineSet{i}{j}.e(2), 'r','.','LineWidth', 3);
        else
            scatter(lineSet{i}{j}.e(1), lineSet{i}{j}.e(2), 'g','.','LineWidth', 3);
        end

    end
    
    scatter(ref_junc{i}(:,2), ref_junc{i}(:,1),'.','r', 'LineWidth', 4);
    
end

%%
