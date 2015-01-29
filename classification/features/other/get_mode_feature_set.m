function train_set = get_mode_feature_set(segm_image, feature_mask, labels)
    
    % m -> height
    % n -> width
    % p -> feature count(e.g. 3 for RGB)
    [m n p] = size(segm_image);
    train_set = zeros(m*n,p);
    rgn_cnt = max(max(labels));
    
    ind_road = 1;
    
    for i=1:rgn_cnt
        
        [x y] = find(labels == (i-1), 1, 'first');
        
        if(feature_mask(x,y)== 1)
            train_set(ind_road,:) = segm_image(x,y,:);
            ind_road = ind_road + 1;
        end
    end
    
    %resize
    train_set = train_set(1:(ind_road-1),:);
     
end