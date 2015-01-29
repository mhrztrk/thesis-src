function train_set = get_dist_feature_set(image, feature_mask, labels)
    
    % m -> height
    % n -> width
    % p -> feature count(e.g. 3 for RGB)
    [m n p] = size(image);
    train_set = zeros(m*n,p);
    rgn_cnt = max(max(labels));
    
    ind_road = 1;
    
    for i=1:(rgn_cnt+1)
        
        [x y] = find(labels == (i-1));
        
        if(feature_mask(x(1),y(1))== 1)
            data_set = zeros(size(x,1),p);
        
            for k=1:size(x,1)
                data_set(k, :) = image(x(k),y(k), :);
            end
        
            sdev = std(data_set) / 4;     
            means = mean(data_set);
            
            for k=1:size(x,1)
                pix = zeros(1,p);
                pix(1,:) = double(image(x(k),y(k),:)); 
                d = abs(pix - means);
                conds = d < sdev;
                if(min(conds) == 1)
                    train_set(ind_road,:) = image(x(k),y(k),:);
                    ind_road = ind_road + 1;
                end
            end
        end
    end
    
    %resize
    train_set = train_set(1:(ind_road-1),:);
     
end