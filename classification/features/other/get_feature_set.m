function train_set = get_feature_set(image, mask)
    
    % m -> height
    % n -> width
    % p -> feature count(e.g. 3 for RGB)
    [m n p] = size(image);
    train_set = zeros(m*n,p);
    
    ind_road = 1;
    
    for i=1:m
        for j=1:n
            if(mask(i,j)== 1)
                if(image(i,j,1) ~= 0 || image(i,j,2) ~= 0 || image(i,j,3) ~= 0)
                    train_set(ind_road,:) = image(i,j,:);
                    ind_road = ind_road + 1;
                end
            end
        end
    end
    
    %resize
    train_set = train_set(1:(ind_road-1),:);
     
end