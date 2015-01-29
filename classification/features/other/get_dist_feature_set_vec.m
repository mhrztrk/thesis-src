function train_set = get_dist_feature_set_vec(vec)
    
    % m -> length
    % p -> feature count(e.g. 3 for RGB)
    [m p] = size(vec);
    train_set = zeros(m,p);
    
    ind_road = 1;
        
    mins = min(vec);
    maxs = max(vec);
    means = mean(vec);
    tol = (maxs - mins)*0.005;
    
    for k=1:m
        pix = zeros(1,p);
        pix(1,:) = double(vec(k,:)); 
        d = min(abs(pix - mins),abs(pix - maxs));
        conds = d > tol;
        if(min(conds) == 1)
            train_set(ind_road,:) = vec(k,:);
            ind_road = ind_road + 1;
        end
    end

    
    %resize
    train_set = train_set(1:(ind_road-1),:);
     
end 