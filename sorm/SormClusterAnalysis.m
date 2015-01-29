function [] = SormClusterAnalysis(CM, CC)

    th_min = 10;
    
    % SM : Scatter Matrix 
    
    % Remove clusters which have pixels less than threshold th_min 
    for i=1:size(CM,1)
        N = size(CM{i},1);
        X = (CM{1}(:,1) - CC(1,1)); 
        Y = (CM{1}(:,2) - CC(1,2)); 
        
        m_x = sum(X(:))/N;
        m_y = sum(Y(:))/N;
        
        var_x = ((X-m_x)'*(X-m_x))/(N*N);
        var_y = ((Y-m_y)'*(Y-m_y))/(N*N);        
        
        cov_xy = ((X-m_x)'*(Y-m_y))/(N*N);
        
        C = [var_x cov_xy; cov_xy var_y];
        
        d = eig(C);
        
        R(i) = max(d)/min(d);
        
    end
    
    
end

    