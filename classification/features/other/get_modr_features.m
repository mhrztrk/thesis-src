function M = get_modr_features(I)

    [H W D] = size(I);

    % Modulation Ratios of MS bands
    M = zeros(H,W, D*(D-1)/2);

    k = 1;                
    for i=1:(D-1)
        for j=(i+1):D
            M(:,:,k) = (I(:,:,i)-I(:,:,j))./(I(:,:,i)+I(:,:,j));
            k = k +1;
        end
    end
        
end