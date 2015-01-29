function evi = get_evi(image)
    
    G = 2.5;
    L = 0.5;
    C1 = 6;
    C2 = 7.5;
    
    evi = G * (image(:,:,4) - image(:,:,1)) ./ ...
            (L + image(:,:,4) + C1 * image(:,:,1) - C2 * image(:,:,3));
  
end