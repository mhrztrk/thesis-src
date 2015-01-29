function savi = get_savi(image)
    
    L = 0.5;

    savi = (image(:,:,4) - image(:,:,1)) * (1 + L)./((image(:,:,4) + image(:,:,1) + L));
  
end