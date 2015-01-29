function ndvi = get_ndvi(image)

    ndvi = (image(:,:,4) - image(:,:,1))./(image(:,:,4) + image(:,:,1));

end