function ndvi = get_vec_ndvi(vec)

    ndvi = zeros(size(vec,1), 1); 

    for i = 1:size(vec,1)
        ndvi(i,1) = (vec(i,4) - vec(i,1))/(vec(i,4) + vec(i,1));
    end

    ndvi = (ndvi + 1) * 0.5;

end