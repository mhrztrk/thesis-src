function mask = generate_mask(imRGB)

    gray = rgb2gray(imRGB);
    mask = im2bw(gray, 0);
    mask = imcomplement(mask);
    
end
