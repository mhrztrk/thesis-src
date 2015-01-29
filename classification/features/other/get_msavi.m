function msavi = get_msavi(image)

    msavi = zeros(size(image,1), size(image,2)); 

    for i = 1:size(image,1)
        for j = 1:size(image,2)
            msavi(i,j) = (2*image(i,j,4) + 1 - sqrt((2*image(i,j,4)+1)*(2*image(i,j,4)+1) - 8*(image(i,j,4)-image(i,j,1)))/2);
        end
    end

end