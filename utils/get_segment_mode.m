function [modes] = get_segment_mode(image, labels) 

    modes = zeros(rgn_cnt,3);
    for m = 1:rgn_cnt

        dim = size(find(labels == m), 1);
        hist = zeros(dim,3);

        k = 1;
        for i = 1:imH
            for j = 1:imW
                if (labels(i,j) == m)
                    hist(k,:) = image(i,j,:);
                    k = k + 1;
                end
            end
        end
        tmp = mode(hist);
        modes(m,:) = mean(tmp); 
    end
    
end