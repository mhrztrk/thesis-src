function  [bw] = cc_threshold(bw, ccmin, majmin)

    CC = bwconncomp(bw);

    stats = regionprops(CC, 'MajorAxisLength', 'MinorAxisLength');
    
    if(majmin)
        majAxLen = [stats.MajorAxisLength];
        minAxLen = [stats.MinorAxisLength];
        ratio = majAxLen ./ minAxLen;
    end
    
    numPixels = cellfun(@numel,CC.PixelIdxList);

    for i=1:size(numPixels,2)
        if(numPixels(i) <= ccmin)
            bw(CC.PixelIdxList{i}) = 0;
        elseif(majmin)
            if((minAxLen(i) > 3) || (ratio(i) < 15))
                bw(CC.PixelIdxList{i}) = 0;
            end
        end
    end 
    
end