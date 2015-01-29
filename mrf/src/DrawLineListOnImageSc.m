function DrawLineListOnImageSc(img, linelist, labels)

    % display image 
    figure, imshow(img, []), hold on;
  
    ls_min = 0; %min([linelist.prob]);
    ls_max = 1; %max([linelist.prob]);
    ls_dst = (ls_max - ls_min) / 2;
    ls_mid = ls_min + ls_dst;
    
    
    for i=1:size(linelist,2)
        if(labels(i) == 1)
            
            val = linelist(i).prob;
                
            if (val > ls_mid)
                b = (val - ls_mid) / ls_dst;
            else
                b = 0;
            end

            if (val < ls_mid)
                g = (val - ls_min) / ls_dst;
            else
                g = (ls_max - val) / ls_dst;
            end

            if (val < ls_mid)
                r = (ls_mid - val) / ls_dst;
            else
                r = 0;
            end
            
            line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],'Color',[r g b],'LineWidth',2);
            scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
            
        end
    end 
end