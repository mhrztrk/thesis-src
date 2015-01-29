function [] = DrawTwoLabelsDiff(img, linelist, labels1, labels2)

    figure; imshow(img,[]); hold on;

    for i=1:size(linelist,2)
        if(labels1(i) == 1)
            if(labels2(i)==1)
                line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
                    'Color',[0 0 1],'LineWidth',2);
                scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
            else
                line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
                    'Color',[0 1 0],'LineWidth',2);
                scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');                
            end
            
        else
            if(labels2(i) == 1)
                line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
                    'Color',[1 0 0],'LineWidth',2);
                scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');            
            end
        end
    end 

end