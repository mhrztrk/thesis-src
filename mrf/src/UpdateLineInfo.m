function [linelist] = UpdateLineInfo(linelist, Cliques, pmap)

    % calculate edge angles & lengths
    for i=1:size(linelist, 2)
    %     linelist(i).ang = atand((linelist(i).e(2) - linelist(i).s(2))/...
    %                                 (linelist(i).e(1) - linelist(i).s(1))); % link orientation   
    %          
        linelist(i).ang = GetAngle((linelist(i).e(1) - linelist(i).s(1)), ...
                                    (linelist(i).e(2) - linelist(i).s(2))); % link orientation   

        linelist(i).len = sqrt((linelist(i).s(1)-linelist(i).e(1))^2 + ...
                                (linelist(i).s(2)-linelist(i).e(2))^2);

%         [~,~,~,X,Y] = bresenham(pmap,...
%                      [linelist(i).s(2) linelist(i).s(1);linelist(i).e(2) linelist(i).e(1)],0);                        
%         ind = sub2ind(size(pmap), X, Y); 
% 
%         linelist(i).prob = mean(pmap(ind));
        
    end

    linelist = UpdateLineProbInfo(pmap, linelist);
    
    % normalize lengths
    max_len = max([linelist(:).len]);

    for i=1:size(linelist, 2)

        linelist(i).len = linelist(i).len / max_len;

    end

    
    
    for i=1:size(linelist, 2)
    
        linelist(i).nclq = 0;
        linelist(i).clq = [];

        for j=1:size(Cliques,1)
            for k=1:size(Cliques{j},2)
                if(Cliques{j}(k)==linelist(i).nr)
                    linelist(i).nclq = linelist(i).nclq + 1;
                    linelist(i).clq(linelist(i).nclq) = j;
                end
            end
        end
    end
    
end

function linelist = UpdateLineProbInfo(pmap, linelist)

    % calculate edge angles & lengths
    for i=1:size(linelist, 2)
        tick = 2;
        dx = tick*cosd(linelist(i).ang);
        dy = tick*sind(linelist(i).ang);

        % Find rectangle(i.e. tick line)'s corner points
        pt(1,1) = linelist(i).s(1) + dy;
        pt(1,2) = linelist(i).s(2) - dx;

        pt(2,1) = linelist(i).s(1) - dy;
        pt(2,2) = linelist(i).s(2) + dx;

        pt(4,1) = linelist(i).e(1) + dy;
        pt(4,2) = linelist(i).e(2) - dx;

        pt(3,1) = linelist(i).e(1) - dy;
        pt(3,2) = linelist(i).e(2) + dx;


        BW = poly2mask(pt(:,2),pt(:,1), size(pmap,1), size(pmap,1));
        linelist(i).prob = mean(pmap(BW==1));

        
    end
end

