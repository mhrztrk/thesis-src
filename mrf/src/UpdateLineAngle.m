function [linelist] = UpdateLineLengthAngle(linelist)


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

end