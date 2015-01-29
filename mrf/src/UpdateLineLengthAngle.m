function [linelist] = UpdateLineLengthAngle(linelist)


    % calculate edge angles & lengths
    for i=1:size(linelist, 2)
        
        linelist(i).ang = GetAngle((linelist(i).e(1) - linelist(i).s(1)), ...
                                    (linelist(i).e(2) - linelist(i).s(2))); % link orientation   

        linelist(i).len = sqrt((linelist(i).s(1)-linelist(i).e(1))^2 + ...
                                (linelist(i).s(2)-linelist(i).e(2))^2);

        
    end

end