function [linelist] = ApplyLineCuttings(linelist)

    for i=1:size(linelist,2)
        if(linelist(i).ncut==1)

            dist_sc = sqrt((linelist(i).s(1)-linelist(i).cutpt(1,1))^2 + ...
                    (linelist(i).s(2)-linelist(i).cutpt(1,2))^2);

            dist_ec = sqrt((linelist(i).e(1)-linelist(i).cutpt(1,1))^2 + ...
                            (linelist(i).e(2)-linelist(i).cutpt(1,2))^2);

            if(min(dist_sc,dist_ec)>5)

                fprintf('line%d divided by line%d\n', linelist(i).nr, linelist(i).cutby(1));

                [linelist, newLineNr] = DivideLineIntoTwo(linelist, linelist(i).nr, ...
                                                        linelist(i).cutpt(1,:));

                linelist = AddNewLine(linelist, linelist(i).nr, linelist(i).cutby(1,1), ...
                                            1, linelist(i).cutby(1,2));

            else
                line1 = linelist(i);
                line2 = linelist(linelist(i).cutby(1,1));    

                dist_cp2sp = sqrt((line1.cutpt(1,1) - line1.s(1))^2 + ...
                                (line1.cutpt(1,2) - line1.s(2))^2);
                dist_cp2ep = sqrt((line1.cutpt(1,1) - line1.e(1))^2 + ...
                                (line1.cutpt(1,2) - line1.e(2))^2);


                % find all lines connecting line1 and line2  
                [Y X] = find((repmat(line1.conn, [size(line2.conn,2) 1]) == ...
                              repmat(line2.conn, [size(line1.conn,2) 1])') ==1);


                if(dist_cp2ep < dist_cp2sp)
                    % check if there is already a line connecting these
                    % endpoints
                    if(~isempty(X))
                        if(~any((line1.connpt(X) == 1) & (line2.connpt(Y) == line1.cutby(1,2))))
                            linelist = AddNewLine(linelist, line1.nr, line2.nr, 1, line1.cutby(1,2));
                        end
                    else
                        linelist = AddNewLine(linelist, line1.nr, line2.nr, 1, line1.cutby(1,2));
                    end
                else
                    if(~isempty(X))
                        if(~any((line1.connpt(X) == 0) & (line2.connpt(Y) == line1.cutby(1,2))))
                            linelist = AddNewLine(linelist, line1.nr, line2.nr, 0, line1.cutby(1,2));
                        end
                    else
                        linelist = AddNewLine(linelist, line1.nr, line2.nr, 0, line1.cutby(1,2));
                    end
                end

                fprintf('--line%d and line%d connected through line%d\n', ...
                                        line1.nr, line2.nr, linelist(end).nr);

            end
        end
        
        linelist(i).ncut = 0;
        
    end

    [linelist] = UpdateLineLengthAngle(linelist);
end