function [linelist] = FindPossibleConnections(linelist)

    lcnt = size(linelist,2);

    plcnt = 1;
    for i=1:(lcnt-1)
        for j=(i+1):lcnt  
            already_connected = 0;
            for m=1:linelist(i).nconn
                if(linelist(i).conn(m) == j)
                    % these two nodes already connected
                    fprintf('link%d and link%d are connected!\n', i, j);
                    already_connected = 1;
                end
                % if these
            end

            if(already_connected)
                continue;
            end

            [linelist(i) linelist(j) new_lines state] = LineMinDist(linelist(i),linelist(j));

            if(state == 0)      % no connection

            elseif(state == 1)  % Li & Lj directly connected

            elseif(state == 2)  % Li divided into two line
                % Cut point = {from i to j, from [ix iy] to [jx jy]}

            elseif(state == 3)  % Lj divided into two line 

            elseif(state == 4)  % Li & Lj connected through another line.

                plcnt = plcnt + size(new_lines,2);

                line1 = linelist(i);
                line2 = linelist(j);

                line1_connpt = line1.connpt(end);
                line2_connpt = line2.connpt(end);

                line1.nconn  = line1.nconn - 1;
                line1.conn   = line1.conn(1:line1.nconn);
                line1.connpt = line1.connpt(1:line1.nconn);

                line2.nconn  = line2.nconn - 1;
                line2.conn   = line2.conn(1:line2.nconn);
                line2.connpt = line2.connpt(1:line2.nconn);

                linelist(line1.nr) = line1;
                linelist(line2.nr) = line2;

                linelist = AddNewLine(linelist, line1.nr, line2.nr, line1_connpt, ...
                                                                    line2_connpt);

                fprintf('line%d & line%d connected through line%d\n', line1.nr, line2.nr, linelist(end).nr);
            end

        end
    end

    [linelist] = UpdateLineLengthAngle(linelist);
    
end