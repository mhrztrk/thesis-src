function [linelist] = FindDetectedConnections(rdclseglist)

    lcnt = 0;
    linelist = [];
    nseg = size(rdclseglist,2); % # of segments
    for i=1:nseg
        nline = (size(rdclseglist{i},1) - 1); % # of lines in the segment 

        for j=1:nline
            lcnt = lcnt + 1;
            linelist(lcnt).s = rdclseglist{i}(j,:);                % start point
            linelist(lcnt).e = rdclseglist{i}(j+1,:);              % end point
            linelist(lcnt).c = (linelist(lcnt).s + linelist(lcnt).e)/2;  % center point
            linelist(lcnt).ang = atand((linelist(lcnt).e(2) - linelist(lcnt).s(2))/...
                                (linelist(lcnt).e(1) - linelist(lcnt).s(1))); % link orientation

            linelist(lcnt).nr = lcnt;
            linelist(lcnt).nconn = 0;   % # of links connected to node i.
            linelist(lcnt).conn = 0;    % id of the connected links. 
            linelist(lcnt).connpt = 0;    % connection point (0 = sp, 1 = ep). 

            linelist(lcnt).ncut = 0;
            linelist(lcnt).cutpt(1,:) = [0 0];
            linelist(lcnt).cutby(1) = 0; 

            fprintf('%2d - node(%2d)\n', i, lcnt);

            if(j > 1)
                linelist(lcnt).nconn = linelist(lcnt).nconn + 1;
                linelist(lcnt-1).nconn = linelist(lcnt-1).nconn + 1;

                linelist(lcnt).conn(linelist(lcnt).nconn) = linelist(lcnt-1).nr;
                linelist(lcnt-1).conn(linelist(lcnt-1).nconn) = linelist(lcnt).nr;

                linelist(lcnt).conn(linelist(lcnt).nconn) = linelist(lcnt-1).nr;
                linelist(lcnt-1).conn(linelist(lcnt-1).nconn) = linelist(lcnt).nr;

                linelist(lcnt).connpt(linelist(lcnt).nconn) = 0;    
                linelist(lcnt-1).connpt(linelist(lcnt-1).nconn) = 1;

            end
        end
    end

    for i=1:(lcnt-1)
        for j=(i+1):lcnt

            if(any(linelist(i).conn == j))
                continue;
            end

            if( all(linelist(i).s == linelist(j).s))
                linelist(i).nconn = linelist(i).nconn + 1;
                linelist(j).nconn = linelist(j).nconn + 1;

                linelist(i).conn(linelist(i).nconn) = j;
                linelist(j).conn(linelist(j).nconn) = i;

                linelist(i).connpt(linelist(i).nconn) = 0;
                linelist(j).connpt(linelist(j).nconn) = 0;

            elseif (all(linelist(i).s == linelist(j).e))
                linelist(i).nconn = linelist(i).nconn + 1;
                linelist(j).nconn = linelist(j).nconn + 1;

                linelist(i).conn(linelist(i).nconn) = j;
                linelist(j).conn(linelist(j).nconn) = i;

                linelist(i).connpt(linelist(i).nconn) = 0;
                linelist(j).connpt(linelist(j).nconn) = 1;

            elseif (all(linelist(i).e == linelist(j).s))
                linelist(i).nconn = linelist(i).nconn + 1;
                linelist(j).nconn = linelist(j).nconn + 1;

                linelist(i).conn(linelist(i).nconn) = j;
                linelist(j).conn(linelist(j).nconn) = i;

                linelist(i).connpt(linelist(i).nconn) = 1;
                linelist(j).connpt(linelist(j).nconn) = 0;

            elseif (all(linelist(i).e == linelist(j).e))
                linelist(i).nconn = linelist(i).nconn + 1;
                linelist(j).nconn = linelist(j).nconn + 1;

                linelist(i).conn(linelist(i).nconn) = j;
                linelist(j).conn(linelist(j).nconn) = i;

                linelist(i).connpt(linelist(i).nconn) = 1;
                linelist(j).connpt(linelist(j).nconn) = 1;

            end
        end
    end
    
    [linelist] = UpdateLineLengthAngle(linelist);
    
end
