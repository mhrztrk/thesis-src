
% state = 0 -> no connection
% state = 1 -> line1 & line2 connected
% state = 2 -> line1 divided
% state = 3 -> line2 divided
% state = 4 -> line1 & line2 connected over another line 
function [line1 line2 new_lines state] = LineMinDist(line1,line2)
    

    % parameters 
    th_dist = 60;   % max allowable distance to connect to line.
    th_angl = 90;  % min angle between lines to connect them 
    
    sp1 = line1.s;
    ep1 = line1.e;
    sp2 = line2.s;
    ep2 = line2.e;

    state = 0;

    
    d_ss = sqrt((sp1(1)-sp2(1))^2+(sp1(2)-sp2(2))^2);   % distance from sp1 to sp2
    d_se = sqrt((sp1(1)-ep2(1))^2+(sp1(2)-ep2(2))^2);   
    d_es = sqrt((ep1(1)-sp2(1))^2+(ep1(2)-sp2(2))^2);   
    d_ee = sqrt((ep1(1)-ep2(1))^2+(ep1(2)-ep2(2))^2);   
    
    if(min([d_ss d_se d_es d_ee])>th_dist)
        new_lines = [];
        return;
    elseif(min([d_ss d_se d_es d_ee])==0) % lines are connected  
        fprintf('line1(%d) and line2(%d) are connected!!!\n', line1.nr, line2.nr);
        
        line1.nconn = line1.nconn + 1;
        line1.conn(line1.nconn) = line2.nr;

        new_lines = [];
        return;
    end
    
    [x y] = FindIntersectionPoint(sp1,ep1,sp2,ep2); 
    if(isempty(x) || isempty(y))
        fprintf('line%d and line%d are parallel!\n', line1.nr, line2.nr);
        
            new_lines = [];
        return;    
        
    end
    
    
    len1 = sqrt((ep1(1)-sp1(1))^2+(ep1(2)-sp1(2))^2);
    d11 = sqrt((x-sp1(1))^2+(y-sp1(2))^2);
    d12 = sqrt((x-ep1(1))^2+(y-ep1(2))^2);
    
    len2 = sqrt((ep2(1)-sp2(1))^2+(ep2(2)-sp2(2))^2);
    d21 = sqrt((x-sp2(1))^2+(y-sp2(2))^2);
    d22 = sqrt((x-ep2(1))^2+(y-ep2(2))^2);
    
    if(len1 > max(d12,d11))
        % intersection point on line1
        fprintf('intersection point on line1(%d)\n',line1.nr);
        
        new_lines(1).s = sp1;
        new_lines(1).e = [x y]; % intersection point
        
        new_lines(2).s = [x y];
        new_lines(2).e = ep1; % intersection point
        
        new_lines(3).s = [x y];
        if(d21 < d22)
            new_lines(3).e = sp2; 
        else
            new_lines(3).e = sp1; 
        end
        
        state = 2; % delete line 1 
        
        if(isempty(line1.ncut))
            line1.ncut = 0;
        end
        line1.ncut = line1.ncut + 1;
        line1.cutpt(line1.ncut,:) = [x y];
        line1.cutby(line1.ncut,1) = line2.nr;
        
        if(d21 < d22)
            line1.cutby(line1.ncut,2) = 0;  % start point 
        else
            line1.cutby(line1.ncut,2) = 1;  % end point
        end
        
    elseif(len2 > max(d22,d21))
        % intersection point on line2
        fprintf('intersection point on line2(%d)\n',line2.nr);
       
        new_lines(1).s = [x y];
        if(d11 < d12)
            new_lines(1).e = sp1; 
        else
            new_lines(1).e = ep1; 
        end
       
        new_lines(2).s = [x y];
        new_lines(2).e = ep2; % intersection point
        
        new_lines(3).s = sp2;
        new_lines(3).e = [x y]; % intersection point
        
        state = 3; % delete line 2 
        
        line2.ncut = line2.ncut + 1;
        line2.cutpt(line2.ncut,1) = x;
        line2.cutpt(line2.ncut,2) = y;
        line2.cutby(line2.ncut,1) = line1.nr;
        
        if(d11 < d12)
            line2.cutby(line2.ncut,2) = 0;  % start point 
        else
            line2.cutby(line2.ncut,2) = 1;  % end point
        end
        
    else

        % if connection already exist, 
        nconn_pts = 0;
        conn_pts = [];
        conn_exist = 0;
        for m=1:line1.nconn
            for n=1:line2.nconn
                if(line1.conn(m) == line2.conn(n))
                    % there is already a link between line i and j
                    conn_exist = 1;
                    
                    fprintf('%d connected to %d through line %d\n', line1.nr, line2.nr, line1.conn(m));
                    nconn_pts = nconn_pts + 1;
                    conn_pts(nconn_pts,:) = [line1.connpt(m) line2.connpt(n)];                    
                end
            end
        end
        
        % these two line may be connected through another line.
        
        if(d11 < d12)
            ep1_new = ep1;
            sp1_new = sp1;
            dist_joint2conn1 = d11; 
        else
            ep1_new = sp1;
            sp1_new = ep1;
            dist_joint2conn1 = d12;
        end
       
        if(d21 < d22)
            ep2_new = ep2;
            sp2_new = sp2;
            dist_joint2conn2 = d21; 
        else
            ep2_new = sp2;
            sp2_new = ep2;
            dist_joint2conn2 = d22; 
        end

%         len1_new = sqrt((x - ep1_new(1))^2 + (y - ep1_new(2))^2);
%         len2_new = sqrt((x - ep2_new(1))^2 + (y - ep2_new(2))^2);
% 
%         dist_ep2ep = sqrt((ep1_new(1) - ep2_new(1))^2 + (ep1_new(2) - ep2_new(2))^2); 
% 
%         if(dist_ep2ep^2 < len1_new^2 + len2_new^2 - 2*len1_new*len2_new*cosd(th_angl))
%             state = 0; % no connection
%             new_lines = [];
%             return;
%         end
%                   

        dss = sqrt((sp1_new(1)-sp2_new(1))^2 + (sp1_new(2)-sp2_new(2))^2);

%         if(min(dist_joint2conn1, dist_joint2conn2) > ...
%                 max(dist_joint2conn1, dist_joint2conn2)*cosd(abs(line1.ang-line2.ang)))
            
        if(((dist_joint2conn1^2 + dss^2)>dist_joint2conn2^2) && ...
                ((dist_joint2conn2^2 + dss^2)>dist_joint2conn1^2) )  

            nline_sp = sp1_new;
            nline_ep = sp2_new;
            
            connpt = [0 1];
            
            line1.nconn = line1.nconn + 1; 
            line1.conn(line1.nconn) = -1;   % will be set later
            if(d11 < d12)
                line1.connpt(line1.nconn) = 0;  
            else
                line1.connpt(line1.nconn) = 1;
            end

            line2.nconn = line2.nconn + 1; 
            line2.conn(line2.nconn) = -1;   % will be set later
            if(d21 < d22)
                line2.connpt(line2.nconn) = 0;  
            else
                line2.connpt(line2.nconn) = 1;
            end
            
        elseif(dist_joint2conn1 > (dist_joint2conn2 + len2))
           
            if((max(d21,d22)^2 + min([d_ss d_se d_es d_ee])^2) < dist_joint2conn1^2)
            
                nline_sp = ep2_new;
                nline_ep = sp1_new;

                connpt = [1 0];

                line1.nconn = line1.nconn + 1; 
                line1.conn(line1.nconn) = -1;   % will be set later
                if(d11 < d12)
                    line1.connpt(line1.nconn) = 0;  
                else
                    line1.connpt(line1.nconn) = 1;
                end

                line2.nconn = line2.nconn + 1; 
                line2.conn(line2.nconn) = -1;   % will be set later
                if(d21 < d22)
                    line2.connpt(line2.nconn) = 1;  
                else
                    line2.connpt(line2.nconn) = 0;
                end   
            else
                new_lines = [];
                return;    
            end
            
        elseif(dist_joint2conn2 > (dist_joint2conn1 + len1))    
            
            if((max(d11,d12)^2 + min([d_ss d_se d_es d_ee])^2) < dist_joint2conn2^2)
            
                nline_sp = ep1_new;
                nline_ep = sp2_new;

                connpt = [0 1];

                line1.nconn = line1.nconn + 1; 
                line1.conn(line1.nconn) = -1;   % will be set later
                if(d11 < d12)
                    line1.connpt(line1.nconn) = 1;  
                else
                    line1.connpt(line1.nconn) = 0;
                end

                line2.nconn = line2.nconn + 1; 
                line2.conn(line2.nconn) = -1;   % will be set later
                if(d21 < d22)
                    line2.connpt(line2.nconn) = 0;  
                else
                    line2.connpt(line2.nconn) = 1;
                end  
            else
                new_lines = [];
                return;    
            end
            
        else
            % otherwise, no connection
            new_lines = [];
            return;
        end

        
        if(conn_exist == 1) 
            for m=1:nconn_pts
                if((conn_pts(m,1) == line1.connpt(line1.nconn)) && ...
                        (conn_pts(m,2) == line2.connpt(line2.nconn)))
                    line1.nconn = line1.nconn - 1;
                    line1.conn = line1.conn(1:line1.nconn);
                    line1.connpt = line1.connpt(1:line1.nconn);

                    line2.nconn = line2.nconn - 1;
                    line2.conn = line2.conn(1:line2.nconn);
                    line2.connpt = line2.connpt(1:line2.nconn);

                    % connection already exist
                    new_lines = [];
                    return;
                end
            end
        end
        
        
        % if length of new line larger than th_dist, remove new line
        if(sqrt((nline_sp(1) - nline_ep(1))^2 + (nline_sp(2) - nline_ep(2))^2) > th_dist)
            % new connection ignored
            line1.nconn = line1.nconn - 1;
            line1.conn = line1.conn(1:line1.nconn);
            line1.connpt = line1.connpt(1:line1.nconn);

            line2.nconn = line2.nconn - 1;
            line2.conn = line2.conn(1:line2.nconn);
            line2.connpt = line2.connpt(1:line2.nconn);
            new_lines = [];
            return;            
        end
        
        % create a new line structure
        new_lines(1).s = nline_sp; 
        new_lines(1).e = nline_ep; 
        
        new_lines(1).c = (new_lines(1).s + new_lines(1).e)/2;  % center point 
        new_lines(1).ang = atand((new_lines(1).e(2) - new_lines(1).s(2))/...
                            (new_lines(1).e(1) - new_lines(1).s(1))); % link orientation 
        
        new_lines(1).nr = -1;
        new_lines(1).nconn = 2;   % # of links connected to node i.
        new_lines(1).conn = [line1.nr line2.nr];    % id of the connected links. 
        new_lines(1).connpt = connpt;    % connection point (0 = sp, 1 = ep). 

        new_lines(1).ncut = 0;
        new_lines(1).cutpt(1,:) = [0 0];
        new_lines(1).cutby(1,:) = [0 0]; 
        
        
        state = 4; % line1 & line2 connected through another line
        
        disp('intersection point outside');   
    end
end


function [x y] = FindIntersectionPoint(sp1,ep1,sp2,ep2) 
    m1 = (sp1(2)-ep1(2))/(sp1(1)-ep1(1));
    m2 = (sp2(2)-ep2(2))/(sp2(1)-ep2(1));
   
    c1 = sp1(2) - m1*sp1(1);
    c2 = sp2(2) - m2*sp2(1);
    
    
    if(m1 == m2) % parallel lines
        x = [];
        y = [];
    elseif(isinf(m1))
        x = sp1(1);
        y = c2 + m2*x;
    elseif(isinf(m2)) 
        x = sp2(1);
        y = c1 + m1*x;
    else 
        x = (c1 - c2)/(m2 - m1);
        y = c1 + m1*x;
    end
    
end

