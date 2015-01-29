function [x, y, intersect] = FindTwoSegmentsIntersection(sp1,ep1,sp2,ep2)

    intersect = 0;
    tolerans = 0.001;
    
    % find lines intersection
    % line equation y = Mx + N
    dx = sp1(1)-ep1(1);
    if(dx)
        M1 = (sp1(2)-ep1(2))/dx;
        N1 = (sp1(1)*ep1(2) - sp1(2)*ep1(1))/dx;
        
        dx = sp2(1)-ep2(1);
        if(dx)
            M2 = (sp2(2)-ep2(2))/dx;
            N2 = (sp2(1)*ep2(2) - sp2(2)*ep2(1))/dx; 
            
            if(M1 == M2)
                x = -1;
                y = -1;
                intersect = 0;
                %disp('parallel lines!!!\n');  
                return;
            end
                
            x = (N2 - N1) / (M1 - M2);
            y = M1*x + N1;
        else
            x = sp2(1);
            y = M1*x + N1;
        end
        
    else
        
        dx = sp2(1)-ep2(1);
        if(dx)
            M2 = (sp2(2)-ep2(2))/dx;
            N2 = (sp2(1)*ep2(2) - sp2(2)*ep2(1))/dx; 

            x = sp1(1);
            y = M2*x + N2;
        else
            x = -1;
            y = -1;
            intersect = 0;
            %disp('parallel lines!!!\n');
            return;
        end      
        
    end
    
    % determine if segments actually intersect
    d1_1 = dist(sp1, [x y]');
    d2_1 = dist(ep1, [x y]');
    len_1 = dist(sp1, ep1');

    d1_2 = dist(sp2, [x y]');
    d2_2 = dist(ep2, [x y]');
    len_2 = dist(sp2, ep2');

    if(((d1_1 + d2_1) <= len_1 + tolerans) && ((d1_2 + d2_2) <= len_2 + tolerans))
        intersect = 1;
    end
    
end