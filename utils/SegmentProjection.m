function [prj1,prj2] = SegmentProjection(L1ep1, L1ep2, L2ep1, L2ep2)

    vec1 =  L1ep2 - L1ep1;
    
    vec2 =  L2ep1 - L1ep1;
    prj1 = (dot(vec1,vec2)/norm(vec1)^2) * vec1 + L1ep1;
    
    vec2 =  L2ep2 - L1ep1;
    prj2 = (dot(vec1,vec2)/norm(vec1)^2) * vec1 + L1ep1;
    
    d1 = norm(prj1-L1ep1);
    d2 = norm(prj1-L1ep2);
    
    if(max(d1,d2) > norm(L1ep2-L1ep1))
        if (d1 > d2)
            prj1 = L1ep2;
        else
            prj1 = L1ep1;
        end
    end
    
    d1 = norm(prj2-L1ep1);
    d2 = norm(prj2-L1ep2);
    
    if(max(d1,d2) > norm(L1ep2-L1ep1))
        if (d1 > d2)
            prj2 = L1ep2;
        else
            prj2 = L1ep1;
        end
    end
    
end