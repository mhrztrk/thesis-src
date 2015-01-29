function prob = GetLineProb(pmap, pt1, pt2)

    % calculate edge angles & lengths
        tick = 2;
        ang = GetAngle((pt1(1)-pt2(1)),(pt1(2)-pt2(2)));
        
        dx = tick*cosd(ang);
        dy = tick*sind(ang);

        % Find rectangle(i.e. tick line)'s corner points
        pt(1,1) = pt1(1) + dy;
        pt(1,2) = pt1(2) - dx;

        pt(2,1) = pt1(1) - dy;
        pt(2,2) = pt1(2) + dx;

        pt(4,1) = pt2(1) + dy;
        pt(4,2) = pt2(2) - dx;

        pt(3,1) = pt2(1) - dy;
        pt(3,2) = pt2(2) + dx;

        BW = poly2mask(pt(:,2),pt(:,1), size(pmap,1), size(pmap,1));
        prob = mean(pmap(BW==1));
        
end