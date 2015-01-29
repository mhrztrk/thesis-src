function d = DistancePoint2Line(pt, line_ep1, line_ep2)

    d1 = dist(pt,line_ep1');
    d2 = dist(pt,line_ep2');

    len = dist(line_ep1, line_ep2');

    s = (d1+d2+len)/2;
    d = 2 * sqrt(s*(s-d1)*(s-d2)*(s-len))/len;

end