function ang = GetAngle(x, y)

    if x<0 && y<0
        ang = 180 + atand(double(abs(y)/abs(x)));
    elseif x<0 && y>=0
        ang = 180 - atand(double(    y /abs(x)));
    elseif x>=0 && y<0
        ang = 360 - atand(double(abs(y)/    x));
    else
        ang =       atand(double(    y /    x));
    end    
end