function defl = defl_angle(line1, line2)
    connpt_1 = line1.connpt(line1.conn == line2.nr);
    connpt_2 = line2.connpt(line2.conn == line1.nr);

    if(isempty(connpt_1) || isempty(connpt_2))
        error('no connection between line%d to line%d\n', line1.nr, line2.nr);
    end

    if(connpt_1 == connpt_2)
        defl = abs(line1.ang-line2.ang);
    else
        defl = abs((line1.ang - 180) - line2.ang);
    end

    if(defl >= 360)
        defl = defl - 360; 
    end

    defl = min(defl, 360 - defl);    
end