function linelist = AddNewLine(linelist, ln1nr, ln2nr, ln1connpt, ln2connpt)

    nline = size(linelist, 2);
    
    nline = nline + 1;
    linelist(nline).nr = nline;
    
    if(ln1connpt == 0)
        linelist(nline).s = linelist(ln1nr).s;
    else
        linelist(nline).s = linelist(ln1nr).e;
    end
    
    if(ln2connpt == 0)
        linelist(nline).e = linelist(ln2nr).s;
    else
        linelist(nline).e = linelist(ln2nr).e;
    end
    
    linelist(nline).c = (linelist(nline).s + linelist(nline).e)/2;
    
    linelist(nline).nconn = 0;
    
    linelist(nline).nconn = linelist(nline).nconn + 1;
    linelist(nline).conn(linelist(nline).nconn) = ln1nr;
    linelist(nline).connpt(linelist(nline).nconn) = 0;
    
    linelist(nline).nconn = linelist(nline).nconn + 1;
    linelist(nline).conn(linelist(nline).nconn) = ln2nr;
    linelist(nline).connpt(linelist(nline).nconn) = 1;

    if(linelist(ln1nr).nconn > 0)
        conns_line1 = linelist(ln1nr).conn(linelist(ln1nr).connpt == ln1connpt);

        for i=1:size(conns_line1,2)
            linelist(conns_line1(i)).nconn = linelist(conns_line1(i)).nconn + 1;
            linelist(conns_line1(i)).conn(linelist(conns_line1(i)).nconn) = nline;
            linelist(conns_line1(i)).connpt(linelist(conns_line1(i)).nconn) = ...
                linelist(conns_line1(i)).connpt(linelist(conns_line1(i)).conn == ln1nr);

            linelist(nline).nconn = linelist(nline).nconn + 1;
            linelist(nline).conn(linelist(nline).nconn) = conns_line1(i);
            linelist(nline).connpt(linelist(nline).nconn) = 0;
        end
    end
    
    if(linelist(ln2nr).nconn > 0)
        conns_line2 = linelist(ln2nr).conn(linelist(ln2nr).connpt == ln2connpt);

        for i=1:size(conns_line2,2)
            linelist(conns_line2(i)).nconn = linelist(conns_line2(i)).nconn + 1;
            linelist(conns_line2(i)).conn(linelist(conns_line2(i)).nconn) = nline;
            linelist(conns_line2(i)).connpt(linelist(conns_line2(i)).nconn) = ...
                linelist(conns_line2(i)).connpt(linelist(conns_line2(i)).conn == ln2nr);

            linelist(nline).nconn = linelist(nline).nconn + 1;
            linelist(nline).conn(linelist(nline).nconn) = conns_line2(i);
            linelist(nline).connpt(linelist(nline).nconn) = 1;
        end
    end
    
    linelist(ln1nr).nconn = linelist(ln1nr).nconn + 1;
    linelist(ln1nr).conn(linelist(ln1nr).nconn) = nline;
    linelist(ln1nr).connpt(linelist(ln1nr).nconn) = ln1connpt;

    linelist(ln2nr).nconn = linelist(ln2nr).nconn + 1;
    linelist(ln2nr).conn(linelist(ln2nr).nconn) = nline;
    linelist(ln2nr).connpt(linelist(ln2nr).nconn) = ln2connpt;    
    
end