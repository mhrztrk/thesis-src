function [V] = SormCliquePotential2(linelist, labels, Clique, mrf)

    % find lines labelled as road 
    plist = Clique((labels(Clique)==1));

    Vmin = -1 * (-mrf.params.Kl*2 + mrf.params.Kc*sind(180));
    
    switch(size(plist,2))
        case 0
            V = 0;

        case 1  % only one connection labelled as road
            
            V = linelist(plist(1)).prob * (mrf.params.Ke ...
                - mrf.params.Kl*linelist(plist(1)).len);

        case 2  % two connection labelled as road 

            defl = GetDeflAngle(linelist(plist(1)), linelist(plist(2)));
            
            if(defl > 90)
                V = -linelist(plist(1)).prob * mrf.params.Kl * linelist(plist(1)).len + ...
                    -linelist(plist(2)).prob * mrf.params.Kl * linelist(plist(2)).len + ...
                   ((linelist(plist(1)).prob + linelist(plist(2)).prob) / 2) * mrf.params.Kc * sind(defl);
            else
                V = ((linelist(plist(1)).prob + linelist(plist(2)).prob) / 2) * mrf.params.Ki*2;
            end
         
        case 3  % Junctions

            V = mrf.params.Ki*3;
            
            if (mrf.params.allowjunc == 1)
                defl(1) = GetDeflAngle(linelist(plist(1)), linelist(plist(2)));
                defl(2) = GetDeflAngle(linelist(plist(1)), linelist(plist(3)));
                defl(3) = GetDeflAngle(linelist(plist(2)), linelist(plist(3)));

                [~,ind] = min([defl(1) defl(2) defl(3)]);

                if(ind == 3)
                    %line_m  = linelist(plist(1));
                    line_b1 = linelist(plist(2));
                    line_b2 = linelist(plist(3));
                elseif (ind == 2)
                    %line_m  = linelist(plist(2));
                    line_b1 = linelist(plist(1));
                    line_b2 = linelist(plist(3));
                elseif (ind == 1)
                    %line_m  = linelist(plist(3));
                    line_b1 = linelist(plist(1));
                    line_b2 = linelist(plist(2)); 
                end

                %Clq = line_b1.clq(line_b1.clq == line_b2.clq); 
                
                if ((line_b1.clq(1) == line_b2.clq(1)) || (line_b1.clq(1) == line_b2.clq(2)))
                    Clq = line_b1.clq(1); 
                elseif ((line_b1.clq(2) == line_b2.clq(1)) || (line_b1.clq(2) == line_b2.clq(2)))
                    Clq = line_b1.clq(2);
                else
                    error('wtf');
                end
                
                Clq_b1 = line_b1.clq(line_b1.clq ~= Clq);
                Clq_b2 = line_b2.clq(line_b2.clq ~= Clq);
                
                conns_b1 = ceil(find([linelist.clq] == Clq_b1)/2);
                conns_b1 = conns_b1(conns_b1 ~= line_b1.nr);
                
                conns_b2 = ceil(find([linelist.clq] == Clq_b2)/2);
                conns_b2 = conns_b2(conns_b2 ~= line_b2.nr);
                
                ind_b1 = find(labels(conns_b1)==1);
                ind_b2 = find(labels(conns_b2)==1);

                if((size(ind_b1,1) == 1) && (size(ind_b2,1) == 1))

%                     D = sqrt((pi - defl_angle(line_b1, linelist(conns_b1(ind_b1)))/180)^2 + ...
%                              (pi - defl_angle(line_b2, linelist(conns_b2(ind_b2)))/180)^2);

                    if((GetDeflAngle(line_b1, linelist(conns_b1(ind_b1))) >= 120) && ...
                            (GetDeflAngle(line_b2, linelist(conns_b2(ind_b2))) >= 120))
                        D = 0;
                    else
                        D = Inf;
                    end

                    defl = sort(defl, 'descend');
                    sindfl = sind(defl);
                    
                    if(defl(2)>90 && defl(end)>30 && D < 4)
                        V = -mrf.params.Kl*sum([linelist(plist).len]) + ...
                             mrf.params.Kc*(sindfl(1) + (1-(sindfl(2) + sindfl(3))/2));
                    end
                    
                end
            end
        otherwise   % more than 3 connection labelled as road
            V = mrf.params.Ki*size(plist,2);
    end
    
    V = V + Vmin; 
    
end



function defl = GetDeflAngle(line1, line2)

    if(all(line1.s == line2.s) || all(line1.e == line2.e))
        defl = abs(line1.ang - line2.ang);
    elseif(all(line1.s == line2.e) || all(line1.e == line2.s))
        defl = abs(180 - (line1.ang - line2.ang));
    else
        error('wtf!!!');
    end

    if(defl > 360)
        defl = defl - 360;
    end
    
    if(defl > 180)
        defl = 360 - defl;
    end
    
end

