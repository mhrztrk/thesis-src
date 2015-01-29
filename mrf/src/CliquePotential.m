function [V] = CliquePotential(linelist, labels, Clique, mrf)

    % find lines labelled as road 
    plist = Clique((labels(Clique)==1));
    
    switch(size(plist,2))
        case 0
            V = 0;

        case 1  % only one connection labelled as road
            
            V = mrf.params.Ke - mrf.params.Kl*linelist(plist(1)).len;

        case 2  % two connection labelled as road 

            defl = linelist(plist(1)).defl(linelist(plist(1)).conn == plist(2));

            if(defl > 90)
                V = -mrf.params.Kl*(linelist(plist(1)).len + linelist(plist(2)).len) + ...
                     mrf.params.Kc*sind(defl);

            else
                V = mrf.params.Ki*2;
            end
         
        case 3  % Junctions

            V = mrf.params.Ki*3;
            
            if (mrf.params.allowjunc == 1)
                defl(1) = defl_angle(linelist(plist(1)), linelist(plist(2)));
                defl(2) = defl_angle(linelist(plist(1)), linelist(plist(3)));
                defl(3) = defl_angle(linelist(plist(2)), linelist(plist(3)));

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

                connpt_b1 = abs(line_b1.connpt(line_b1.conn == line_b2.nr) - 1);
                connpt_b2 = abs(line_b2.connpt(line_b2.conn == line_b1.nr) - 1);

                conns_b1 = line_b1.conn((line_b1.connpt == connpt_b1));
                conns_b2 = line_b2.conn((line_b2.connpt == connpt_b2));

                ind_b1 = find(labels(conns_b1)==1);
                ind_b2 = find(labels(conns_b2)==1);

                if((size(ind_b1,1) == 1) && (size(ind_b2,1) == 1))

%                     D = sqrt((pi - defl_angle(line_b1, linelist(conns_b1(ind_b1)))/180)^2 + ...
%                              (pi - defl_angle(line_b2, linelist(conns_b2(ind_b2)))/180)^2);

                    if((defl_angle(line_b1, linelist(conns_b1(ind_b1))) >= 120) && ...
                            (defl_angle(line_b2, linelist(conns_b2(ind_b2))) >= 120))
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
    
end



