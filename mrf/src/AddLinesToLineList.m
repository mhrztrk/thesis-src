% 
function [linelist_ex, labels_ex] = AddLinesToLineList(new_lines, linelist, labels, CC, img)


    % first find the new one
    node_new = zeros(size(new_lines,2),2);
    node_exist = node_new;

    % linelist_ex -> new created lines
    % linelist_label_change 
    labels_ex = labels;

    linelist_ex = linelist;
    nr_new_line = size(linelist,2);

    for j=1:size(new_lines,2)

         Ds = dist(CC,new_lines(j).s');
         [val_s, ~] = min(Ds);

         De = dist(CC,new_lines(j).e');
         [val_e, ~] = min(De);

         new_node_available = 1;

         if(val_e == 0 && val_s == 0)
            % her iki noktada varolan bir node'un uzerinde,
            % yeni node yaratmaya gerek yok  
            new_node_available = 0;    
         elseif(val_e == 0) 
            node_new(j,:) = new_lines(j).s;
            node_exist(j,:) = new_lines(j).e;
         elseif (val_s == 0)
            node_new(j,:) = new_lines(j).e;
            node_exist(j,:) = new_lines(j).s;
         else
            error('wtf');
         end

         if(new_node_available)
             % find which lines existing nodes belong.
             line_found = 0;
             for k=1:size(linelist,2)
                 if(labels(k)==1)
                    if (all(linelist(k).s == node_exist(j,:)))
                        node_belongto(j,:) = [k 0];
                        line_found = 1;
                    elseif (all(linelist(k).e == node_exist(j,:)))
                        node_belongto(j,:) = [k 1];
                        line_found = 1;
                    end
                 end
             end

             % find the line which is closest to new node.  
             for k=1:size(linelist,2)
                if(labels(k) ==  1)
                    d1 = dist(node_new(j,:),linelist(k).s');
                    d2 = dist(node_new(j,:),linelist(k).e');


                    len = dist(linelist(k).s, linelist(k).e');
                    s = (d1+d2+len)/2;
                    d = 2 * sqrt(s*(s-d1)*(s-d2)*(s-len))/len;

                    % check if projection on the line 
                    if((d1^2+len^2 < d2^2) || (d2^2+len^2 < d1^2))
                        % projection not on the line, take the minimum 
                        if(d1 < d2)
                            D(k) = d1;
                            intersection_point(k,:) = linelist(k).s;
                        else
                            D(k) = d1;
                            intersection_point(k,:) = linelist(k).e; 
                        end
                    else

                        k1 = sqrt(d1^2 - d^2);
                        k2 = sqrt(d2^2 - d^2);
                        ratio = k1/(k1+k2); 

                        int_pt_x = linelist(k).s(1) - ratio*(linelist(k).s(1) - linelist(k).e(1));
                        int_pt_y = linelist(k).s(2) - ratio*(linelist(k).s(2) - linelist(k).e(2));

                        intersection_point(k,:) = [int_pt_x int_pt_y];

                        D(k) = d; 
                    end
                else
                    D(k) = Inf;
                end

            end
            [~,ind] = min(D);
            node_closeto(j) = ind;
            node_intersection(j,:) = intersection_point(ind,:);

            nr_new_line = nr_new_line + 1;
            labels_ex(nr_new_line) = 1;
            
            if(line_found) 
                if(node_belongto(j,2) == 0)
                    linelist_ex(nr_new_line).s = linelist(node_belongto(j,1)).s;
                else
                    linelist_ex(nr_new_line).s = linelist(node_belongto(j,1)).e;
                end
            else
                % image boundary
                linelist_ex(nr_new_line).s = node_exist(j,:);
            end
            linelist_ex(nr_new_line).e = node_intersection(j,:);

         else    
             line_found = 0;
             % Check if there is such a line exist
             for k=1:size(linelist,2)

                 if((all(linelist(k).s == new_lines(j).s) && all(linelist(k).e == new_lines(j).e)) || ...
                         (all(linelist(k).s == new_lines(j).e) && all(linelist(k).e == new_lines(j).s)))
                     labels_ex(k) = 1; 
                     line_found = 1;
                     break;
                 end

             end

             if(line_found == 0)
                % create new line
                nr_new_line = nr_new_line + 1;
                labels_ex(nr_new_line) = 1;
                linelist_ex(nr_new_line).s = new_lines(j).s;
                linelist_ex(nr_new_line).e = new_lines(j).e;
             end

         end
    end

end