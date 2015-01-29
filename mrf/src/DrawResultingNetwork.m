function [completeness correctness rmse] = DrawResultingNetwork(refLine, extLine, labels, img)

% DrawLabelledLineList(img, extLine, labels);


%%

LineWidth = 3.5;

figure;imshow(img,[]); hold on;

matched_reference = cell(1, size(refLine,2));


for j=1:size(refLine,2)    
    %
    % refLine ile extLine'in [x,y] siralamasi farkli, duzeltmeyi simdilik
    % burada yapiyoruz.
    %
    ref_line.s = [refLine{j}.s(2) refLine{j}.s(1)];
    ref_line.e = [refLine{j}.e(2) refLine{j}.e(1)];
    
    pt = GenRectWindow(ref_line.s, ref_line.e, 5);
         
    %DrawRectWindow(pt);

    n_matched_reference = 0;
    
    for k=1:size(extLine,2)
        
        if(labels(k) == 1)
            
            [mline] = GetMatchingLine(extLine(k), ref_line, pt);
            
            if(~isempty(mline))
                n_matched_reference = n_matched_reference + 1;
                matched_reference{j}(n_matched_reference) = mline;
            end
            
            %line([extLine(k).s(2) extLine(k).e(2)], [extLine(k).s(1) extLine(k).e(1)],'Color',[1 0 0],'LineWidth',3);
                       
        end
    end     
    
    line([ref_line.s(2) ref_line.e(2)],[ref_line.s(1) ref_line.e(1)],'Color',[1 0 0],'LineWidth',LineWidth);
%     scatter(ref_line.s(2), ref_line.s(1), 'r','.','LineWidth', 4);
%     scatter(ref_line.e(2), ref_line.e(1), 'r','.','LineWidth', 4); 
    
    if(~isempty(matched_reference{j}))
       
        
        for k=1:size(matched_reference{j},2)
            
            mline = matched_reference{j}(k);
            
%             line([mline.s(2) mline.e(2)],[mline.s(1) mline.e(1)],'Color',[1 0 0], 'LineWidth',3);
%             scatter(mline.s(2), mline.s(1), 'b','o','LineWidth', 3);
%             scatter(mline.e(2), mline.e(1), 'b','o','LineWidth', 3);
            
            [prj1,prj2] = SegmentProjection(ref_line.s, ref_line.e, mline.s, mline.e);
     
            line([prj1(2) prj2(2)], [prj1(1) prj2(1)],'Color',[0 1 0],'LineWidth',LineWidth);
            
%             h = scatter(prj1(2), prj1(1), 'o', 'filled', 'LineWidth', 1);
%             set(h, 'MarkerFaceColor', [0 0.5 0]);
%             h = scatter(prj2(2), prj2(1), 'o', 'filled', 'LineWidth', 1);
%             set(h, 'MarkerFaceColor', [0 0.5 0]);
           
        end
        
    else
        fprintf('no matched ref %d\n',j);
    end
end

%%

% figure;imshow(img,[]); hold on;

matched_extraction = cell(1, size(extLine,2));

for j = 1:size(extLine,2)
        
    if(labels(j) == 1)
                
        extLine(j).matched = 0;
        
        pt = GenRectWindow(extLine(j).s, extLine(j).e, 5);

        %DrawRectWindow(pt);

        n_matched_extraction = 0;

        for k=1:size(refLine,2)

            %
            % refLine ile extLine'in [x,y] siralamasi farkli, duzeltmeyi simdilik
            % burada yapiyoruz.
            %
            ref_line.s = [refLine{k}.s(2) refLine{k}.s(1)];
            ref_line.e = [refLine{k}.e(2) refLine{k}.e(1)];

            [mline] = GetMatchingLine(ref_line, extLine(j), pt);

            if(~isempty(mline))
                n_matched_extraction = n_matched_extraction + 1;
                matched_extraction{j}(n_matched_extraction) = mline;
                extLine(j).matched = extLine(j).matched + 1;
                extLine(j).mline(extLine(j).matched) = mline;
            end

        end
        
        if(~isempty(matched_extraction{j}))

%             line([extLine(j).s(2) extLine(j).e(2)],[extLine(j).s(1) extLine(j).e(1)],'Color',[0 0.5 0],'LineWidth',3);

            for k=1:size(matched_extraction{j},2)

                mline = matched_extraction{j}(k);

%                 line([mline.s(2) mline.e(2)],[mline.s(1) mline.e(1)],'Color',[1 0 0], 'LineWidth',3);
%                 scatter(mline.s(2), mline.s(1), 'b','.','LineWidth', 3);
%                 scatter(mline.e(2), mline.e(1), 'b','.','LineWidth', 3);
                
            end

        else
            fprintf('no matched ext %d\n',j);
        end        
        
    end 
   
end

for j = 1:size(extLine,2)
    if(labels(j) == 1)
        if(extLine(j).matched == 0)
            line([extLine(j).s(2) extLine(j).e(2)],[extLine(j).s(1) extLine(j).e(1)],'Color',[0 0 1],'LineWidth',LineWidth);
            
            scatter(extLine(j).s(2), extLine(j).s(1), 'b','.','LineWidth', LineWidth);
            scatter(extLine(j).e(2), extLine(j).e(1), 'b','.','LineWidth', LineWidth);
                
        else
            pts = [];
            tot_len = 0;
            pts(1,:) = [extLine(j).s 0 0];
            pts(2,:) = [extLine(j).e 1 dist(extLine(j).e, extLine(j).s')];
            for k=1:extLine(j).matched
                tot_len = tot_len + dist(extLine(j).mline(k).s, extLine(j).mline(k).e');
                pts(2*k+1,:) = [extLine(j).mline(k).s 0 dist(extLine(j).mline(k).s, extLine(j).s')];
                pts(2*k+2,:) = [extLine(j).mline(k).e 1 dist(extLine(j).mline(k).e, extLine(j).s')];
            end
            pts = sortrows(pts,4);
            
            if(pts(2,4)- pts(1,4) > 5)
                line([pts(2,2) pts(1,2)],[pts(2,1) pts(1,1)],'Color',[0 0 1],'LineWidth',LineWidth);
                
                scatter(pts(2,2), pts(2,1), 'b','.','LineWidth', LineWidth);
                scatter(pts(1,2), pts(1,1), 'b','.','LineWidth', LineWidth);
            
            end
            
            if(pts(end,4)- pts(end-1,4) > 5)
                line([pts(end,2) pts(end-1,2)],[pts(end,1) pts(end-1,1)],'Color',[0 0 1],'LineWidth',LineWidth);
                
                scatter(pts(end,2), pts(end,1), 'b','.','LineWidth', LineWidth);
                scatter(pts(end-1,2), pts(end-1,1), 'b','.','LineWidth', LineWidth);
                
            end          
            
        end
    end
end

return;
%%
length_of_matched_extraction = 0;
length_of_matched_reference = 0;
total_length_of_reference = 0;
total_length_of_extraction = 0;


%% completeness

for j = 1:size(refLine,2)
    
    %
    % refLine ile extLine'in [x,y] siralamasi farkli, duzeltmeyi simdilik
    % burada yapiyoruz.
    %
    ref_line.s = [refLine{j}.s(2) refLine{j}.s(1)];
    ref_line.e = [refLine{j}.e(2) refLine{j}.e(1)];
    
    total_length_of_reference = total_length_of_reference + dist(ref_line.s, ref_line.e');
    
    if(~isempty(matched_reference{j}))
        for k=1:size(matched_reference{j},2)
            mline = matched_reference{j}(k);
            length_of_matched_reference = length_of_matched_reference + dist(mline.s, mline.e');
        end
    end    
end

completeness = length_of_matched_reference / total_length_of_reference;

%% correctness 

for j = 1:size(extLine,2)
    
    if(labels(j) == 1)
        
        loe = dist(extLine(j).s, extLine(j).e');
        total_length_of_extraction = total_length_of_extraction + loe;
        

        if(~isempty(matched_extraction{j}))
            lome = 0; 
            
            for k=1:size(matched_extraction{j},2)
                mline = matched_extraction{j}(k);
                lome = lome + dist(mline.s, mline.e');
            end
            
            % donemclerde matched segment boyutu reference segmentten buyuk
            % olabiliyor. Bu durum correctness'i oldugunda daha buyuk cikartiyor.
            % Simdilik buyuk oldugu durumlarda matched
            % extraction'i reference segment boyutuna esitliyoruz.
            if(lome > loe)
                lome = loe;
            end
            
            length_of_matched_extraction = length_of_matched_extraction + lome;
            
        end
    end
end

correctness = length_of_matched_extraction / total_length_of_extraction;

%% rmse
     
K = 0;
d_tot_squared = 0;

for j = 1:size(refLine,2)
    
    %
    % refLine ile extLine'in [x,y] siralamasi farkli, duzeltmeyi simdilik
    % burada yapiyoruz.
    %
    ref_line.s = [refLine{j}.s(2) refLine{j}.s(1)];
    ref_line.e = [refLine{j}.e(2) refLine{j}.e(1)];
    

    if(~isempty(matched_reference{j}))
        
        for k=1:size(matched_reference{j},2)
            
            K = K + 1;
            
            [prj1,prj2] = SegmentProjection(ref_line.s, ref_line.e, matched_reference{j}(k).s, matched_reference{j}(k).e);

            d1 = dist(matched_reference{j}(k).s, prj1');
            d2 = dist(matched_reference{j}(k).e, prj2');

            d_tot_squared = d_tot_squared + (d1+d2)/2;
            
        end
    end    
end

rmse = sqrt(d_tot_squared / K);

%%

end


function wind_pt = GenRectWindow(sp, ep, tickness)
    
    ang = GetAngle((sp(2) - ep(2)),(sp(1) - ep(1)));

    dx = tickness*cosd(ang);
    dy = tickness*sind(ang);

    % Find rectangle(i.e. tick line)'s corner points
    wind_pt(1,1) = sp(1) + dx;
    wind_pt(1,2) = sp(2) - dy;

    wind_pt(2,1) = sp(1) - dx;
    wind_pt(2,2) = sp(2) + dy;

    wind_pt(3,1) = ep(1) - dx;
    wind_pt(3,2) = ep(2) + dy;
    
    wind_pt(4,1) = ep(1) + dx;
    wind_pt(4,2) = ep(2) - dy;

end


function DrawRectWindow(pt)
    line([pt(1,2) pt(2,2)], [pt(1,1) pt(2,1)],'Color',[0 0.5 1],'LineWidth',LineWidth);
    line([pt(2,2) pt(3,2)], [pt(2,1) pt(3,1)],'Color',[0 0.5 1],'LineWidth',LineWidth);
    line([pt(3,2) pt(4,2)], [pt(3,1) pt(4,1)],'Color',[0 0.5 1],'LineWidth',LineWidth);
    line([pt(4,2) pt(1,2)], [pt(4,1) pt(1,1)],'Color',[0 0.5 1],'LineWidth',LineWidth);   
end



function [mline] = GetMatchingLine(extLine, refLine, wind)
    
    sp_inside = 0;
    ep_inside = 0;
    tolerans = 0.001; 

    mline = [];
    
    % check if sp of linelist is inside the buffer
    d1 = DistancePoint2Line(extLine.s, wind(1,:), wind(2,:)); 
    d2 = DistancePoint2Line(extLine.s, wind(3,:), wind(4,:)); 
    if((d1 + d2) <= dist(wind(2,:), wind(3,:)') + tolerans) 
        d1 = DistancePoint2Line(extLine.s, wind(2,:), wind(3,:)); 
        d2 = DistancePoint2Line(extLine.s, wind(4,:), wind(1,:));
        if((d1 + d2) <= dist(wind(1,:), wind(2,:)') + tolerans)
            sp_inside = 1;
        end
    end

    % check if ep of linelist is inside the buffer
    d1 = DistancePoint2Line(extLine.e, wind(1,:), wind(2,:)); 
    d2 = DistancePoint2Line(extLine.e, wind(3,:), wind(4,:)); 
    if((d1 + d2) <= dist(wind(2,:), wind(3,:)') + tolerans)
        d1 = DistancePoint2Line(extLine.e, wind(2,:), wind(3,:)); 
        d2 = DistancePoint2Line(extLine.e, wind(4,:), wind(1,:));
        if((d1 + d2) <= dist(wind(1,:), wind(2,:)') + tolerans)
            ep_inside = 1;
        end
    end

    if (sp_inside > 0 && ep_inside > 0)
        % both point inside
        % check if orientations are similar
        intersect_1 = extLine.s;
        intersect_2 = extLine.e;

    elseif (sp_inside > 0 || ep_inside > 0)

       if(sp_inside > 0) 
            intersect_1 = extLine.s;
        else
            intersect_1 = extLine.e;
        end

        % find intersection point
        matched_reference_ep_found = 0;
        for m=1:4
            ep1 = wind(m,:);
            ep2 = wind(mod(m,4)+1,:);
            [x, y, intersect] = FindTwoSegmentsIntersection(ep1,ep2,extLine.s,extLine.e);
            if(intersect)
                intersect_2 = [x y];
                matched_reference_ep_found = 1;
                break;
            end
        end

        if(~matched_reference_ep_found)
            error('wtf!!!');
        end

    else
        % endpoints of line not in the buffer zone
        n_intersect = 0;

        matched_reference_found = 0;

        for m=1:4
            ep1 = wind(m,:);
            ep2 = wind(mod(m,4)+1,:);
            [x, y, intersect] = FindTwoSegmentsIntersection(ep1,ep2,extLine.s,extLine.e);
            if(intersect)
                n_intersect = n_intersect + 1;
                if(n_intersect == 1)
                    intersect_1 = [x y];
                else 
                    intersect_2 = [x y];
                    matched_reference_found = 1;
                    break;
                end
            end
        end

        if(~matched_reference_found)
            return;
        end
    end

    [x, y, intersect] = FindTwoSegmentsIntersection(intersect_1, intersect_2, refLine.s, refLine.e);

    if(intersect)
        line_orientation       = GetAngle((x - refLine.e(1)), (y - refLine.e(2)));
        matched_reference_orientation = GetAngle((x - intersect_2(1))    , (y - intersect_2(2)));
    else
        % parallel lines
        line_orientation = GetAngle((refLine.s(1) - refLine.e(1)), (refLine.s(2) - refLine.e(2)));
        matched_reference_orientation = GetAngle((intersect_1(1) - intersect_2(1)), (intersect_1(2) - intersect_2(2)));        
    end
    
    defl = abs(line_orientation - matched_reference_orientation);
    if (defl > 270)
        defl = 360 - defl;
    elseif (defl > 180)
        defl = defl - 180;
    elseif (defl > 90)
        defl = 180 - defl;
    end

    if(defl < 30)
        % orientation difference is acceptable
        mline.s = intersect_1;
        mline.e = intersect_2;
    end
    
end
