i = 8;

h = figure;imshow(img{i})
hold on;
ref{i} = [];

j = 0;
while (1)
    [x,y] = ginput(2);

    s = [x(1) y(1)];
    e = [x(2) y(2)];
    
    line([s(1) e(1)], [s(2) e(2)], 'Color',[0 0 1],'LineWidth',3);
    ref{i}(s(1), s(2), 'g','.','LineWidth', 3);
    ref{i}(e(1), e(2), 'r','.','LineWidth', 3);

    j = j + 1;
    ref{i}{j}.s = s;
    ref{i}{j}.e = e;
    
    w = waitforbuttonpress;

    if(w ~= 0)
        break;
    end
end

%%

for i=5:8
    figure;imshow(img{i})
    hold on
    for j=1:size(ref{i},2) 
        line([ref{i}{j}.s(1) ref{i}{j}.e(1)], [ref{i}{j}.s(2) ref{i}{j}.e(2)],'Color',[0 0 1],'LineWidth',3);
        scatter(ref{i}{j}.s(1), ref{i}{j}.s(2), 'r','.','LineWidth', 3);
        scatter(ref{i}{j}.e(1), ref{i}{j}.e(2), 'r','.','LineWidth', 3);
    end

%     if(~isempty(junc{i}))
%         for j=1:size(junc{i},1)    
%             scatter(junc{i}(j, 1), junc{i}(j, 2), 'g','x','LineWidth', 4);
%         end
%     end
% 
%     if(~isempty(endp{i}))
%         for j=1:size(endp{i},1)    
%             scatter(endp{i}(j, 1), endp{i}(j, 2), 'r','x','LineWidth', 4);
%         end
%     end
end

clear i j

%%
h = figure;imshow(img1)
i = 0;
lineSet = [];

%%
hline=gline();


%%
i = i + 1;
lineSet{i}.s = get(hline,'Xdata');
lineSet{i}.e = get(hline,'Ydata');

  %%      