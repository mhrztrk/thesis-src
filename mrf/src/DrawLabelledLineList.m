% DrawLabelledLineList 
%   DrawLabelledLineList(img,linelist,labels) 
%
%
function [] = DrawLabelledLineList(varargin)

    if (nargin <= 1)
        error('Invalid argument!');
    else
        img = varargin{1};
        linelist = varargin{2};

        if( nargin == 2 )
            labels = ones(size(linelist));
        else
            labels = varargin{3};
        end
    end

   figure; imshow(img,[]); hold on;

    for i=1:size(linelist,2)
        if(labels(i) == 1)
            line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
                'Color',[0 0.2 1],'LineWidth',4);
            %scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
            scatter(linelist(i).s(2), linelist(i).s(1), 'r','o','filled','LineWidth', 0.5);
            scatter(linelist(i).e(2), linelist(i).e(1), 'r','o','filled','LineWidth', 0.5);
        end
    end 

end