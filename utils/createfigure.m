function createfigure(X1, Y1, S1, C1)
%CREATEFIGURE(X1,Y1,S1,C1)
%  X1:  scatter x
%  Y1:  scatter y
%  S1:  scatter s
%  C1:  scatter c

%  Auto-generated by MATLAB on 27-Oct-2012 16:02:46

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on');

xlim(axes1,[0 255]);
ylim(axes1,[0 255]);
hold(axes1,'all');

% Create scatter
scatter(b,n,'Marker','.');

