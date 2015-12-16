function [density] = make_jointDensity(data)

% generate a Gaussian mixture with distant modes
%data=[randn(100,1), randn(100,1)/4;
%randn(100,1)+18, randn(100,1);
%randn(100,1)+15, randn(100,1)/2-18;];

%call the routine
[bandwidth,density,X,Y]=kde2d(data, 2^5);
% plot the data and the density estimate
surf(X,Y,density,'LineStyle','none'), view([0,60])
colormap hot, hold on, alpha(.8)
set(gca, 'color', 'blue');
plot(data(:,1),data(:,2),'w.','MarkerSize',5) 

% plot contour
figure(2)
contour3(X,Y,density,15), hold on
plot(data(:,1),data(:,2),'r.','MarkerSize',5)

end

