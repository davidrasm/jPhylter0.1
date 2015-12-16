function [void] = testKDE2()

% Example (Gaussian mixture with distant modes):
clear all
% generate a Gaussian mixture with distant modes
data=[randn(100,1), randn(100,1)/4;
randn(100,1)+18, randn(100,1);
randn(100,1)+15, randn(100,1)/2-18;];
% call the routine
[bandwidth,density,X,Y]=kde2d(data);
% plot the data and the density estimate
surf(X,Y,density,'LineStyle','none'), view([0,60])
colormap hot, hold on, alpha(.8)
set(gca, 'color', 'blue');
plot(data(:,80),data(:,2),'w.','MarkerSize',5) 

% Example (simple Gaussian mixture)
% clear all
% % generate a Gaussian mixture with distant modes
% data=[randn(500,2);
% randn(500,1)+3.5, randn(500,1);];
% % call the routine
% [bandwidth,density,X,Y]=kde2d(data);
% % plot the data and the density estimate
% contour3(X,Y,density,1), hold on
% plot(data(:,1),data(:,2),'r.','MarkerSize',5)

% Example (Sinusoidal density):
% clear all
% X=rand(1000,1); Y=sin(X*10*pi)+randn(size(X))/3; data=[X,Y];
% % apply routine
% [bandwidth,density,X,Y]=kde2d(data);
% % plot the data and the density estimate
% surf(X,Y,density,'LineStyle','none'), view([0,70])
% colormap hot, hold on, alpha(.8)
% set(gca, 'color', 'blue');
% plot(data(:,1),data(:,2),'w.','MarkerSize',5) 
