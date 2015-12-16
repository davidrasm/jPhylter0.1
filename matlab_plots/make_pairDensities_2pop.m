function [void] = make_pairDensities_2pop()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

theta_samples = load('PISTest_params_2PopStochRhoH_tree1');
theta_samples = theta_samples';
burnin = 1;
estParams = 2;
iterations = length(theta_samples(1,:));

figure(1);

beta = theta_samples(1, 1:1:iterations);
rho = theta_samples(2, 1:1:iterations);
true_beta = beta(1);
true_rho = rho(1);

mSize = 6.0;

%betaELim = [0.25,2]; % was 1.25
%betaELimAlter = [0.5,1.1];
%betaCLim = [0.0,0.6];
%betaALim = [0.0,5]; %was 1.15
%betaALimAlter = [0.0,1];

%Plot smoothed joint densities
colormap hot, hold on, alpha(.8)

data = [beta', rho'];
[bandwidth,density,X,Y]=kde2d(data, 2^5);
contourf(X,Y,density,15, 'LineStyle', 'none') %xlim(betaCLim), ylim(betaELimAlter)
set(gca, 'color', 'k');
hold on; plot(true_beta, true_rho, 'ob', 'MarkerSize', mSize, 'MarkerFaceColor', 'b');