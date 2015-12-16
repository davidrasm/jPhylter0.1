function [] = make_freqhists_spatial(theta_samples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%clear all;
%load PHYLterThetaSamples041312;

%theta_samples = PHYLterThetaSamples';

burnin = 501;
MCMC_params.iterations = length(theta_samples(1,:));

% beta
[beta_summary] = quantile(theta_samples(1, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(1,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,1),bar(x,n,1),xlim([0.0,6.0])
xlabel('\beta','FontSize',14) 
ylabel('Posterior Density','FontSize',14)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
%line([1.0002, 1.0002], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% alpha
[beta_summary] = quantile(theta_samples(2, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(2,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,2),bar(x,n,1),xlim([0,4])
xlabel('Seasonality (\alpha)','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
%line([0.08, 0.08], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% Delta
[beta_summary] = quantile(theta_samples(3, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(3,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,3),bar(x,n,1),xlim([0.0,370])
xlabel('Seasonal phase (\delta)','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
%line([0.001, 0.001], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% M
[beta_summary] = quantile(theta_samples(4, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(4,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,4),bar(x,n,1),xlim([0,0.4])
xlabel('Migration (M)','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
%line([1565, 1565], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% Ige
[beta_summary] = quantile(theta_samples(5, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(5,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,5),bar(x,n,1),xlim([3000,7000])
xlabel('Effective I_g','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')


end

