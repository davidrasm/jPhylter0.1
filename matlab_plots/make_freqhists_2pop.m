function [] = make_freqhists_2pop(theta_samples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%clear all;
%load PHYLterThetaSamples041312;

%theta_samples = PHYLterThetaSamples';

burnin = 101;
estParams = 2;
MCMC_params.iterations = length(theta_samples(1,:));

% beta
[beta_summary] = quantile(theta_samples(1, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(1,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,1),bar(x,n,1),xlim([0.38,0.48])
xlabel('\beta','FontSize',12) 
ylabel('Posterior Density','FontSize',12)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([0.4287, 0.4287], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% beta_b
[beta_summary] = quantile(theta_samples(2, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(2,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,2),bar(x,n,1),xlim([0,0.5])
xlabel('\rho','FontSize',12) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([0.2, 0.2], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

end
