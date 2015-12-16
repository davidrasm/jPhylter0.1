function [] = make_freqhists(theta_samples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%clear all;
%load PHYLterThetaSamples041312;

%theta_samples = PHYLterThetaSamples';

burnin = 1;
estParams = 6;
MCMC_params.iterations = length(theta_samples(1,:));

% beta_f
[beta_summary] = quantile(theta_samples(1, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(1,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,1),bar(x,n,1),xlim([0.2,0.5])
xlabel('\beta_f','FontSize',14) 
ylabel('Posterior Density','FontSize',10)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([0.48, 0.48], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% beta_b
[beta_summary] = quantile(theta_samples(2, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(2,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,2),bar(x,n,1),xlim([0,0.01])
xlabel('\beta_b','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([0.0048, 0.0048], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% beta_g
[beta_summary] = quantile(theta_samples(3, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(3,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,3),bar(x,n,1),xlim([0.1,0.4])
xlabel('\beta_g','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([0.24, 0.24], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% alpha
[beta_summary] = quantile(theta_samples(4, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(4,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,4),bar(x,n,1),xlim([0,0.15])
xlabel('\alpha','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
%line([1565, 1565], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% delta
[beta_summary] = quantile(theta_samples(5, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(5,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,5),bar(x,n,1),xlim([0.2,0.8])
xlabel('\delta','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')

% S init
[beta_summary] = quantile(theta_samples(6, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(6,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,6),bar(x,n,1),xlim([0.0,1.0])
xlabel('M','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')

end

