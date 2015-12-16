function [] = make_freqhists_2PatchFig(theta_samples)
%For six params in two-patch model
%Top row is for pop1 and bottom row is for pop2

%clear all;
%load PHYLterThetaSamples041312;

%theta_samples = PHYLterThetaSamples';

burnin = 1;
MCMC_params.iterations = length(theta_samples(1,:));

% beta
[beta_summary] = quantile(theta_samples(1, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(1,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,1),bar(x,n,1),xlim([0.40,0.46])
xlabel('\beta','FontSize',12) 
ylabel('Posterior Density','FontSize',12)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
trueVal = 0.4287;
line([trueVal, trueVal], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% rho
[beta_summary] = quantile(theta_samples(5, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(5,burnin:MCMC_params.iterations), 18);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,2),bar(x,n,1),xlim([0.0,0.1])
xlabel('\rho','FontSize',12) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
trueVal = 0.012;
line([trueVal, trueVal], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% alpha1
[beta_summary] = quantile(theta_samples(2, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(2,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,3),bar(x,n,1),xlim([0.00,0.20])
xlabel('\alpha','FontSize',12) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([0.08, 0.08], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% alpha2
% [beta_summary] = quantile(theta_samples(3, burnin:MCMC_params.iterations), [.025 .5 .975]);
% [n,x]=hist(theta_samples(3,burnin:MCMC_params.iterations), 14);
% n = n/(MCMC_params.iterations-burnin);
% subplot(2,3,5),bar(x,n,1),xlim([0.02,0.14])
% xlabel('\alpha_2','FontSize',14) 
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
% line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
% line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
% line([0.08, 0.08], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% Delta1
[beta_summary] = quantile(theta_samples(3, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(3,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,4),bar(x,n,1),xlim([0.0,1.0])
xlabel('\delta_1','FontSize',12) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([0.50, 0.50], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% Delta2 - recentered around 1
for n = 1:MCMC_params.iterations
    if theta_samples(4,n) <= 0.5
        theta_samples(4,n) = theta_samples(4,n) + 1;
    end
end
[beta_summary] = quantile(theta_samples(4, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(4,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,5,5),bar(x,n,1),xlim([0.5,1.5])
xlabel('\delta_2','FontSize',12) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
line([1.0, 1.0], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

end
