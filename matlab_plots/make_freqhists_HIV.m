function [void] = make_freqhists_HIV(theta_samples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

theta_samples = load('PISTest_params_HIVStoch_tree25');
theta_samples = theta_samples';
burnin = 201;
estParams = 3;
MCMC_params.iterations = length(theta_samples(1,:));

% betaChronic
theta_samples(1,:) = theta_samples(1,:) * 365.25; %rescale transmission rates to years
[beta_summary] = quantile(theta_samples(1, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(1,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,1),bar(x,n,1),xlim([0.0,0.4])
xlabel('\beta_C','FontSize',14) 
ylabel('Posterior Density','FontSize',12)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
betaC = 1.2943e-04 * 365.25;
line([betaC, betaC], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% betaEarly
theta_samples(2,:) = theta_samples(2,:)*365.25; %rescale transmission rates to years
[beta_summary] = quantile(theta_samples(2, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(2,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,2),bar(x,n,1),xlim([0.25,1.75])
xlabel('\beta_E','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
betaE = 20*betaC;
line([betaE, betaE], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% betaAIDS
theta_samples(3,:) = theta_samples(3,:) *365.25; %rescale transmission rates to years
[beta_summary] = quantile(theta_samples(3, burnin:MCMC_params.iterations), [.025 .5 .975]);
[n,x]=hist(theta_samples(3,burnin:MCMC_params.iterations), 14);
n = n/(MCMC_params.iterations-burnin);
subplot(1,estParams,3),bar(x,n,1),xlim([0.0,1.25])
xlabel('\beta_A','FontSize',14) 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
betaA = 5*betaC;
line([betaA, betaA], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% alpha
% [beta_summary] = quantile(theta_samples(4, burnin:MCMC_params.iterations), [.025 .5 .975]);
% [n,x]=hist(theta_samples(4,burnin:MCMC_params.iterations), 14);
% n = n/(MCMC_params.iterations-burnin);
% subplot(1,estParams,4),bar(x,n,1),xlim([0.0000,0.0015])
% xlabel('\alpha','FontSize',14) 
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
% line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
% line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
% alpha = 0.00015;
%line([alpha, alpha], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% % Init I
% [beta_summary] = quantile(theta_samples(5, burnin:MCMC_params.iterations), [.025 .5 .975]);
% [n,x]=hist(theta_samples(5,burnin:MCMC_params.iterations), 14);
% n = n/(MCMC_params.iterations-burnin);
% subplot(1,estParams,5),bar(x,n,1),xlim([0,5.0])
% xlabel('I_1_9_7_5','FontSize',14) 
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
% line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
% line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
% initI = 0.00015;
% line([alpha, alpha], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

% % Init time
% [beta_summary] = quantile(theta_samples(5, burnin:MCMC_params.iterations), [.025 .5 .975]);
% [n,x]=hist(theta_samples(5,burnin:MCMC_params.iterations), 14);
% n = n/(MCMC_params.iterations-burnin);
% subplot(1,estParams,5),bar(x,n,1),xlim([720990,723181])
% xlabel('InitTime','FontSize',14) 
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
% line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
% line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')
% datetick('x','yyyy');
%initI = 0.00015;
%line([alpha, alpha], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')

%R0
% figure(2)
% gammaE = 1 / (2.55 * 365.25);
% gammaC = 1 / (6.31 * 365.25);
% gammaA = 1 / (2.55 * 365.25);
% mu = 1 / (40.28*365.25);
% betaCTrue = 1.2943e-04;
% betaATrue = 5.0 * betaCTrue;
% eTerm = theta_samples(2,:) ./ (mu + gammaE);
% probSurviveE = gammaE / (mu + gammaE);
% cTerm = probSurviveE .* theta_samples(1,:) .* (1/(mu + gammaC));
% probSurviveC = gammaC / (mu + gammaC);
% aTerm = probSurviveE * probSurviveC .* (theta_samples(3,:)) .* (1/gammaA);
% %aTerm = probSurviveE * probSurviveC .* betaATrue .* (1/gammaA);
% estR0 = eTerm + cTerm + aTerm;
% hist(estR0);

end

