function [] = make_prevplot(X)
    
burnin = 102;
time_indexes = 1:73;
serial_times = X(1,time_indexes);
iterations = length(X(:,1));
focalPopIndexes = burnin:2:iterations;
for i = 1:length(time_indexes)
   index = time_indexes(i);
   upper_traj(i) = quantile(X(focalPopIndexes, index), .975);
   median_traj(i) = quantile(X(focalPopIndexes, index), .5);
   lower_traj(i) = quantile(X(focalPopIndexes, index), .025);
end
plot(serial_times, median_traj, '-k','LineWidth', 2)
hold on
plot(serial_times, upper_traj, '--r','LineWidth', 2)
plot(serial_times, lower_traj, '--r','LineWidth', 2)
set(gca, 'XTick', serial_times);
datetick('x','yyyy');
%n_times = length(prevalence);
%plot(prevalence(1:(1/(MCMC_params.dt*10)):n_times-1), 'b', 'LineWidth', 2)
ylabel('Prevalence','FontSize',14)
xlabel('Year', 'FontSize', 14)
    
    
    