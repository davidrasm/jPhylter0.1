function [] = make_incplot(X)
    
burnin = 5 + (4*25);

time_indexes = 1:163;
times = length(time_indexes);
serial_times = X(1,time_indexes);
iterations = length(X(:,1));
focalPopIndexes = burnin:4:iterations;
focalN = length(focalPopIndexes);

%If incidence is data is cumulative
I = zeros(focalN,times);
for p = 1:focalN
    currIndex = focalPopIndexes(p);
    I(p,2:times) = diff(X(currIndex,:));
end
I = I ./ (7*12); %to get incidence per day

%If incidence data is cumulative
%I = zeros(iterations,times);
%for p = 2:iterations
    %I(p,2:times) = diff(X(p,:));
%end

for i = 1:length(time_indexes)
   index = time_indexes(i);
   upper_traj(i) = quantile(I(:, index), .975);
   median_traj(i) = quantile(I(:, index), .5);
   lower_traj(i) = quantile(I(:, index), .025);
end

bounds = 2:162;
plot(serial_times(bounds), median_traj(bounds), '-k','LineWidth', 2)
hold on
plot(serial_times(bounds), upper_traj(bounds), '--r','LineWidth', 2)
plot(serial_times(bounds), lower_traj(bounds), '--r','LineWidth', 2)
set(gca, 'XTick', serial_times);
datetick('x','yyyy');
%n_times = length(prevalence);
%plot(prevalence(1:(1/(MCMC_params.dt*10)):n_times-1), 'b', 'LineWidth', 2)
ylabel('Incidence per day','FontSize',14)
xlabel('Year', 'FontSize', 14)

