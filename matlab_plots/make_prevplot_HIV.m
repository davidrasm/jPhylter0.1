function [] = make_prevplot_HIV()
    
    
load PISTest_xTrajs_HIVStoch_tree5; %load MCMC samples
X = PISTest_xTrajs_HIVStoch_tree5;

burnin = 2; %plus one because first row entry is time
iterations = length(X(:,1));
var1_indexes = burnin:5:iterations;
var2_indexes = (burnin+1):5:iterations;
var3_indexes = (burnin+2):5:iterations;
var4_indexes = (burnin+3):5:iterations;
time_indexes = X(1,:);
serial_times = (X(1,:)) ./ 365.25; %rescale time axes to years and shift so origin is year 0

%Transform cumulative incidence back to incidence
I = diff(X(var4_indexes,:),1,2);
[rI, cI] = size(I);
I = [zeros(rI,1), I];

%Compute quantiles for trajectories
for i = 1:length(time_indexes)
   upper_traj_var1(i) = quantile(X(var1_indexes, i), .975);
   median_traj_var1(i) = quantile(X(var1_indexes, i), .5);
   lower_traj_var1(i) = quantile(X(var1_indexes, i), .025);
   
   upper_traj_var2(i) = quantile(X(var2_indexes, i), .975);
   median_traj_var2(i) = quantile(X(var2_indexes, i), .5);
   lower_traj_var2(i) = quantile(X(var2_indexes, i), .025);
   
   upper_traj_var3(i) = quantile(X(var3_indexes, i), .975);
   median_traj_var3(i) = quantile(X(var3_indexes, i), .5);
   lower_traj_var3(i) = quantile(X(var3_indexes, i), .025);
   
%    upper_traj_var4(i) = quantile(X(var4_indexes, i), .975);
%    median_traj_var4(i) = quantile(X(var4_indexes, i), .5);
%    lower_traj_var4(i) = quantile(X(var4_indexes, i), .025);
   
   upper_traj_var4(i) = quantile(I(:,i), .975);
   median_traj_var4(i) = quantile(I(:, i), .5);
   lower_traj_var4(i) = quantile(I(:, i), .025);
   
end
%Plot var1 quantiles
plot(serial_times, median_traj_var1, '-c','LineWidth', 2)
%plot(serial_times, upper_traj_var1, '--c','LineWidth', 2)
hold on;
%plot(serial_times, lower_traj_var1, '--c','LineWidth', 2)
jbfill(serial_times, upper_traj_var1, lower_traj_var1, 'c', 'c', 1, 0.4);

%Plot var2 quantiles
hold on;
plot(serial_times, median_traj_var2, '-m','LineWidth', 2)
%plot(serial_times, upper_traj_var2, '--r','LineWidth', 2)
%plot(serial_times, lower_traj_var2, '--r','LineWidth', 2)
jbfill(serial_times, upper_traj_var2, lower_traj_var2, 'm', 'm', 1, 0.4);

%Plot var3 quantiles
hold on;
plot(serial_times, median_traj_var3, '-r','LineWidth', 2)
%plot(serial_times, upper_traj_var3, '--r','LineWidth', 2)
%plot(serial_times, lower_traj_var3, '--r','LineWidth', 2)
jbfill(serial_times, upper_traj_var3, lower_traj_var3, 'r', 'r', 1, 0.4);

box off;
ylabel('Prevalence','FontSize',14)
xlabel('Year', 'FontSize', 14)
%xlim([1970, 2012]);

%Plot incidence in figure 2
%figure(2);
%plot(serial_times(1:(end-1)), (median_traj_var4(1:(end-1))/28), 'k', 'LineWidth', 2.0)
%hold on;
%jbfill(serial_times(1:(end-1)), (upper_traj_var4(1:(end-1))/28), (lower_traj_var4(1:(end-1))/28), 'k', 'k', 1, 0.4);
%box off;
%ylabel('Incidence per day','FontSize',14)
%xlabel('Year', 'FontSize', 14)
%xlim([1970, 2012]);

%figure(3)
%totalPrevalence = median_traj_var1 + median_traj_var2 + median_traj_var3;
%plot(serial_times, totalPrevalence, 'k');

%Superimpose true prevalence
load PISTestPopData_HIVStoch_sim5.txt %load sim data
S = PISTestPopData_HIVStoch_sim5;
%simDataTimeVec = 0.0:7.0:(37*365.25);
%simDataTimeIndexes = 1:1:(length(S(2,:)-2)); %zeros(1,length(time_indexes));
%for i = 1:length(time_indexes)
    %simDataTimeIndexes(i) = find(simDataTimeVec == time_indexes(i)); 
%end
%Plot true prevalence
serial_times = (0.0:(7.0):(365.25*37))/365.25;
serial_times(end+1) = 37.0;
hold on;
plot(serial_times, S(2,:), 'c', 'LineWidth', 2);
plot(serial_times, S(3,:), 'm', 'LineWidth', 2);
plot(serial_times, S(4,:), 'r', 'LineWidth', 2);


%n_times = length(prevalence);
%plot(prevalence(1:(1/(MCMC_params.dt*10)):n_times-1), 'b', 'LineWidth', 2)
%ylabel('Prevalence','FontSize',14)
%xlabel('Year', 'FontSize', 14)

%Convert times on x-axis to calander times
%set(gca, 'XTick', serial_times);
%datetick('x','mmmyy');
