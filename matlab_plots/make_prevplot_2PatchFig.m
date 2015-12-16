function [] = make_prevplot_2PatchFig()
    
    
load FigSimXTraj2PatchLowMAsync; %load MCMC samples
X = FigSimXTraj2PatchLowMAsync;
%load FigSimData2PatchLowMAsync.txt %load sim data
%S = FigSimData2PatchLowMAsync;

burnin = 501+1; %plus one because first row entry is time
iterations = length(X(:,1));
var1_indexes = burnin:2:iterations;
var2_indexes = (burnin+1):2:iterations;
time_indexes = X(1,:);
serial_times = (X(1,:) - (16*365.25)) ./ 365.25; %rescale time axes to years and shift so origin is year 0

%Compute quantiles for trajectories
for i = 1:length(time_indexes)
   upper_traj_var1(i) = quantile(X(var1_indexes, i), .975);
   median_traj_var1(i) = quantile(X(var1_indexes, i), .5);
   lower_traj_var1(i) = quantile(X(var1_indexes, i), .025);
   
   upper_traj_var2(i) = quantile(X(var2_indexes, i), .975);
   median_traj_var2(i) = quantile(X(var2_indexes, i), .5);
   lower_traj_var2(i) = quantile(X(var2_indexes, i), .025);
end
%Plot var1 quantiles
%plot(serial_times, median_traj_var1, '-k','LineWidth', 2)
%plot(serial_times, upper_traj_var1, '--b','LineWidth', 2)
%hold on;
%plot(serial_times, lower_traj_var1, '--b','LineWidth', 2)
jbfill(serial_times, upper_traj_var1, lower_traj_var1, 'b', 'b', 1, 0.4);

%Plot var2 quantiles
%plot(serial_times, median_traj_var2, '-k','LineWidth', 2)
%plot(serial_times, upper_traj_var2, '--r','LineWidth', 2)
%plot(serial_times, lower_traj_var2, '--r','LineWidth', 2)
jbfill(serial_times, upper_traj_var2, lower_traj_var2, 'r', 'r', 1, 0.4);

%Superimpose true prevalence
load FigSimData2PatchLowMAsync.txt %load sim data
S = FigSimData2PatchLowMAsync;
simDataTimeVec = (365.25*16.0):0.25:(365.25*20.0);
simDataTimeIndexes = zeros(1,length(time_indexes));
for i = 1:length(time_indexes)
    simDataTimeIndexes(i) = find(simDataTimeVec == time_indexes(i)); 
end
%Plot true prevalence
hold on;
plot(serial_times, S(3,simDataTimeIndexes), 'b', 'LineWidth', 2);
plot(serial_times, S(4,simDataTimeIndexes), 'r', 'LineWidth', 2);


%n_times = length(prevalence);
%plot(prevalence(1:(1/(MCMC_params.dt*10)):n_times-1), 'b', 'LineWidth', 2)
ylabel('Prevalence','FontSize',14)
xlabel('Year', 'FontSize', 14)

%Convert times on x-axis to calander times
%set(gca, 'XTick', serial_times);
%datetick('x','mmmyy');