function [] = make_prevplot_2PopDengue()
    
    
load DENV1_subMixed1130_xTrajs_structGlobalTied_060513; %load MCMC samples
allX = DENV1_subMixed1130_xTrajs_structGlobalTied_060513;

X = allX ;%(1:399,:);
%X = [X; allX(700:end,:)];

burnin = 2; %plus one because first row entry is time
iterations = length(X(:,1));
var1_indexes = burnin:2:iterations;
var2_indexes = (burnin+1):2:iterations;
time_indexes = X(1,:);
%serial_times = (X(1,:) - (16*365.25)) ./ 365.25; %rescale time axes to years and shift so origin is year 0

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
%subplot(2,1,1)
plot(time_indexes, median_traj_var1, 'r','LineWidth', 2)
%plot(serial_times, upper_traj_var1, '--b','LineWidth', 2)
hold on;
%plot(serial_times, lower_traj_var1, '--b','LineWidth', 2)
jbfill(time_indexes, upper_traj_var1, lower_traj_var1, 'r', 'r', 1, 0.4);


%Plot var2 quantiles
%subplot(2,1,2)
hold on;
plot(time_indexes, median_traj_var2, 'b','LineWidth', 2)
%plot(serial_times, upper_traj_var2, '--r','LineWidth', 2)
%plot(serial_times, lower_traj_var2, '--r','LineWidth', 2)
hold on;
jbfill(time_indexes, upper_traj_var2, lower_traj_var2, 'b', 'b', 1, 0.4);


ylabel('Monthly Incidence','FontSize',14)
xlabel('Year', 'FontSize', 14)

%Convert times on x-axis to calander times
set(gca, 'XTick', time_indexes);
datetick('x','yyyy');

figure(2)
plot(time_indexes, (median_traj_var1/10e06), 'r','LineWidth', 2)
hold on;
plot(time_indexes, (median_traj_var2/25e06), 'b','LineWidth', 2)
set(gca, 'XTick', time_indexes);
datetick('x','yyyy');