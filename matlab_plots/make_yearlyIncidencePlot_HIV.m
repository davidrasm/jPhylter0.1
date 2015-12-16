function [] = make_yearlyIncidencePlot_HIV()
    
    
load PISTest_xTrajs_DetroitMSMHIVStoch_tree8; %load MCMC samples
X = PISTest_xTrajs_DetroitMSMHIVStoch_tree8;

Y = load('eriksHIVIncidenceCounts.txt');

%Add zeros for 1973 and 1974
row1973 = [1973, 0];
row1974 = [1974, 0];
Y = [row1973; row1974; Y];
sampleDates = X(1,:);

nYears = 39;
yearSerialDates = zeros(1,nYears);
yearIncidenceCounts = zeros(1,nYears);
yearIncidenceEstimate = zeros(1,nYears);
%yearEarlyEstimate = zeros(1,nYears);
yearLowerEstimate = zeros(1,nYears);
yearUpperEstimate = zeros(1,nYears);
yearLength = zeros(1,nYears);

iterations = length(X(:,1));
var4_indexes = 5:5:iterations;
%var5_indexes = 6:5:iterations;

%Get data for each year from 1973 to 2011
for yr = 1:nYears;
    year = 1972+yr;
    yearString = strcat('Jan-01-', num2str(year));
    nextYearString = strcat('Jan-01-', num2str(year+1));
    serialDate = datenum(yearString);
    serialDateNext = datenum(nextYearString);
    yearSerialDates(yr) = serialDate;
    yearIncidenceCounts(yr) = Y(yr,2);
    
    startLoc = find(abs(sampleDates - serialDate) == min(abs(sampleDates-serialDate)));
    endLoc = find(abs(sampleDates - serialDateNext) == min(abs(sampleDates-serialDateNext)));
    yearLength(yr) = sampleDates(endLoc(1)) - sampleDates(startLoc(1));
    cumIStart = X(var4_indexes,startLoc(1));
    cumIEnd = X(var4_indexes,endLoc(1));
    diffs = cumIEnd - cumIStart;
    medianIncidence = median(diffs);
    yearIncidenceEstimate(yr) = medianIncidence;
    %cumIStart = X(var5_indexes,startLoc(1));
    %cumIEnd = X(var5_indexes,endLoc(1));
    %diffs = cumIEnd - cumIStart;
    %medianIncidence = median(diffs);
    %yearEarlyEstimate(yr) = medianIncidence;
    yearLowerEstimate(yr) = quantile(diffs, .975);
    yearUpperEstimate(yr) = quantile(diffs, .025);
end

plot(yearSerialDates, yearIncidenceCounts, 'k--'); hold on;
plot(yearSerialDates, yearIncidenceEstimate, 'r');
plot(yearSerialDates, yearLowerEstimate, 'r--');
plot(yearSerialDates, yearUpperEstimate, 'r--');
%plot(yearSerialDates, yearLength, 'b');
box off;
ylabel('Incidence','FontSize',14)
xlabel('Year', 'FontSize', 14)
datetick('x','yyyy');
%xlim([1970, 2012]);    


% burnin = 2; %plus one because first row entry is time
% iterations = length(X(:,1));
% var1_indexes = burnin:4:iterations;
% var2_indexes = (burnin+1):4:iterations;
% var3_indexes = (burnin+2):4:iterations;
% var4_indexes = (burnin+3):4:iterations;
% time_indexes = X(1,:);
% serial_times = (X(1,:)) ./ 365.25; %rescale time axes to years and shift so origin is year 0
% 
% %Transform cumulative incidence back to incidence
% I = diff(X(var4_indexes,:),1,2);
% [rI, cI] = size(I);
% I = [zeros(rI,1), I];
% 
% %Compute quantiles for trajectories
% for i = 1:length(time_indexes)
%    upper_traj_var1(i) = quantile(X(var1_indexes, i), .975);
%    median_traj_var1(i) = quantile(X(var1_indexes, i), .5);
%    lower_traj_var1(i) = quantile(X(var1_indexes, i), .025);
%    
%    upper_traj_var2(i) = quantile(X(var2_indexes, i), .975);
%    median_traj_var2(i) = quantile(X(var2_indexes, i), .5);
%    lower_traj_var2(i) = quantile(X(var2_indexes, i), .025);
%    
%    upper_traj_var3(i) = quantile(X(var3_indexes, i), .975);
%    median_traj_var3(i) = quantile(X(var3_indexes, i), .5);
%    lower_traj_var3(i) = quantile(X(var3_indexes, i), .025);
%    
% %    upper_traj_var4(i) = quantile(X(var4_indexes, i), .975);
% %    median_traj_var4(i) = quantile(X(var4_indexes, i), .5);
% %    lower_traj_var4(i) = quantile(X(var4_indexes, i), .025);
%    
%    upper_traj_var4(i) = quantile(I(:,i), .975);
%    median_traj_var4(i) = quantile(I(:, i), .5);
%    lower_traj_var4(i) = quantile(I(:, i), .025);
%    
% end
% %Plot var1 quantiles
% plot(serial_times, median_traj_var1, '-c','LineWidth', 2)
% %plot(serial_times, upper_traj_var1, '--c','LineWidth', 2)
% hold on;
% %plot(serial_times, lower_traj_var1, '--c','LineWidth', 2)
% jbfill(serial_times, upper_traj_var1, lower_traj_var1, 'c', 'c', 1, 0.4);
% 
% %Plot var2 quantiles
% hold on;
% plot(serial_times, median_traj_var2, '-m','LineWidth', 2)
% %plot(serial_times, upper_traj_var2, '--r','LineWidth', 2)
% %plot(serial_times, lower_traj_var2, '--r','LineWidth', 2)
% jbfill(serial_times, upper_traj_var2, lower_traj_var2, 'm', 'm', 1, 0.4);
% 
% %Plot var3 quantiles
% hold on;
% plot(serial_times, median_traj_var3, '-r','LineWidth', 2)
% %plot(serial_times, upper_traj_var3, '--r','LineWidth', 2)
% %plot(serial_times, lower_traj_var3, '--r','LineWidth', 2)
% jbfill(serial_times, upper_traj_var3, lower_traj_var3, 'r', 'r', 1, 0.4);
% 
% box off;
% ylabel('Prevalence','FontSize',14)
% xlabel('Year', 'FontSize', 14)
% xlim([1970, 2012]);
% 
% %Plot incidence in figure 2
% figure(2);
% plot(serial_times(1:(end-1)), (median_traj_var4(1:(end-1))/56), 'k', 'LineWidth', 2.0)
% hold on;
% jbfill(serial_times(1:(end-1)), (upper_traj_var4(1:(end-1))/56), (lower_traj_var4(1:(end-1))/56), 'k', 'k', 1, 0.4);
% box off;
% ylabel('Incidence per day','FontSize',14)
% xlabel('Year', 'FontSize', 14)
% xlim([1970, 2012]);

%figure(3)
%totalPrevalence = median_traj_var1 + median_traj_var2 + median_traj_var3;
%plot(serial_times, totalPrevalence, 'k');

%Superimpose true prevalence
%load HIVSimDataDNoise_tree2.txt %load sim data
%S = HIVSimDataDNoise_tree2;
%simDataTimeVec = 0.0:7.0:(35*365.25);
%simDataTimeIndexes = 1:1:(length(S(2,:)-2)); %zeros(1,length(time_indexes));
%for i = 1:length(time_indexes)
    %simDataTimeIndexes(i) = find(simDataTimeVec == time_indexes(i)); 
%end
%Plot true prevalence
%serial_times = (0.0:(7.0):(365.25*35))/365.25;
%serial_times(end+1) = 35.0;
%hold on;
%plot(serial_times, S(2,:), 'c', 'LineWidth', 2);
%plot(serial_times, S(3,:), 'm', 'LineWidth', 2);
%plot(serial_times, S(4,:), 'r', 'LineWidth', 2);


%n_times = length(prevalence);
%plot(prevalence(1:(1/(MCMC_params.dt*10)):n_times-1), 'b', 'LineWidth', 2)
%ylabel('Prevalence','FontSize',14)
%xlabel('Year', 'FontSize', 14)

%Convert times on x-axis to calander times
%set(gca, 'XTick', serial_times);
%datetick('x','mmmyy');

