function [void] = make_LTTPlotByState(S)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

times = length(S(1,:));
entries = length(S(:,1));

% %Plot true LTTs
% plot(1:times, S(1,:), 'c', 'LineWidth', 2.0);
% hold on;
% plot(1:times, S(2,:), 'm', 'LineWidth', 2.0);
% plot(1:times, S(3,:), 'r', 'LineWidth', 2.0);
% 
% %Plot sampled LTTs
% plot(1:times, S(4,:), 'c--');
% plot(1:times, S(5,:), 'm--');
% plot(1:times, S(6,:), 'r--');


%As subplots
%Plot true LTTs
subplot(3,1,1), plot(1:times, S(1,:), 'c', 'LineWidth', 2.0);
subplot(3,1,2), plot(1:times, S(2,:), 'm', 'LineWidth', 2.0);
subplot(3,1,3), plot(1:times, S(3,:), 'r', 'LineWidth', 2.0);

%Plot sampled LTTs
v1 = 4:3:entries;
v2 = 5:3:entries;
v3 = 6:3:entries;

meansState1 = mean(S(v1,:),1);
meansState2 = mean(S(v2,:),1);
meansState3 = mean(S(v3,:),1);

subplot(3,1,1), hold on, plot(1:times, meansState1, 'c--');
subplot(3,1,2), hold on, plot(1:times, meansState2, 'm--');
subplot(3,1,3), hold on, plot(1:times, meansState3, 'r--');

end

