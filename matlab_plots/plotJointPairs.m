function [ output_args ] = plotJointPairs(data)
%UNTITLED Summary of this function goes here
%  Plot the joint pairwise densities of each pair of params as scatter
%  plots

nParams = length(data(1,:));

count = 1;
for i = 1:nParams
    for j = 1:nParams
        
        subplot(nParams, nParams, count), scatter(data(:,j),data(:,i));
        count = count + 1;
        
    end
end
        



end

