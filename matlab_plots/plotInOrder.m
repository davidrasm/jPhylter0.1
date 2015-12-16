function [ output_args ] = plotInOrder(x)

rows = 3:2:length(x);

for n = 1:length(rows)
    
    plot(x(1,:),x(rows(n),:));
    pause

end

