function [void] = make_likeSlice()

load LikeSurfLikeMatrixM02S1000;
load LikeSurfParam2MatrixM02S1000;
likes = LikeSurfLikeMatrixM02S1000;
params = LikeSurfParam2MatrixM02S1000;
%load LikeSurfParam2Matrix2Stage;

plot(params, likes,  'k', 'LineWidth', 2)
line([0.2 0.2], [ylim], 'LineStyle','-', 'Color', 'r', 'LineWidth', 2)

end

