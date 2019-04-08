function[diagonal_means] = getScalingPlotData(mymatrix)
[nrows, ncols] = size(mymatrix);
diagonal_means = zeros(1, nrows);


for i = 1:nrows
    this_sum = 0;
    for k = 1:(ncols-i) + 1
        this_number = mymatrix(i + (k-1) ,k);
        this_sum = this_sum + this_number;
    end
    this_average = this_sum / ( nrows - i + 1);
    diagonal_means(i) = this_average
end
end
