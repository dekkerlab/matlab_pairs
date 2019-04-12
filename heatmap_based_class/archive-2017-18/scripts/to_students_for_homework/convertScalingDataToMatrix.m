function[scaling_matrix] = convertScalingDataToMatrix(scaling_data)
  [nrows, ncols] = size(scaling_data);
  scaling_matrix = zeros(ncols, ncols);
  for i = 1:ncols
    for k = 1:(ncols - i + 1)
      scaling_matrix(i + (k-1) ,k) = scaling_data(i);
      scaling_matrix(k, i + (k-1)) = scaling_data(i);
      end
  end
end
