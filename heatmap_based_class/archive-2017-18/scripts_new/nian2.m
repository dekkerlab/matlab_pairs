matrix_file = '/Users/ozadamh/Documents/PSB_course/C-500000/C-500000_raw.txt';
raw_matrix = dlmread(matrix_file)

HeatMap( log(raw_matrix) , 'Colormap', hot)

filipped_hot = makeUpsideDownMatrix(hot)

HeatMap( log(raw_matrix) , 'Colormap', filipped_hot)

raw_sums = sum(raw_matrix)

clipped_matrix = raw_matrix(39:end, 39:end)

HeatMap( log(clipped_matrix) , 'Colormap', filipped_hot)

flipped_matrix = makeUpsideDownMatrix(clipped_matrix)

HeatMap( log(flipped_matrix) , 'Colormap', filipped_hot)

scaling_data = getScalingPlotData(clipped_matrix)
plot(log(1:nrows), log(scaling_data) )



decay_matrix = convertScalingDataToMatrix(scaling_data)
flipped_decay_matrix = makeUpsideDownMatrix(decay_matrix)
HeatMap( log(flipped_decay_matrix) , 'Colormap', filipped_hot)

observed_over_expected = clipped_matrix ./ decay_matrix
flipped_observed_over_expected = makeUpsideDownMatrix(observed_over_expected)
HeatMap( log(flipped_observed_over_expected) , 'Colormap', filipped_hot)

correlation_matrix = corr(observed_over_expected)
correlation_matrix_flipped = makeUpsideDownMatrix(correlation_matrix)
HeatMap( (correlation_matrix_flipped) , 'Colormap', filipped_hot)

imshow(correlation_matrix_flipped);
edge_image = edge(correlation_matrix, 'canny');
imshow(edge_image)

edge_sums = sum(transpose(edge_image))
plot(edge_sums)