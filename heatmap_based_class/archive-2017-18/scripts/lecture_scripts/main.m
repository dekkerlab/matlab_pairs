matrix_file = '/Users/ozadamh/Documents/PSB_course/C-500000/C-500000_raw.txt';
raw_matrix = dlmread(matrix_file);

mycolormap = makeUpsideDownMatrix(colormap(hot));

%row_sums = sum(raw_matrix);
clipped_matrix = raw_matrix(39:215 ,39:215 );

upside_down_matrix = makeUpsideDownMatrix(clipped_matrix);

HeatMap( log(upside_down_matrix), 'Colormap', mycolormap );

[nrows, ncols] = size(clipped_matrix);

scaling_data = getScalingPlotData(clipped_matrix);
plot( log(1:nrows) , log(scaling_data) );

scaling_matrix = convertScalingDataToMatrix(scaling_data);
HeatMap( log(scaling_matrix), 'Colormap', mycolormap );

observed_over_expected_matrix = clipped_matrix ./ scaling_matrix;
HeatMap( log(observed_over_expected_matrix), 'Colormap', mycolormap );

correlation_matrix_pre = corr(observed_over_expected_matrix);
correlation_matrix = makeUpsideDownMatrix(correlation_matrix_pre);
HeatMap(correlation_matrix,  'Colormap', mycolormap);
