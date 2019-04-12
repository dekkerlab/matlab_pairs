matrix_file = '/Users/ozadamh/Documents/PSB_course/C-500000/C-500000_raw.txt';
raw_matrix = dlmread(matrix_file)

HeatMap( log(raw_matrix) , 'Colormap', hot )

filtered_matrix = raw_matrix(39:end , 39:end)

HeatMap( log(filtered_matrix) , 'Colormap', hot )

mymatrix = makeUpsideDownMatrix(filtered_matrix)

myscale = makeUpsideDownMatrix(colormap(hot))


HeatMap( log(mymatrix) , 'Colormap', myscale )



[nrows, ncols] = size(mymatrix);

