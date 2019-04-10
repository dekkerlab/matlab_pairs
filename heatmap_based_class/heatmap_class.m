% let's start with the previous homework
% uploading it here is sort of "strange", but
% if there are student who can browse the repository
% to find prepared homework and understand it
% it's already something - mostly good, I think


%  1.
% read pairs from a text file as a table.
filename = '../U54_HFF_plate_subset.txt';
table = readtable(filename,"Format","%s%d%d%s%s");
table.Properties.VariableNames = {'chrom','pos1','pos2','str1','str2'};
% check to see how it looks like
head(table)

% 2.
% define bins ..., e.g., bins=0:binsize:chromlen
% use "discretize" with pre-set bins on pos1 pos2 ...
% there is bins functions - sounds intriguing what does it do?
% they woudl turn into 'bin1' 'bin2'
% define the bins

% chromlen as the largest position ...
chromlen = max(max(table.pos1),max(table.pos2))
% round it up - to 1000-sands
% chromlen = round(chromlen,-3);
binsize=100000; % 100kb bins -- feel free to try different size bins
chromlen_in_bins = int32(ceil(double(chromlen)/binsize));
chromlen_rounded = chromlen_in_bins * binsize;
bin_edges=0:binsize:chromlen_rounded;

% Let bin pos1 into "bins" i.e.
% assign each position to corresponding bin range
% e.g: pos1=25000 will belong to bin_index1=3 and
% corresponding to bin_ranges1=(20000-30000)

bin_index1 = discretize(table.pos1,bin_edges);

bin_index2 = discretize(table.pos2,bin_edges);
% seems like a column already ...

% append column of bin index back to our initial table

table_binned = addvars(table,bin_index1,'After','pos1');
table_binned = addvars(table_binned,bin_index2,'After','pos2');

% 3.

num_bins=length(bin_edges)-1; % compare with max(bin_index1)

% then initialize empty/zero matrix of num_bins*num_bins
heatmap_mat = zeros(num_bins, num_bins);

% loop though pairs and increment mat(bin1,bin2) += 1
for row=1:height(table)
	% row and column of the heatmap:
	i = int32(table_binned.bin_index1(row));
	j = int32(table_binned.bin_index2(row));
	% "put" a count in a corresponding "cell"
	% of a heatmap:
	heatmap_mat(i,j) = heatmap_mat(i,j) + 1;
	heatmap_mat(j,i) = heatmap_mat(j,i) + 1;
end

% display a heatmap of the matrix
m  = heatmap(heatmap_mat,'Colormap',flip(autumn),'ColorScaling','log','ColorLimits',[-1 5])
m.GridVisible = 0


% Draft for tomorrow

% scaling plot
% distance decay of the contact probability

[nrows, ncols] = size(heatmap_mat);

mean_list=[]
for i = 1:nrows
    % touch a bit upon what exactly "diag" does
    % depending on the input mat,vec, no index ...
   mean_list=[mean_list,mean(diag(heatmap_mat,i-1))];
end

% 1:nrows -> turn into basepairs *binsize
plot(1:nrows, mean_list )
plot(1:nrows, log(mean_list) )

% loglog

% Create the expected matrix by taking the means of the diahonals

% retire this
exp_matrix = zeros(ncols, ncols);
for i = 1:nrows
%    m=mean(diag(heatmap_mat,i-1)); 
   meansvals=logical(diag(heatmap_mat,i-1)).*mean(diag(heatmap_mat,i-1));
   upper=diag(meansvals,i-1);
   lower=diag(meansvals,1-i);
   exp_matrix = exp_matrix + upper;
   exp_matrix = exp_matrix + lower;
end

% simple way to create an expected matrix
% touch up on this toeplitz matrix concept
exp_matrix = toeplitz(mean_list);

h=heatmap(exp_matrix, 'Colormap',flip(autumn),'ColorScaling','log')
h.GridVisible = 'off'

% observed over expected ?
observed_over_expected=heatmap_mat./exp_matrix;
h=heatmap( observed_over_expected ,'Colormap',jet,'ColorScaling','log')
h.GridVisible = 'off'



% discuss Pearson correlation - linear regression
% and what corr does here ...
observed_over_expected(isnan(observed_over_expected)) = 0;
correlation_matrix_pre = 1+corr(observed_over_expected);
h = heatmap(correlation_matrix_pre,  'Colormap', jet,'ColorScaling','log');
% imagesc(correlation_matrix_pre);
% colorbar
h.GridVisible = 'off'
 
imshow(correlation_matrix_pre);
edge_image = edge(correlation_matrix_pre, 'canny');
imshow(edge_image);
 
edge_sums = sum(transpose(edge_image));
plot(edge_sums);

