%  1.
% read pairs from a text file as a table.
% Explain table vs Matrix along the way:
%  - homogeneous vs hetergeneous
%  - every matrix is a "table", not every "table" is a matrix
filename = 'U54_HFF_plate_subset.txt';
table = readtable(filename);
% give those columns human readable names and types
table.Properties.VariableNames = {'chrom','pos1','pos2','str1','str2'};
table.Properties.VariableUnits = {'string', 'int32', 'int32', 'char', 'char'};
% check to see how it looks like
head(table)
% show how to access a given column:
table.pos1(1:10)

% 2.
% define bins ..., e.g., bins=0:binsize:chromlen
% use "discretize" with pre-set bins on pos1 pos2 ...
% there is bins functions - sounds intriguing what does it do?
% they woudl turn into 'bin1' 'bin2'

% define the bins
chromlen = round(max(max(table.pos1),max(table.pos2)),-3);
binsize=10000; % 10kb bins -- feel free to try different size bins
bin_edges=0:binsize:chromlen;

% Let bin pos1 into "bins" i.e.
% assign each position to corresponding bin range
% e.g: pos1=25000 will belong to bin_index1=3 and
% corresponding to bin_ranges1=(20000-30000)

bin_index1 = discretize(table.pos1,bins);
binned_pos1 = transpose(bin_index1); % turn row into column

% repeat the same for the second position

% append column of bin index back to our initial table

table_binned = addvars(table,binned_pos1,'After','pos1');
table_binned = addvars(table_binned,binned_pos2,'After','pos2');

% 3.

num_bins=length(bin_edges)-1; % compare with max(bin_index1)

% then initialize empty/zero matrix of num_bins*num_bins

% loop though pairs and increment mat(bin1,bin2) += 1

% display a heatmap of the matrix

% you've got your heatmap!
% it represents frequency of interactions between
% pairs of loci ...
