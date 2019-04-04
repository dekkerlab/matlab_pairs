

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


%  2.
% explain what the data is:
% columns in the table - what they are?
% strands - understand them better !


% one little thing would be handy for us
% going forward - estimate chrom-size - how?
% we did't give you reference genome - so how'd
% you find chrom_len ?! chr19 human - google ?!
% chrom_len = max(table.pos1)
% chrom_len = max(table.pos2)
% chrom_len = max(max(table.pos2),max(table.pos1))


%  3.
% what can learn from this data ? ...
% let's just look at it !
scatter(table.pos2,table.pos1,0.8,'red','filled')
pbaspect([1 1 1])
xlabel("pos2, bp")
ylabel("pos1, bp")
title("Hi-C interactions for chrom 19")
%  uncomment next 2 lines in order to see symmetrical scatter:
% hold
% scatter(table.pos1,table.pos2,0.8,'blue','filled')

% what does this scatter-plot tells us ?

% can you guys, tall us what is the relationship between
% pos1 and pos2 ?

% where are the pairs/interactions that
% have the same genomic separation ?

% small genomic separations vs large genomic
% separations - any obvious differences ?

% uncomment 2 lines to see a diagonal:
% p = refline(slope, intercept)
% p.LineWidth = 5


% 4.
% what else was in the data-file - strands !
%  what do strands tell us ?

% how many combinations of strand orientation
% there would be ?

% let's count and see:
groupsummary(table,{'str1','str2'})
% one group is overrepresented
% why is that ?
% where do these ~120000 extra pairs are ?!

% 5.
% let's group the data by the strand orientation
% and explore it using scatterplots, again!

% groupby:
[G, Gkey1, Gkey2] = findgroups(table.str1,table.str2);
% what this is ?
G(1:15)
% just indices of those 4 combinations of strand orientations
% we can check it real quick (uncomment below if you want to check):
% hist(G)

%  num groups would be handy below
num_groups = max(G);

% groupped key-combinations:
horzcat(Gkey1,Gkey2)

% now, let's split table into 4:
table_pp = table(G==1,:);
table_pm = table(G==2,:);
table_mp = table(G==3,:);
table_mm = table(G==4,:);

% will have all 4 plots in one
subplot(2,2,1)
scatter(table_pp.pos1,table_pp.pos2,0.7,'red', 'filled')
subplot(2,2,2)
scatter(table_pm.pos1,table_pm.pos2,0.7,'blue', 'filled')
subplot(2,2,3)
scatter(table_mp.pos1,table_mp.pos2,0.7,'magenta', 'filled')
subplot(2,2,4)
scatter(table_mm.pos1,table_mm.pos2,0.7,'green', 'filled')

%  or we could do it in a for loop if we have time ...
for i = 1:num_groups
	subplot(2,2,i);
	scatter(table(G==i,:).pos2,table(G==i,:).pos1,0.7,'red', 'filled')
	title(join([Gkey1(i),Gkey2(i)],''))
end

% and even loop through colors ...
CM = jet(2*num_groups);
for i = 1:num_groups
	subplot(2,2,i);
	scatter(table(G==i,:).pos2,table(G==i,:).pos1,2.5,CM(i,:),'filled')
	title(join([Gkey1(i),Gkey2(i)],''))
	xlabel("pos2")
	ylabel("pos1")
	pbaspect([1 1 1])
end

% can you guys spot any differences ?!
% I certainly cannot!

% but where are this ~120'000 extra pairs anyways ?!

%  6.

%  query table and check for sizes ...

% pairs from PP that are more than 10'000 bp apart:
table_pp(table_pp.pos2-table_pp.pos1>10000,:)

% pairs from PM that are between 1000 and 10000 bp apart:
sum(abs(table_pm.pos2-table_pm.pos1)>1000 & abs(table_pm.pos2-table_pm.pos1)<10000)



% do some querying to find out where the difference is !
% maybe do a hist of a short range of pairs if we have time ...

% group, extract groups and do 4 "scalings"-counts-distributions ...
% there must be a difference in short-range interactions
% what does that mean ? - refer to the paper - and known HiC artifacts
% also mention that we do not know the enzyme cut sites, yet short pairs
% vs long pairs can be distinguished ...

% discussion , understanding all of that


%  7.
%  binning - homework ?!

% it ultimately leads to the question of
% interaction frequency - which pairs themselves do not
% represent - it's their "density" ...


% go to binning - decide if it's homework or not
% anyways introduce neccessary tools/functions
% show it briefly ...
% discretize function returns the indicies of the bin. 

% define bins ..., e.g., bins=0:binsize:chromlen
% use "discretize" with pre-set bins on pos1 pos2 ...
% there is bins functions - sounds intriguing what does it do?
% they woudl turn into 'bin1' 'bin2'
% then initialize empty/zero matrix of num_bins*num_bins
% loop though pairs and increment mat(bin1,bin2) += 1
chromlen = round(max(max(table.pos1),max(table.pos2)),-3)
binsize=1000;
bins=0:binsize:chromlen;

%show a pos2, pos1 scatter again and do a grid on, grid minor
% make sure grid is nicely visible
% discuss what would binning give us on top of the scatter ...
%we really want to estimate the contact frequency, which scatter plot
% nicely demonstrates, but it is not "that thing" yet ...

% concept of a contact frequency - how often does chr:1M-2M and chr:5M-6M
% in comparison with following 2 loci: chr:1M-2M and chr:2M-3M .

%todo - deal with scatter/grid etc
% draft a better binning code - and maybe keep some public, while other
% stuff "hidden" to make it a homework

discretize(table_pp(:,'pos1'),bins)


file1_binned = horzcat(transpose(bins()),transpose(bins(discretize(file1(:,2),bins))));
file2_binned = horzcat(transpose(bins(discretize(file2(:,1),bins))),transpose(bins(discretize(file2(:,2),bins))));
file3_binned = horzcat(transpose(bins(discretize(file3(:,1),bins))),transpose(bins(discretize(file3(:,2),bins))));
file4_binned = horzcat(transpose(bins(discretize(file4(:,1),bins))),transpose(bins(discretize(file4(:,2),bins))));

% these are the sep for individula files
file1_sep=file1_binned(:,2)-file1_binned(:,1);
file2_sep=file2_binned(:,2)-file2_binned(:,1);
file3_sep=file3_binned(:,2)-file3_binned(:,1);
file4_sep=file4_binned(:,2)-file4_binned(:,1);

% are we doinf this part or no?
hold on;
[n, xout] = hist(file1_sep);
plot(xout, n,'red');
hold on;
[n, xout] = hist(file2_sep);
plot(xout, n,'blue');
hold on;
[n, xout] = hist(file3_sep);
plot(xout, n,'magenta');
hold on;
[n, xout] = hist(file4_sep);
plot(xout, n,'green');

% % RETIRED CODE

% % let's just look at the overall distribution
% % of interactions with genomic separation
% % see how much it depends on the type of bins (unfortunately)
% % but it's fine
% % distance decay is still there - just depends on how you look
% %  at it

% % calculate genomic separation for the pairs
% % using nice named columns of the table:
% sep = table.pos2 - table.pos1;
% table = addvars(table,sep,'After','pos2');

% % let's plot overall histogramm of
% % genomic separations in linear bins
% [n, xout] = hist(table.sep);
% plot(xout, n);
% % doesn't look very informative
% loglog(xout, n)
% % a bit better ...
% % set(gca,'YScale','log');
% % set(gca,'XScale','log');

% %  how about log-spaced bins ...
% bins=logspace(1,7,20); % create bin edges with logarithmic scale
% [n, xout] = hist(table.sep,bins);
% plot(xout, n);
% % loglog(xout, n);

