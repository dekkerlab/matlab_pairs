
%  1.
% read some pair data in the form of a table
filename = 'U54_HFF_plate_subset.txt';
table = readtable(filename);
% give those columns human readable names and types
table.Properties.VariableNames = {'chrom','pos1','pos2','str1','str2'};
table.Properties.VariableUnits = {'string', 'int32', 'int32', 'char', 'char'};
% check to see how it looks like
head(table)

% what is a chrom length
% from the data ?
% chromlen = max(max(pos1),max(pos2)) - best you could do from this data 

%  2.
% explain what the data is:
% columns in the table - what they are?
% strands - understand them better !

%  3.
% after we show the table with pais
% let's go straight ahead to the scatter
% of pos2 pos1, which would demonstrate
% distance decay, i.e. more interactions
% between nearby loci (small genomic separation)

% how does this even look like ... 
scatter(table.pos2,table.pos1,0.5,'red')
% also it would ultimately lead to the question of
% interaction frequency - which pairs themselves do not
% represent - it's their "density" ...

%  4.
% let's just look at the overall distribution
% of interactions with genomic separation
% see how much it depends on the type of bins (unfortunately)
% but it's fine
% distance decay is still there - just depends on how you look
%  at it

% calculate genomic separation for the pairs
% using nice named columns of the table:
sep = table.pos2 - table.pos1;
table = addvars(table,sep,'After','pos2');

% let's plot overall histogramm of
% genomic separations in linear bins
[n, xout] = hist(table.sep);
plot(xout, n);
% doesn't look very informative
loglog(xout, n)
% a bit better ...
% set(gca,'YScale','log');
% set(gca,'XScale','log');

%  how about log-spaced bins ...
bins=logspace(1,7,20); % create bin edges with logarithmic scale
[n, xout] = hist(table.sep,bins);
plot(xout, n);
% loglog(xout, n);

%  5.
% what else is in the data ?
% we could group pairs by strands - and
% immediately interesting thing pops up:
groupsummary(table,{'str1','str2'})
% one group is overrepresented
% we could group table by strands and plot
% 4 different scatters ... - should not be much difference ...

G=findgroups(table.str1,table.str2});
% do something like that to extract groups ...

file1 = horzcat(table.pos1(G==1),table.pos2(G==1))
file2 = horzcat(table.pos1(G==2),table.pos2(G==2))
file3 = horzcat(table.pos1(G==2),table.pos2(G==2))
file4 = horzcat(table.pos1(G==2),table.pos2(G==2))

% will have all 4 plots in one
subplot(2,2,1)
scatter(table.pos1(G==1),table.pos2(G==1),0.5,'red')
subplot(2,2,2)
scatter(table.pos1(G==2),table.pos2(G==2),0.5,'blue')
subplot(2,2,3)
scatter(table.pos1(G==3),table.pos2(G==3),0.5,'magenta')
subplot(2,2,4)
scatter(table.pos1(G==4),table.pos2(G==4),0.5,'green')

% i think G -> 1,2,3,4 shoudl correspond to whatever was in
% groupsummary - but we need to check that.

%  you learn about groupsummary from G as well ...


% what about interaction distribution with genomic separation
% can we narrow down the differences between 4 groups to a certain
% separation-scale ?

% group, extract groups and do 4 "scalings"-counts-distributions ...
% there must be a difference in short-range interactions
% what does that mean ? - refer to the paper - and known HiC artifacts
% also mention that we do not know the enzyme cut sites, yet short pairs
% vs long pairs can be distinguished ...

% discussion , understanding all of that

% go to binning - decide if it's homework or not
% anyways introduce neccessary tools/functions
% show it briefly ...
% discretize function returns the indicies of the bin. 
chromlen =round(max(max(table.pos1),max(table.pos2)),-3)
binsize=1000;
bins=0:binsize:chromlen;

file1_binned = horzcat(transpose(bins(discretize(file1(:,1),bins))),transpose(bins(discretize(file1(:,2),bins))));
file2_binned = horzcat(transpose(bins(discretize(file2(:,1),bins))),transpose(bins(discretize(file2(:,2),bins))));
file3_binned = horzcat(transpose(bins(discretize(file3(:,1),bins))),transpose(bins(discretize(file3(:,2),bins))));
file4_binned = horzcat(transpose(bins(discretize(file4(:,1),bins))),transpose(bins(discretize(file4(:,2),bins))));

% these are the sep for individula files
file1_sep=file1_bin(:,2)-file1_bin(:,1);
file2_sep=file2_bin(:,2)-file2_bin(:,1);
file3_sep=file3_bin(:,2)-file3_bin(:,1);
file4_sep=file4_bin(:,2)-file4_bin(:,1);

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

% define bins ..., e.g., bins=0:binsize:chromlen
% use "discretize" with pre-set bins on pos1 pos2 ...
% there is bins functions - sounds intriguing what does it do?
% they woudl turn into 'bin1' 'bin2'
% then initialize empty/zero matrix of num_bins*num_bins
% loop though pairs and increment mat(bin1,bin2) += 1

% unedited - to be reused:


% let's start grouping ...
% strands=table(:,4:5);
% % strand2=table{5};
% G=findgroups(strands);
% len=length(G)
% % k=splitapply(table,k)







% 
% 
% % figure out grouping better - something like that:
% sep1 = table.Var3(G==2) - table.Var2(G==2);
% sep1 = table.Var3(G==1) - table.Var2(G==1);
% sep2 = table.Var3(G==2) - table.Var2(G==2);
% sep3 = table.Var3(G==3) - table.Var2(G==3);
% 
% strands=table(:,4:5);
% % strand2=table{5};
% G=findgroups(strands);
% len=length(G)
% % k=splitapply(table,k)
% for i=1:len
%     if G(i)=='1' 
%         A_=table(i,:);
%         A=vertcat(A,A_);
%     end
%     if G(i)=='2' 
%         B_=table(i,:);
%         B=vertcat(B,B_);
%     end
%     if G(i)=table='3' 
%         C_=table(i,:);
%         C=vertcat(C,C_);
%     end
%     if G(i)=='4' 
%         D_=table(i,:);
%         D=vertcat(D,D_); 
%     end
% end
