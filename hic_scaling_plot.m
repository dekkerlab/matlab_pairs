% read some pair data in the form of a table
filename = 'U54_HFF_plate_subset.txt';
table = readtable(filename);
% give those columns human readable names and types
table.Properties.VariableNames = {'chrom','pos1','pos2','str1','str2'};
table.Properties.VariableUnits = {'string', 'int32', 'int32', 'char', 'char'};
% check to see how it looks like
head(table)

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

% how does this even look like ... 
scatter(table.pos1,table.pos2,0.5,'red')

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
