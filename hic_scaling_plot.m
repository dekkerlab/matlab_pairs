filename = 'U54_HFF_plate_subset.txt';
table = readtable(filename);
head(table)


% way to etract columns from table:
% table.Var3 - 3rd column

% genomic separation here:
sep = table.Var3 - table.Var2;
% aiming to do :
hist(sep);

% log-spaced histogram matlab ...
y=100*rand(1000,1); % create random data
x=logspace(-1,2,20); % create bin edges with logarithmic scale
histogram(y,x); % create the plot
set(gca,'xscale','log'); % scale the x-axis


% also demonstrate 
scatter(table.Var3,table.Var2,0.5,'red')


% figure out grouping better - something like that:
sep1 = table.Var3(G==2) - table.Var2(G==2);
sep1 = table.Var3(G==1) - table.Var2(G==1);
sep2 = table.Var3(G==2) - table.Var2(G==2);
sep3 = table.Var3(G==3) - table.Var2(G==3);

strands=table(:,4:5);
% strand2=table{5};
G=findgroups(strands);
len=length(G)
% k=splitapply(table,k)
for i=1:len
    if G(i)=='1' 
        A_=table(i,:);
        A=vertcat(A,A_);
    end
    if G(i)=='2' 
        B_=table(i,:);
        B=vertcat(B,B_);
    end
    if G(i)=='3' 
        C_=table(i,:);
        C=vertcat(C,C_);
    end
    if G(i)=='4' 
        D_=table(i,:);
        D=vertcat(D,D_); 
    end
end
