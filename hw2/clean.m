% cd './problem_set_2_sample/group'
filePattern = fullfile('*.dat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(baseFileName);
  %fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  % such as reading it in, plotting it, using it, or whatever.
  A{k} = readtable(fullFileName);
end

%%

cd '../network'
filePattern = fullfile('*.dat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(baseFileName);
  %fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  % such as reading it in, plotting it, using it, or whatever.
  B{k} = importdata(fullFileName);
end

%%

% distinguish GPA from the group data
% Convert the data from table to double
for i=1:76
    A{1,i} =removevars(A{1,i},'Var7');
    A{1,i} =removevars(A{1,i},'Var15');
    A{1,i} =removevars(A{1,i},'Var16');
    A{1,i} =removevars(A{1,i},'Var17');
    A{1,i} =removevars(A{1,i},'Var19');
end

for k=1:76
    g = strcat('A',num2str(k));
    C.(g) = A{1,k}{:,:} ;
end

for k=1:76
    g = strcat('B',num2str(k));
    h = strcat('A',num2str(k));
    D.(g) = B{1,k}{:,:} ;
    F.(g) = C.(h)(:,14);
end


for i = 1:76
    a = ['A',num2str(i)];
    b = ['B',num2str(i)];
    c = ['c',num2str(i)];
    E.(c) = D.(b) * C.(a);
end








