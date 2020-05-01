networks = dir('./problem_set_2_sample/network/*.dat');
covariates = dir('./problem_set_2_sample/group/*.dat');

pre1 = './problem_set_2_sample/network/';
pre2 = './problem_set_2_sample/group/';

W = {};
%%
for i= 1:length(networks)
    
    filepath = fullfile(pre1, networks(i).name);
    W{i} = importdata(filepath);
    
end

Y = {}; % GPA
X = {}; % Covariates
G = {}; % mid-variable to read all columns

for i= 1:length(covariates)
    
    filepath = fullfile(pre2, covariates(i).name);
    G{i} = importdata(filepath);
    
    
    Y{i} = G{i}(:, 18);
    X{i} = G{i}(:, [1:6, 8:14]);
    
end
%%
save data W X Y;






