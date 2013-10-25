clear all,close all,clc;

%% Step 1: Load data
tmp = load('data.mat','data');
data = tmp.data;
n_toxicants = 20;
n_concentrations = 6;
n_genes = 39;
n_times = 25;

n_conditions = n_toxicants*n_concentrations;

data_use = zeros(n_conditions,n_genes,n_times);
for i = 1:(n_toxicants*n_concentrations)
    for j = 1:n_genes
       data_use(i,j,:) = data((i-1)*n_conditions+j);
    end
end

index_gene = 20;
figure 
hold on

for i = 1:n_conditions
    plot(data((i-1)*n_genes+index_gene,:))
end

% figure
% hold on
% index_gene = 1;
% 
% for i = 1:n_conditions
%     tmp = zeros(n_times,1);
%     for j = 1:n_times
%         tmp(j,1) = data_use(i,index_gene,j);
%     end
%     plot(tmp);
% end