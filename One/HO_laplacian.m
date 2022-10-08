function [L_record] = HO_laplacian(L_record)
L_record = diag(diag(L_record)) - L_record;%拉普拉斯矩阵变为对角元素为0，其他位置变为负数

% imagesc(L_record)
L_record = L_record * L_record;%L=L*L
L_record_P = L_record(L_record > 0);%L_record_P将L中大于0得数排列成一列
L_record_org = L_record;
mean_value = mean(L_record_P);%L_record_P中这一列的均值
std_value = std(L_record_P);%标准差
L_record(L_record < mean_value - std_value/2) = 0;%将L中小于L_record < mean_value - std_value/2改成0
[L_record] = LaplicanGeneration(L_record);%产生二阶拉普拉斯矩阵

L_record = kernel_completion(L_record, L_record_org);