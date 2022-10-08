function [laplacian] = LaplicanGeneration(K)%K=A为相似度矩阵
eye_matrix = 1 - eye(size(K));%对角线为0，非对角线为1
K = K .* eye_matrix;%对角线元素为0，其他元素不变
c_diag = sum(K, 2);%每行的和
c_diag(c_diag == 0) = 1;
c_diag(c_diag < 10^(-10)) = 10^(-10);
c_diag = diag(sqrt(c_diag.^(-1)));%c_diag为度矩阵D，D=D^(-1/2)
laplacian = eye(size(K)) - c_diag * K * c_diag;