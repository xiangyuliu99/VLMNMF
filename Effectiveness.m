%This is  a  sample demo
%Test Digits dataset
addpath('tools/');
addpath('print/');
options = [];
options.maxIter = 200;
options.error = 1e-6;
options.nRepeat = 30;
options.minIter = 50;
options.meanFitRatio = 0.1;
options.rounds = 60;
options.K=10;
% options.K=6;
options.Gaplpha=100;
options.WeightMode='Binary';
options.bro=6;%6
options.nu1=0.7;%0.1
options.nu2=0.3;%0.9


% options.kmeans means whether to run kmeans on v^* or not
% options alpha is an array of weights for different views
% options.alphas =0.01* ones(1,27);
options.alphas = [0.01 0.06];%0.01 0.06
% options.alphas = [0.01 0.03 0.01 0.01 0.01 0.01 0.04 0.01 0.01 0.01 0.01 0.06 0.03];%λ=0.01,0.05
options.kmeans = 1;%代表需不需要在V上需不需要进行kmeans
options.beta=11;%μ=11




%% read dataset
%   load 3sources.mat;
load handwritten.mat;
% load BBCSport.mat
%  load proteinFold_Kmatrix.mat;
%     load 20newsgroups.mat;
data=cell(1,2);
 data{1} = fourier';
 data{2} = pixel';   %2个视角
% K = 6;
X=data;
K = options.K;
N = size(data{1},2);
K_complement = zeros(N,N);
H = ones(N,N)*(1/N)*(-1) + eye(N);

%% normalize data matrix
for i = 1:length(data)
%     dtemp=computeDistMat(data{i},2);
%     W{i}=constructW(dtemp,20);
%     data{i} = data{i} / sum(sum(data{i}));
    options.WeightMode='Binary';
    W{i}=constructW_cai(data{i}',options);
    data{i} = data{i} / sum(sum(data{i}));
end


%%
% [W] = V9_LocalKernelCalculation(KH , 1, K);%KHN为相似度矩阵A
% W1=cell(1,size(W,3));
% for i=1:size(W,3)
%     for j=1:size(W,3)
%         if i==j
%             W1{i}=W(:,:,j);
%         end
%     end
% end
% W=W1;
%%
% for i = 1:length(data)
% %     dtemp=computeDistMat(data{i},2);
% %     W{i}=constructW(dtemp,20);
% %     data{i} = data{i} / sum(sum(data{i}));
%     data{i} = abs(data{i}) / sum(sum (abs(data{i})));
%     
% end

% run 20 times
% U_final = cell(1,3);
% V_final = cell(1,3);
% V_centroid = cell(1,3);
U_final = cell(1,3);
V_final = cell(1,3);
V_centroid = cell(1,3);
gnd=(gnd+1);
% truelabel=truelabel{1}';
for i = 1:4
    if i==1
        options.beta=0;
    end
    if i==2
        options.bro=0;
    end
    if i==3
        options.nu1=0;%0.1
    end
    if i==4
        options.nu2=0;%0.9
    end
    
    [U_final{i}, V_final{i}, V_centroid{i} log] = GMultiNMF(data, K, W,gnd, options);
%       [U_final{i}, V_final{i}, V_centroid{i} log] = GMultiNMF(data, K, W,Y, options);
%       printResult( V_centroid{i}, Y, K, options.kmeans);
       a(i)=printResult( V_centroid{i}, gnd, K, options.kmeans);
   fprintf('\n');
    options.bro=6;%10
   options.beta=11;
   options.nu1=0.7;%0.1
   options.nu2=0.3;%0.9
end
