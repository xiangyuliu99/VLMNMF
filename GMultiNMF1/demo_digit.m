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
options.K=15;
% options.K=6;
options.Gaplpha=100;%100
options.WeightMode='Binary';
options.bro=100;%10

options.nu1=0.9;%0.1
options.nu2=0.1;%0.9



% options.kmeans means whether to run kmeans on v^* or not
% options alpha is an array of weights for different views
% options.alphas =0.01* ones(1,27);
options.alphas = [1 1 8];%[0.01 0.01 0.08];
% options.alphas = [0.01 0.03 0.01 0.01 0.01 0.01 0.04 0.01 0.01 0.01 0.01 0.06 0.03];%λ=0.01,0.05
options.kmeans = 1;%代表需不需要在V上需不需要进行kmeans
options.beta=0.5;%μ=0.5   69.09   73.38




%% read dataset
%   load 3sources.mat;
% load handwritten.mat;
% load BBCSport.mat
 load yale_mtv.mat;
%     load 20newsgroups.mat;
%  data{1} = fourier';
%  data{2} = pixel';   %2个视角
% K = 6;
K = options.K;
data=X;
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
    data{i} = abs(data{i}) / sum(sum (abs(data{i})));
    
end

%save('handwrittenW','W');



% run 20 times
% U_final = cell(1,3);
% V_final = cell(1,3);
% V_centroid = cell(1,3);
U_final = cell(1,4);
V_final = cell(1,4);
V_centroid = cell(1,4);
% truelabel=truelabel{1}';
for i = 1:1
   
%     [U_final{i}, V_final{i}, V_centroid{i} log] = GMultiNMF(data, K, W,gnd, options);
      [U_final{i}, V_final{i}, V_centroid{i} log] = GMultiNMF(data, K, W,gt, options);
      a(i)=printResult( V_centroid{i}, gt, K, options.kmeans);
%         printResult( V_centroid{i}, gnd, K, options.kmeans);
   fprintf('\n');
end
