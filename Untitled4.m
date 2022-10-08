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
options.K=100;
% options.K=6;
options.Gaplpha=100;
options.WeightMode='Binary';
options.bro=100;%0
options.nu1=0.5;%0.5
options.nu2=0.5;%0.5

options.alphas = [0.02 0.01 0.01];%0.02 0.01 0.01
options.kmeans = 1;%代表需不需要在V上需不需要进行kmeans
options.beta=1.2;%μ=1.2.
options.theta =100;%100
options.eta = 0.001;%0.001  0.7769   
options.Nu = 0.001;%0.001   0.9715
options.Xi = 0;

%% read dataset
%   load 3sources.mat;
load 100leaves.mat;

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
  data{i} = abs(data{i}) / sum(sum (abs(data{i})));
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
gnd=truelabel{1};
% truelabel=truelabel{1}';
iter = 1;    %0.7819
for i = 1:7
    options.theta = 0.001;
%     for j = 1:7
    [U_final{1}, V_final{1}, V_centroid{1},T,log] = GMultiNMF(data, K, W,gnd, options);
        [~,pred_label] = max(T,[],1);   
        pred_label=pred_label';
         result = ClusteringMeasure(gnd, pred_label);                           
         [f,p,r] = compute_f(gnd,pred_label);                                   
         [A nmi avgent] = compute_nmi(gnd,pred_label);
         result=[result f p r];
         fprintf('\nAll view results: ACC = %.4f, NMI = %.4f, Purity = %.4f, F-score = %.4f, Precision = %.4f and , Recall = %.4f\n',result(1),result(2),result(3),result(4),result(5),result(6));
        ac(iter) = result(1); 

       iter =iter+1;
        
        fprintf('\n');
        options.Nu = options.Nu*10;
end
