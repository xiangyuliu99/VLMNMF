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
options.bro=0;%0.1
options.nu1=0.66;%0.66
options.nu2=0.34;%0.34

options.mu = [0.05 0.09];%0.05 0.09   
options.kmeans = 1;%代表需不需要在V上需不需要进行kmeans
options.alpha=11;%μ=11
options.theta =0.1;%0.1 
options.eta = 0;%1  
options.beta = 0.001;%0.001   0.9715


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
%      data{i}=abs(mapminmax(data{i},-1,1));
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
iter = 1;
options.alpha=0.001;%μ=0.01
for i = 1:7
    [U_final{1}, V_final{1}, V_centroid{1},T, log] = GMultiNMF(data, K, W,gnd, options);
        [~,pred_label] = max(T,[],1);   
        pred_label=pred_label';
%          result = ClusteringMeasure(gnd, pred_label);                           
%          [f,p,r] = compute_f(gnd,pred_label);                                   
%          [A nmi avgent] = compute_nmi(gnd,pred_label);
         
         [ac, nmi_value, cnt,AR,F,P,R] = CalcMetrics(gnd, pred_label);
         result=[ac nmi_value AR F P R];
disp(sprintf('ac: %0.4f\t%d/%d\tnmi:%0.4f\t  AR:  %0.4f\t fscore: %0.4f\t  P: %0.4f\t  R:  %0.4f\t', ac, cnt, length(gnd), nmi_value,AR,F,P,R));
        
        HWacc(iter) = result(1); 
        HWnmi(iter) = result(2);
        HWAR(iter) = result(3);
        HWF(iter) = result(4);
         HWP(iter) = result(5);
        HWR(iter) = result(6);
        iter =iter+1;
        options.alpha=options.alpha*10;
        fprintf('\n');
end
