clear
clc
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
options.Gaplpha=100;
options.WeightMode='Binary';
options.bro=0;%10
options.nu1=0.61;%0.61
options.nu2=0.39;%0.39

options.mu =  [0.01 0.01 0.08];%0.01 0.01 0.08   
options.kmeans = 1;%代表需不需要在V上需不需要进行kmeans
options.alpha=0.01;%μ=0.01
options.theta =0.001;%0.001 
options.beta = 0.001;%0.001   0.7091


%% read dataset

 load yale_mtv.mat;

K = options.K;
data=X;
N = size(data{1},2);
K_complement = zeros(N,N);
H = ones(N,N)*(1/N)*(-1) + eye(N);
%% normalize data matrix
for i = 1:length(data)
    options.WeightMode='Binary';
    W{i}=constructW_cai(data{i}',options);
    data{i} = abs(data{i}) / sum(sum (abs(data{i})));
end
U_final = cell(1,3);
V_final = cell(1,3);
V_centroid = cell(1,3);
gnd = im2double(gt);
for i = 1:165
    gnd(i) = gt(i);
end
    [U_final{1}, V_final{1}, V_centroid{1},T, log] = GMultiNMF(data, K, W,gnd, options);
        [~,pred_label] = max(T,[],1);   
        pred_label=pred_label';
         [ac, nmi_value, cnt,AR,F,P,R] = CalcMetrics(gnd, pred_label);
         result=[ac nmi_value AR F P R];
disp(sprintf('ac: %0.4f\t%d/%d\tnmi:%0.4f\t  AR:  %0.4f\t fscore: %0.4f\t  P: %0.4f\t  R:  %0.4f\t', ac, cnt, length(gnd), nmi_value,AR,F,P,R));
        Yaleacc = result(1); 
        Yalenmi = result(2);
        YaleAR = result(3);

