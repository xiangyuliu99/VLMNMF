function [finalU, finalV, finalcentroidV, log] = GMultiNMF(X, K, W, label,options)
%	Notation:
% 	X ... a cell containing all views for the data
% 	K ... number of hidden factors
% 	W ... weight matrix of the affinity graph 
% 	label ... ground truth labels

%	Writen by Jialu Liu (jliu64@illinois.edu)
% 	Modified by Zhenfan Wang (zfwang@mail.dlut.edu.cn)

%	References:
% 	J. Liu, C.Wang, J. Gao, and J. Han, “Multi-view clustering via joint nonnegative matrix factorization,” in Proc. SDM, Austin, Texas, May 2013, vol. 13, pp. 252–260.
% 	Zhenfan Wang, Xiangwei Kong, Haiyan Fu, Ming Li, Yujia Zhang, FEATURE EXTRACTION VIA MULTI-VIEW NON-NEGATIVE MATRIX FACTORIZATION WITH LOCAL GRAPH REGULARIZATION, ICIP 2015.


tic;
viewNum = length(X);
Rounds = options.rounds;
nSmp=size(X{1},2);
U_ = [];
V_ = [];

U = cell(1, viewNum);
V = cell(1, viewNum);

j = 0;
log = 0;
ac=0;
bro=options.bro;
nu1=options.nu1;
nu2=options.nu2;
N = size(X{1},2);
K_complement = zeros(N,N);
H = ones(N,N)*(1/N)*(-1) + eye(N);
% initialize basis and coefficient matrices
tic;
while j < 3
    j = j + 1;
    Goptions.alpha=options.Gaplpha;
     Goptions.nu1=options.nu1;
     Goptions.nu2=options.nu2;
    if j == 1
        
%         [U{1}, V{1}] = NMF1(X{1}, K,  options, U_, V_);
        rand('twister',5489);
        [U{1}, V{1}] = GNMF(X{1}, K, W{1}, Goptions);
        rand('twister',5489);
%         printResult(V{1}, label{1}, options.K, options.kmeans,bro);
          printResult(V{1}, label, options.K, options.kmeans);%BBCSport
    else
%         [U{1}, V{1}] = NMF1(X{1}, K, options, U_, V{viewNum});
        rand('twister',5489);
        [U{1}, V{1}] = GNMF(X{1}, K, W{1}, options, U_, V{viewNum});%BBCSport
        rand('twister',5489);
% %        printResult(V{1}, label{1}, options.K, options.kmeans);  
         printResult(V{1}, label, options.K, options.kmeans);
    end
    for i = 2:viewNum
%         [U{i}, V{i}] = NMF1(X{i}, K, options, U_, V{i-1});
        rand('twister',5489);
        [U{i}, V{i}] = GNMF(X{i}, K, W{i},Goptions, U_, V{i-1});
        rand('twister',5489);
         printResult(V{i}, label, options.K, options.kmeans);
    end
end
toc;
optionsForPerViewNMF = options;
oldac=0;
maxac=0;
j = 0;
sumRound=0;
while j < Rounds%Rounds当j《30
    sumRound=sumRound+1;%Rouund++;
    j = j + 1;%j++;
    if j==1
        centroidV = V{1};
    else
        centroidV = options.alphas(1) * V{1};%第一个视图权重；
        for i = 2:viewNum
            centroidV = centroidV + options.alphas(i) * V{i};%centroidV =0.1*V1+0.1*V2;
        end
        centroidV = centroidV / sum(options.alphas);%%(centroidV =0.1*V1+0.1*V2)/0.2;
    end
    logL = 0;
    for vv=1:viewNum
            K1{vv}=V{vv}*V{vv}';
    end
        tmp1 = 0;
        tmp2 = 0;
        tmp3 = 0;  
    for i = 1:viewNum
        alpha=options.alphas(i);
        if alpha > 0
            Wtemp = options.beta*alpha*W{i};%   Tr(VLV)
            DCol = full(sum(Wtemp,2));
            D = spdiags(DCol,0,nSmp,nSmp);
%              L= D - Wtemp;
            L_record= D - Wtemp;
             
            W2=Wtemp*Wtemp/(options.beta*alpha);
           
            DD = full(sum(W2,2));
            D2 = spdiags(DD,0,nSmp,nSmp);
           
            Wtemp=(nu1*Wtemp+nu2*W2);
            D=(nu1*D+nu2*D2);
            highorder_L = cell(1, 2);
            highorder_L{1} = L_record;
            highorder_L{2}=HO_laplacian(L_record);
             
            L=(nu1*highorder_L{1}+ nu2*highorder_L{2});%
            if isfield(options,'NormW') && options.NormW
                D_mhalf = spdiags(DCol.^-.5,0,nSmp,nSmp) ;
%                  L = D_mhalf*L*D_mhalf;
                L_record = D_mhalf*L*D_mhalf;
                highorder_L = cell(1, 2);
                highorder_L{1} = L_record;
                highorder_L{2}=HO_laplacian(L_record);
                
                L=(0.5*highorder_L{1}+ 0.5*highorder_L{2});
               
            end
        else
            L = [];
        end
        tmp1 = (X{i} - U{i}*V{i}');
        tmp2 = (V{i} - centroidV);
        for j = 1:viewNum
                if (abs(j-i)>0)
                    tmp3 = tmp3 + trace(H*K1{i}*H*K1{j});
                end
        end
%         logL = logL + sum(sum(tmp1.^2)) + alpha* (sum(sum(tmp2.^2)))+sum(sum((V{i}'*L)*V{i}))+;  %修改，加入SampleW和V'*L*V,损失函数
          logL = logL + norm(tmp1,'fro')^2  + alpha* norm(tmp2,'fro')^2+alpha*trace((V{i}'*L)*V{i})+bro*tmp3;
    end
    
    logL
    log(end+1)=logL;
    rand('twister',5489);
    ac = printResult(centroidV, label, options.K, options.kmeans);
    if ac>oldac
        tempac=ac;
        tempU=U;
        tempV=V;
        tempcentroidV=centroidV;

    elseif oldac>maxac
        maxac=oldac;
        maxU=tempU;
        maxV=tempV;
        maxcentroidV=tempcentroidV;
    end

    oldac=ac;
    if(tempac>maxac)
        finalU=tempU;
        finalV=tempV;
        finalcentroidV=tempcentroidV;

    else
        finalU=maxU;
        finalV=maxV;
        finalcentroidV=maxcentroidV;
    end
    
    
    if sumRound==Rounds%30
        break;
    end
   %更新U 和V 
    for i = 1:viewNum
        optionsForPerViewNMF.alpha = options.alphas(i);
        rand('twister',5489);
         K_complement= K_complement*0;
        for k=1:viewNum
            if (abs(k-i)>0) 
            K_complement =  K_complement + H*V{k}*V{k}'*H;
            end
        end
        [U{i}, V{i}] = PerViewNMF(X{i}, K, centroidV, W{i} , optionsForPerViewNMF, finalU{i}, finalV{i},K_complement,bro,nu1,nu2); 
    end

end
toc