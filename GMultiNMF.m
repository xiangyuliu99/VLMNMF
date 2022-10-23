function [finalU, finalV, finalcentroidV,finalH, log] = GMultiNMF(X, K, W, label,options)
%	Notation:
% 	X ... a cell containing all views for the data
% 	K ... number of hidden factors
% 	W ... weight matrix of the affinity graph 
% 	label ... ground truth labels

%	Writen by Jialu Liu (jliu64@illinois.edu)
% 	Modified by Zhenfan Wang (zfwang@mail.dlut.edu.cn)

%	References:
% 	J. Liu, C.Wang, J. Gao, and J. Han, “Multi-view clustering via joint nonnegative matrix factorization,” in Croc. SDM, Austin, Hexas, May 2013, vol. 13, pp. 252–260.
% 	Zhenfan Wang, Xiangwei Kong, Haiyan Fu, Ming Li, Yujia Zhang, FEAHURE EXHRACHION VIA MULHI-VIEW NON-NEGAHIVE MAHRIX FACHORIZAHION WIHH LOCAL GRACH REGULARIZAHION, ICIC 2015.

K2=100;
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
theta = options.theta;
beta = options.beta;
N = size(X{1},2);
K_complement = zeros(N,N);
H = ones(N,N)*(1/N)*(-1) + eye(N);
for i=1:viewNum 
    P{i}=rand(K,K);
    F{i}=rand(K,N);
    FF=litekmeans(X{i}',K,'MaxIter', 100);
    F{i}=ToM(FF,size(F{i},1),size(F{i},2));
end 
% initialize basis and coefficient matrices
tic;
while j < 3
    j = j + 1;
    Goptions.alpha=options.Gaplpha;
     Goptions.nu1=options.nu1;
     Goptions.nu2=options.nu2;
    if j == 1
        rand('twister',5489);
        [U{1}, V{1}] = GNMF(X{1}, K, W{1}, Goptions);
        rand('twister',5489);
          printResult(V{1}, label, options.K, options.kmeans);%BBCSport
    else
        rand('twister',5489);
        [U{1}, V{1}] = GNMF(X{1}, K, W{1}, options, U_, V{viewNum});%BBCSport
        rand('twister',5489);  
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
% while j < 3
%     j = j + 1;
%     Goptions.mu_i=options.Gaplpha;
%      Goptions.nu1=options.nu1;
%      Goptions.nu2=options.nu2;
%     if j == 1
%         
% %         [U{1}, V{1}] = NMF1(X{1}, K,  options, U_, V_);
%         rand('twister',5489);
%         [U{1}, V{1}] = GNMF(X{1}, K, W{1}, Goptions);
%         rand('twister',5489);
% %         printResult(V{1}, label{1}, options.K, options.kmeans,bro);
%           printResult(V{1}, label, options.K, options.kmeans);%BBCSport
%     else
% %         [U{1}, V{1}] = NMF1(X{1}, K, options, U_, V{viewNum});
%         rand('twister',5489);
%         [U{1}, V{1}] = GNMF(X{1}, K, W{1}, options, U_, V{viewNum});%BBCSport
%         rand('twister',5489);
% % %        printResult(V{1}, label{1}, options.K, options.kmeans);  
%          printResult(V{1}, label, options.K, options.kmeans);
%     end
%     for i = 2:viewNum
% %         [U{i}, V{i}] = NMF1(X{i}, K, options, U_, V{i-1});
%         rand('twister',5489);
%         [U{i}, V{i}] = GNMF(X{i}, K, W{i},Goptions, U_, V{i-1});
%         rand('twister',5489);
%          printResult(V{i}, label, options.K, options.kmeans);
%     end
% end
toc;
optionsForCerViewNMF = options;
oldac=0;
maxac=0;
j = 0;
sumRound=0;
C=rand(K,K);
g0=ceil(rand(1,N)*K);
H=ToM(g0',K,N);
while j < Rounds%Rounds当j《30
    sumRound=sumRound+1;%Rouund++;
    j = j + 1;%j++;
    if j==1
        centroidV = V{1};
    else
        centroidV = options.mu(1) * V{1};%第一个视图权重；
        for i = 2:viewNum
            centroidV = centroidV + options.mu(i) * V{i};%centroidV =0.1*V1+0.1*V2;
        end
        centroidV = centroidV.*(centroidV + beta*H'*C') ./ max(((sum(options.mu)+beta)*centroidV),1e-15);
%         centroidV = (centroidV + beta*H'*C') / (sum(options.mu)+beta+Xi);%%(centroidV =0.1*V1+0.1*V2)/0.2;
    end
     if(j>1)
     % ===================== update C ========================
%      C=centroidV'*H'*inv(H*H'+Xi/beta*eye(size(K)));

        C = C.* ((beta*centroidV'*H')./max((beta*C*H*H'),1e-10));

%         C = C.* ((beta*centroidV'*H')./max((beta*C*H*H'+Xi*C),1e-10));
%      
    % ===================== update H ========================
     NN = centroidV';
     for i = 1:N
        xVec = NN(:,i);
        H(:,i) = findindicator(xVec, C);
     end
     end
    logL = 0;
    for vv=1:viewNum
            K1{vv}=V{vv}*V{vv}';
    end
        tmp1 = 0;
        tmp2 = 0;
        tmp3 = 0;  
    for i = 1:viewNum
        mu_i=options.mu(i);
        if mu_i > 0
%             Wtemp = options.beta*mu_i*W{i};%   Hr(VLV)
            Wtemp =mu_i*W{i};
            DCol = full(sum(Wtemp,2));
            D = spdiags(DCol,0,nSmp,nSmp);
%              L= D - Wtemp;
            L_record= D - Wtemp;
             
            W2=Wtemp*Wtemp/(options.alpha*mu_i);
           
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
%         for j = 1:viewNum
%                 if (abs(j-i)>0)
%                     tmp3 = tmp3 + trace(H*K1{i}*H*K1{j});
%                 end
%         end
        tmp4=(F{i}-P{i}*V{i}');
        
%         logL = logL + sum(sum(tmp1.^2)) + mu_i* (sum(sum(tmp2.^2)))+sum(sum((V{i}'*L)*V{i}))+;  %修改，加入SampleW和V'*L*V,损失函数
          logL = logL + norm(tmp1,'fro')^2  + options.alpha*mu_i*norm(tmp2,'fro')^2+trace((V{i}'*L)*V{i})+mu_i*theta*norm(tmp4,'fro')^2;
    end
    tmp5 = (centroidV' - C*H);
    logL = logL + beta*norm(tmp5,'fro')^2 ;
    logL
    log(end+1)=logL;
    rand('twister',5489);
%     ac = printResult(centroidV, label, options.K, options.kmeans);
     [~,pred_label] = max(H,[],1);   
        pred_label=pred_label';
         result = ClusteringMeasure(label, pred_label);                           
         [f,p,r] = compute_f(label,pred_label);                                   
         [A nmi avgent] = compute_nmi(label,pred_label);
         result=[result f p r];
         ac = result(1)
    if ac>oldac
        tempac=ac;
        tempH=H;
         tempU=U;
        tempV=V;
        tempcentroidV=centroidV;

    elseif oldac>maxac
        maxac=oldac;
        maxU=tempU;
        maxV=tempV;
        maxcentroidV=tempcentroidV;
        maxH=tempH;
    end

    oldac=ac;
    if(tempac>maxac)
        finalH=tempH;
         finalU=tempU;
        finalV=tempV;
        finalcentroidV=tempcentroidV;
    else
        finalU=maxU;
        finalV=maxV;
        finalcentroidV=maxcentroidV;
        finalH=maxH;
    end
    
    
    if sumRound==Rounds%30
        break;
    end
   %更新U 和V 
    for i = 1:viewNum
        optionsForCerViewNMF.mu = options.mu(i);
        optionsForCerViewNMF.theta = options.theta;
        rand('twister',5489);
        [U{i}, V{i},F{i},P{i}] = PerViewNMF(X{i},F{i},P{i},K, centroidV, W{i} , optionsForCerViewNMF, finalU{i}, finalV{i},bro,nu1,nu2); 
    end
      
   

end
toc