function [AC, nmi_value, error_cnt,AR,F,P,R] = CalcMetrics(label, result)

result = bestMap(label, result);
error_cnt = sum(label ~= result);
AC = length(find(label == result))/length(label);
% F=fscore(label, K, result);
[AR,RI,MI,HI]=RandIndex(label,result);
[F,P,R] = compute_f(label,result);
nmi_value = nmi(label, result);

