function [ac] = printResult(X, label, K, kmeansFlag)

if kmeansFlag == 1
    indic = litekmeans(X, K, 'Replicates',20);
else
    [~, indic] = max(X, [] ,2);
end
%  result = bestMap(label, indic);
 if size(label)~=size(indic)
    label=label';
end
[ac, nmi_value, cnt,AR,F,P,R] = CalcMetrics(label, indic);
disp(sprintf('ac: %0.4f\t%d/%d\tnmi:%0.4f\t  AR:  %0.4f\t fscore: %0.4f\t  P: %0.4f\t  R:  %0.4f\t', ac, cnt, length(label), nmi_value,AR,F,P,R));