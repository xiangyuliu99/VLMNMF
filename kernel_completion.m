function kernel = kernel_completion(kernel, org_kernel)

Avg_index = kernel ~= 0;
Ker_sum = sum(Avg_index, 2);
index = find(Ker_sum == 0);
threshold = mean(kernel(:));
if ~isempty(index)
    org_kernel = org_kernel - diag(diag(org_kernel));
    Small_samples = org_kernel(index, :);
    [~, smp_indexes] = sort(Small_samples,2,'descend');
    smp_indexes = smp_indexes(:,1);
    
    for ii = 1 : size(Small_samples, 1)
        kernel(index(ii), smp_indexes(ii)) = threshold;
    end
end