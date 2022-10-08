iter = 1;
A = 0.001
for i = 1:7
    B = 0.001;
    for j = 1:7
    
       acc(iter) = iter;
       iter =iter+1;
        B = B*10;
        fprintf('\n');
    end
    A = A*10;
end
