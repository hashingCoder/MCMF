function  [b] = DGHMF(net,init_b,gamma,epsilon,max_iter)
b =init_b;
[n,d] = size(init_b);
S = net;

loss=zeros(max_iter+1,1);

for iter=1:max_iter
    
    q=proj_stiefel_manifold( S' * b);
    I_1 = ones(n,n);
    J = eye(n)-1/n*I_1;
    
    Y = proj_stiefel_manifold( J * b);
    
    
    xx = S*q +  gamma*Y;
    
    idx0 =  find(xx==0);
    xx(idx0) = b(idx0);
    
    b = sign(xx);
    
    loss(iter+1)=norm(S-b*q','fro')+gamma*norm(b-Y,'fro');
    
    if abs(loss(iter+1) - loss(iter)) < epsilon
        break
    end
end


end


