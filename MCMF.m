function  [b] = MCMF(net,init_b,gamma,epsilon,max_iter)
    b =init_b;
    [n,d] = size(init_b);
    Y = b;
    %%make sure n is even
    lo=0;
    if mod(n,2)~=0
        lo = 2 - mod(n,2);
        net = [net;sparse(lo,n-lo)];
        net = [net sparse(n,lo)];
    end

    S = net;
    
    loss=zeros(max_iter+1,1);
    
    for iter=1:max_iter
        
        q=proj_stiefel_manifold( S' * b);
        
        if gamma~=0
            Y = proj_stiefel_manifold( S' * q)*sqrt(n);
            c = - vec((S*q)) - gamma * vec(Y);
        else
            c = - vec((S*q)) ;
        end
        
        
        b_tem=ones(n*d,1);
        
        idx_tem =[];
        
        for i = 1:d
            c_tem = c((i-1)*n +1 : i*n);
            [~,idx]=maxk(c_tem,n/2);
            idx_tem = [ idx_tem ; (i-1)*n + idx];
        end
        
        neg_one = - ones(length(idx_tem) ,1);
        b_tem(idx_tem) = neg_one;
        b = reshape(b_tem,n,d);
        
        
        loss(iter+1)=norm(S-b*q','fro')+gamma*norm(b-Y,'fro');
        
        if abs(loss(iter+1) - loss(iter)) < epsilon
            break
        end
    end
    
    b = b(1:n-lo,:);
    
end


