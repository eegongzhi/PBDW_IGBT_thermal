function [Q, R] = gson(X,dS_field1, dS_field2)
    % Gram-Schmidt orthonormalization which produces the same result as [Q,R]=qr(X,0)
    
%     [d,n] = size(X);
%     m = min(d,n);
%     R = zeros(m,n);
%     Q = zeros(d,m);
%     for i = 1:m
%         R(1:i-1,i) = innerProduct(Q(:,1:i-1),X(:,i),dS_field1, dS_field2); % <Q(:,1:i-1), X(:,i)>
%         v = X(:,i)-Q(:,1:i-1)*R(1:i-1,i);
%         R(i,i) = sqrt(innerProduct(v,v,dS_field1, dS_field2));  % sqrt(<v, v>)
%         Q(:,i) = v/R(i,i);
%     end
%     R(:,m+1:n) = Q'*X(:,m+1:n);

    [m,n] = size(X);
    Q = zeros(m,n);
    R = zeros(n,n);
    for j=1:n
        v = X(:,j);
        for i = 1:j-1
%             R(i,j) = Q(:,i)' * X(:,j);
            R(i,j) = innerProduct(Q(:,i), X(:,j),dS_field1, dS_field2);
            v = v - R(i,j) * Q(:,i);
        end
%         R(j,j) = norm(v);
        R(j,j) = sqrt(innerProduct(v,v,dS_field1, dS_field2));
        Q(:,j) = v/R(j,j);
    end
end
