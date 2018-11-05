function [y,grady] = quadobj(x2,x1,Q,f,c)
[V, D] = eig(Q);
eigvalue=diag(D);
if all(eigvalue>=0)
    y = 1/2*x2'*Q*x2 + f'*x2 + c;
    if nargout > 1
        grady = Q*x2 + f;
    end
else
    Qn=zeros(size(Q));
    % init index
    if any(eigvalue>=0)
        index=find(eigvalue>=0); 
    else
        index=length(Q);
    end
    for i=1:index
        Qn=Qn-eigvalue(i)*V(:,i)*V(:,i)';
    end
    Qp=zeros(size(Q));
    for i=index:length(Q)
        Qp=Qp+eigvalue(i)*V(:,i)*V(:,i)';
    end
    temp=1/2*x1'*Qn*x1;
    y=1/2*x2'*Qp*x2+f'*x2+c - (temp+(Qn*x1)'*(x2-x1));
    if nargout > 1
        grady=Qp*x2+f-(Qn*x1);
    end
      
end
end