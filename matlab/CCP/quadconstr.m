function [y,yeq,grady,gradyeq] = quadconstr(x2,x1,H,k,d)
jj = length(H); % jj is the number of inequality constraints
y = zeros(1,jj);
grady = zeros(length(x2),jj);
for i = 1:jj
    [V, D] = eig(H{i});
    eigvalue=diag(D);
    eigvalue=real(eigvalue);
    if all(eigvalue>=0)
        y(i) = 1/2*x2'*H{i}*x2 + k{i}'*x2 + d{i};
        if nargout > 2
            grady(:,i) = H{i}*x2 + k{i};
        end
    else
        Qn=zeros(size(H{i}));
        % init index
        if any(eigvalue>=0)
            index=find(eigvalue>=0); 
        else
            index=length(H{i});
        end
        for j=1:index
            Qn=Qn-eigvalue(j)*V(:,j)*V(:,j)';
            Qn=real(Qn);
        end
        Qp=zeros(size(H{i}));
        for j=index:length(H{i})
            Qp=Qp+eigvalue(j)*V(:,j)*V(:,j)';
            Qp=real(Qp);
        end
        temp=1/2*x1'*Qn*x1;
        y(i)=1/2*x2'*Qp*x2+k{i}'*x2+d{i} - (temp+(Qn*x1)'*(x2-x1));
        if nargout > 2
            grady(:,i)=Qp*x2+k{i}- (Qn*x1);
        end
    end
end
yeq = [];
gradyeq = [];
end
