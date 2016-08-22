function X=matrix_completion_sparse(A,d,Omega,X0,options)
% matrix completion: 
% A- given matrix, each row only has the nonzeros indices;
% d-diagonal indices of A
% Omega- visible indices:consists of n vectors, must be symmetric;
% X0 initial-all zeros

n=size(A,2);
X=X0;
%R=A-X*X';
R=A;% X is 0 originally
for t=1:n
	size(X(t,:));
	size(X(Omega{t},:));
	R{t}=R{t}-X(Omega{t},:)*X(t,:)';
end

normA=matrix_norm(A);
k=size(X0,2);
for iter=1:options.maxiter
	fprintf('#%d, observed error=%f\n',iter,f(R,normA));
	for i=1:k
	
		x=X(:,i);
		%R=R+x*x';
		for t=1:n
			R{t}=R{t}+x(t)*x(Omega{t});
		end
		for j=1:n
			id=Omega{j};
			p=norm(x(id))^2-x(j)^2-R{j}(d(j));
			q=-x(id)'*R{j}+R{j}(d(j))*x(j);% need to save diagonal part for efficient call
			x(j)=root(p,q);
		end
		%R=R-x*x';
		for t=1:n
			R{t}=R{t}-x(t)*x(Omega{t});
		end
		X(:,i)=x;
	end
end

end

function value=matrix_norm(A)
	value=0;
	for i=1:size(A,2)
		value=value+norm(A{i})^2;
	end

end
function value=root(a,b)
	delta=4*a^3+27*b^2;
	d=(sqrt(delta/27)-b)/2;
	if (delta>0)
		value=nthroot(d,3)+nthroot(-b/2-sqrt(delta/27)/2,3);
		% only one possible result
	else
		r=2*nthroot(norm(d),3);
		theta=angle(d)/3;
		zs=0;ys=inf;
		for k=0:2
			z=r*cos(theta+2*k*pi/3);
			if (z^4/4+a*z^2/2+b*z<ys)
				zs=z;
				ys=z^4/4+a*z^2/2+b*z;
			end
		end
		value=zs;
	end
end

function obj=f(R,normA)
	
	obj=matrix_norm(R)/4/normA;
end


