function [X]=matrix_completion_sparse_mex(A,d,Omega,X0,options)
% matrix completion: 
% A- given matrix, each row only has the nonzeros indices;
% d-diagonal indices of A
% Omega- visible indices:consists of n vectors, must be symmetric;
% X0 initial-all zeros

% preprocessing:
mex exactCDmex.cpp
n=size(A,2);
m=0;
lenA=zeros(n,1);
for i=1:n
	%if (length(A{i})>m)
	%	m=length(A{i});
	%end
	lenA(i)=length(A{i});
end

m=max(lenA)

R=A;%remainder part
nA=zeros(n,m);
nO=ones(n,m)*(n+1);
% those vacant spaces, let the index be n+1=> zero value

for i=1:n
	nA(i,1:length(A{i}))=A{i};
	nO(i,1:length(A{i}))=Omega{i};
end
nR=nA;
ind=sub2ind([n m],1:n,d');% get the 1D index of the position of the matrix
X=X0;
x=rand(n,1);
X(:,1)=x;
xp=[x;0];
tmp=xp(nO);

for i=1:m
	%size(tmp(:,i).*x)
	nR(:,i)=nR(:,i)-tmp(:,i).*x;
end

%R=A-X*X';
%R=A;% X is 0 originally

%for t=1:n
%%	size(X(t,:));
%%	size(X(Omega{t},:));
%	R{t}=R{t}-X(Omega{t},:)*X(t,:)';
%end

normA=matrix_norm(nA);
k=size(X0,2);
tic;
X=exactCDmex(nA,nR,nO,X0,lenA,d,normA,options.maxiter);
toc

end

function value=matrix_norm(A)
	%value=0;
	%for i=1:size(A,2)
	%	value=value+norm(A{i})^2;
	%end
	value=norm(A,'fro');
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

function value=roots(a,b)
	delta=4*a.^3+27*b.^2;
	d=(sqrt(delta/27)-b)/2;
	idx=delta>0;
	value=zeros(size(a));
	%if (delta>0)
		value(idx)=nthroot(d(idx),3)+nthroot(-b(idx)/2-sqrt(delta(idx)/27)/2,3);
		% only one possible result
	%else
%a(1:20)
%b(1:20)
		id2=delta<=0;
		a=a(id2);
		b=b(id2);
		delta=delta(id2);
		d=d(id2);

		r=2*nthroot(abs(d),3);
		theta=angle(d)/3;
		%size(r)
		%size(theta)
		%size(id2)
		%total=sum(id2)
		ys=ones(size(r))*inf;
		zs=zeros(size(r));
		for k=0:2
			z=r.*cos(theta+2*k*pi/3);
			tt=z.^4/4+a.*z.^2/2+b.*z;
			idz=tt<ys;
			zs(idz)=z(idz);
			ys(idz)=tt(idz);
		end
		value(id2)=zs;

%		for i=1:sum(id2)
%			zs=0;ys=inf;
%			for k=0:2
%				z=r(i)*cos(theta(i)+2*k*pi/3);
%				if (z^4/4+a(i)*z^2/2+b(i)*z<ys)
%					zs=z;
%					ys=z^4/4+a(i)*z^2/2+b(i)*z;
%				end
%			end
%			value(i)=zs;
%		end
	%endi
	value(1:20)
end

function obj=f(R,normA)
	
	obj=matrix_norm(R)/4/normA;
end


