%generating the kernel matrix from Euclidean distance
function [D,Omega,d]=Kernel_sparse_Euclid(X,n,m)
	% use the first n users, and generate approximately m pairs among them
	%X=random_TS(n,d);
	D={};
	Omega={};
	d=zeros(n,1);
	idi=randsample(n,m/2,'true');%allow repetition
	idj=randsample(n,m/2,'true');
	v=zeros(m/2,1);
    nrm=zeros(n,1);
	for i=1:n
		nrm(i)=norm(X{i},'fro');
	end
	
	for k=1:m/2
		%v(i)=0;i
		i=idi(k);
		j=idj(k);
		v(k)=(nrm(i)^2+nrm(j)^2-EUCdist(X{i},X{j})^2)/(nrm(i)^2+nrm(j)^2);
	end
    
	size(v)

	idx=[idi;idj;(1:n)'];
	idy=[idj;idi;(1:n)'];
	v=[v;v;ones(n,1)];
	size(idx)
	size(idy)
	size(v)
	S=sparse(idx,idy,v,n,n);	
    % has to do this step, since there might be duplicates
    [row,col,v]=find(S);
	for i=1:n
		id=find(col==i);
		Omega{i}=row(id);
		D{i}=v(id);
		d(i)=find(Omega{i}==i);
	end

    % now D is the kernel matrix required

end

function value=EUCdist(X,Y)
	if (size(X,2)~=size(Y,2))
		print('Error! Two time series not in same space');
		exit(0);
	end
	s=min(size(X,1),size(Y,1));
	value=0;
	for i=1:s
		value=value+norm(X(i,:)-Y(i,:))^2;
	end
	value=sqrt(value);
end
