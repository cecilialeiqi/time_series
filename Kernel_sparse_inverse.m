%generating the kernel matrix
function [D,Omega,d]=Kernel_sparse_inverse(X,n,m)
	% use the first n users, and generate approximately m pairs among them
	%X=random_TS(n,d);
	D={};
	Omega={};
	d=zeros(n,1);
	idi=randsample(n,m/2-n,'true');%allow repetition
	idj=randsample(n,m/2-n,'true');
	v=zeros(m/2-n,1);
  	
	for k=1:m/2-n
		%v(i)=0;i
		i=idi(k);
		j=idj(k);
		v(k)=dtw(X{i},X{j});
	end
    
	size(v)
	zeroidx=v==0;
	posidx=v>0;
	minvalue=min(v(posidx));
	v(zeroidx)=1;
	v(posidx)=minvalue./v(posidx);

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

