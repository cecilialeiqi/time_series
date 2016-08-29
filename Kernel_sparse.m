%generating the kernel matrix
function [D,Omega,d]=Kernel_sparse(X,n,m)
	% use the first n users, and generate approximately m pairs among them
	%X=random_TS(n,d);
	mex dtw_c.c;
	D={};
	Omega={};
	d=zeros(n,1);

	id2d=randsample(n*n,2*m,'false');
	idi=floor((id2d-1)/n)+1;
	idj=id2d-n*(idi-1);
	id=find(idi<idj);
	idi=idi(id);
	idi=idi(1:(m-n)/2);
	idj=idj(id);
	idj=idj(1:(m-n)/2);
	% still need to remove same indices...
	
	v=zeros(floor((m-n)/2),1);
    nrm=zeros(n,1);
	tic;
	for i=1:n
		nrm(i)=dtw_c(X{i},zeros(1,size(X{i},2)),1);
	end
	
	for k=1:(m-n)/2
		%v(i)=0;i
		i=idi(k);
		j=idj(k);
		v(k)=(nrm(i)^2+nrm(j)^2-dtw_c(X{i},X{j},10)^2)/(nrm(i)^2+nrm(j)^2);
		%v(k)=(nrm(i)^2+nrm(j)^2-dtw(X{i},X{j})^2)/2/nrm(i)/nrm(j)^2;
	end
    toc

	size(v)

	col=[idi;idj;(1:n)'];
	row=[idj;idi;(1:n)'];
	v=[v;v;ones(n,1)];
	m=size(col)
	size(row)
	size(v)
	[col,I]=sort(col);
	row=row(I);
	v=v(I);
	%S=sparse(idx,idy,v,n,n);	
    % has to do this step, since there might be duplicates
    %[row,col,v]=find(S);
	start=1;
	nd=1;
	for i=1:n
	    while (nd<=m && col(nd)==i)
		nd=nd+1;
	    end
	    Omega{i}=row(start:nd-1);
	    D{i}=v(start:nd-1);
	    d(i)=find(Omega{i}==i);
	    start=nd;
	end

	%for i=1:n
	%	id=find(col==i);
	%	Omega{i}=row(id);
	%	D{i}=v(id);
	%	d(i)=find(Omega{i}==i);
	%end

    % now D is the kernel matrix required

end

