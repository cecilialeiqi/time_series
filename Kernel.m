%generating the kernel matrix
function [D]=Kernel(X,n)
	% use the first n users, and generate the time series vs time series similarity matrix
	D=zeros(n,n);
	
    nrm=zeros(n,1);
	for i=1:n
		nrm(i)=dtw(X{i},zeros(1,size(X{i},2)));
	end
	
	for i=1:n
		D(i,i)=1;
		for j=i+1:n
			D(i,j)=(nrm(i)^2+nrm(j)^2-dtw(X{i},X{j})^2)/(nrm(i)^2+nrm(j)^2);
			D(j,i)=D(i,j);
		end
	end
    % now D is the kernel matrix required

end

