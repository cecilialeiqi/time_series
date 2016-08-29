%run transform.m
n=size(X,2);
[D,Omega,d]=Kernel_sparse(X,n,n*50);
X0=zeros(n,15);
options.maxiter=15;
Xt=matrix_completion_sparse(D,d,Omega,X0,options);
