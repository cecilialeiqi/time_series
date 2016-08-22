function X=factorize(A,k)
% factorize the matrix A with 2 rank k matrices
[V,D]=eig(A);
n=size(A,1);
X=V(:,n-k+1:n)*sqrt(D(n-k+1:n,n-k+1:n));
fprintf('remaining error is: %f\n',norm(A-X*X','fro')/norm(A,'fro'));

end
