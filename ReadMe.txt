Step 1: run transform.m directly: it read in the dynamic data and transform it
into n matrices: X{i} denotes the i-th time series

	(run transform_multiplelines.m: it only transform the time series with
	 multiple lines into X)

Step 2: for dense similarity generation: call D=Kernel(X,n). n being the
number of users you want to use.
	for sparse similarity generation: call [D,Omega,d]=Kernel_sparse(X,m). m being the
	number of pairs one wants to observe. It outputs a sparse matrix D, and
	the indices of observed pairs in each row. 

Step 3: with dense similarity A, call X=matrix_completion(A,Omega,X0,options);
	with sparse input: call matrix_completion_sparse(A,d,Omega,X0,options)



