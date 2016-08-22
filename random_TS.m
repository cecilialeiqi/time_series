function [X]=random_TS(n,d)
%generating random time series, d is the same for every user, but l is different
for i=1:n
	l=randsample(d,1)+10;
	X{i}=rand(l,d);
end
