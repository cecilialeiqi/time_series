M=csvread('data/time_series_onehot_data.csv');
%M=csvread('data/time_series_numerical_part.csv');
id=int64(M(:,1));
%labels=int64(M(:,105));
%label=int64(M(:,size(M,2)));
% 
head=1;
X={};
count=1;
tail=head+1;
%Xt=[];
%n=size(M,2)
%M=[M(:,2:104),M(:,106:n)];
n=size(M,2)

while tail<=size(M,1)
	tail=head+1;
	while tail<=size(M,1) && id(tail)==id(head)
		tail=tail+1;
	end
	if (tail-head>1)
	    X{count}=M(head:tail-1,2:n);% from 2, the ID is eliminated
	    %label(count)=labels(head);
		%Xt(count,1:n)=mean(X{count});
	    count=count+1
	end
	head=tail;
end
% distinguish static data from dynamic ones
%Xt(1,:)
%csvwrite('mean_all_feature_multiplelines.csv',Xt);


static=ones(size(X{1},2),1);

for i =1:size(X,2)
	if size(X{i},1)>1
		for j=1:size(X{i},2)
			if (length(unique(int64(X{i}(:,j))))>1)
				static(j)=0;
			end
		end
	end
end
dynamic=logical(1-static);

static(86)=0;
dynamic(86)=0; %not fully observed feature, need to elimiate
%dynamic(102)=0;
%dynamic(34)=1;

static_data=zeros(size(X,2),sum(static));
for i=1:size(static_data,1)
	static_data(i,:)=X{i}(1,logical(static));
end

%csvwrite('static_data_multiplelines.csv',static_data);
%csvwrite('label_multiplelines.csv',label);
n=size(X,2);
X{1}=X{1}(:,dynamic);
D=[ones(size(X{1},1),1),X{1}];
for i=2:n
	X{i}=X{i}(:,dynamic);
	D=[D;[i*ones(size(X{i},1),1),X{i}]];
end
csvwrite('time_series_dynamic_multiplelines.csv',D);

% also, we need to normalize each column


totalmax=zeros(1,size(X{1},2));


totalmean=zeros(1,size(X{1},2));


for i=1:n
	totalmean=totalmean+mean(X{i});
end
totalmean=totalmean/n;


for i=1:n
	%for j=1:size(X{i},1)
	%	X{i}(j,:)=X{i}(j,:)-totalmean;
	%end
	currentmaxvalue=max(X{i})
    totalmax=max([currentmaxvalue;totalmax]);
end

%size(totalmax)
totalmax(totalmax==0)=1;

for i=1:n
	for j=1:size(X{i},1)
		X{i}(j,:)=X{i}(j,:)./totalmax;
	end
end

% the calculated X is only dynamic data, and normalized 



