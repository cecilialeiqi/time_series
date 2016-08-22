M=csvread('./data/time_series_onehot_data.csv');
%M=csvread('time_series_numerical_part.csv');
id=int64(M(:,1));

%label=int64(M(:,size(M,2)));
% 
head=1;
X={};
count=1;
tail=head+1;
while tail<=size(M,1)
	tail=head+1;
	while tail<=size(M,1) && id(tail)==id(head)
		tail=tail+1;
	end
	X{count}=M(head:tail-1,2:size(M,2));% from 2, the ID is eliminated
	label(count)=M(head,size(M,2));
	head=tail;
	count=count+1
end
% distinguish static data from dynamic ones


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

static_data=zeros(size(X,2),sum(static));
for i=1:size(static_data,1)
	static_data(i,:)=X{i}(1,logical(static));
end

%csvwrite('static_data.csv',static_data);

n=size(X,2);
for i=1:n
	X{i}=X{i}(:,dynamic);
end


totalmax=zeros(1,size(X{1},2));
for i=1:n 
	
    if (size(X{i},1)==1)
		currentmaxvalue=X{i};
	else
		currentmaxvalue=max(X{i})
	end
	totalmax=max([currentmaxvalue;totalmax]);
end
totalmax(totalmax==0)=1;

for i=1:n
    for j=1:size(X{i},1)
        X{i}(j,:)=X{i}(j,:)./totalmax;
    end
end

