function val=pmean(num)
% This function gives the mean of non zero and non(NaN) numbers in a vector

num=num(:);
for i=1:numel(num)
    if isnan(num(i))
        num(i)=0;
    end
end
ind=num==0;
num(ind)=[];
val=mean(num);