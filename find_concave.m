function [minimum,maximum]=find_concave(vector)

[m,n]=size(vector);
[A,IX]=sort(vector, 'descend');

C1=(1:1:m)'-IX;

Index=find(C1==0);
a=length(Index);
minimum=sum(Index==(1:a)')+1;

maximum=m-sum(Index==((m-a+1):m)'); 




