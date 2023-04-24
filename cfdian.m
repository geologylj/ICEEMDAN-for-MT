function out=cfdian(e,rc1,rc2,jj)

for i=1:rc1*5
    out1(i)=e*[1-exp(-i/rc1)];
end
for i=1:rc2*5
    out2(i)=e*exp(-i/rc2);
end

res=zeros(1,rc1*5+rc2*5+1);
out=zeros(1,(rc1+rc2)*10+jj);

for i=1:rc1*5
    res(i+1)=out1(i);
end
for i=1:rc2*5
    res(i+rc1*5+1)=out2(i);
end

for i=1:(rc1+rc2)*5
    out(i+jj)=res(i);
    out(i+(rc1+rc2)*5+jj)=-res(i);
end

