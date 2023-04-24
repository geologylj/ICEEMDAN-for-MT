function out = addsj( I,fz,rc1,jj,wz )
%Add charge discharge triangular wave¡£IFor the input signal, rc1 is the 1/4
%peak valley distance, ...
% fz is the amplitude percentage, jj is the spacing, ...
% and wz is the position where the charge discharge triangular wave is added

N=length(I);
asj=zeros(1,N);
big=max(I);   

e=fix(big*fz);    %Conversion of amplitude percentage to actual amplitude
len=rc1*40+jj;     %Length of single positive and negative phase charging and discharging waveform
n=fix(N/len);     

sk=fix(cfdian(e,rc1,rc1*3,jj));
sj=zeros(1,len);  
for i=141:len
    sj(i-140)=sk(i);
end
%--------------------------------------

switch wz  %Determine the position of adding charge discharge triangular waves based on wz
    case 0  
        for i=1:n
            asj((i*len-len+1):i*len)=sj(1:len);
        end
    case 1  
        for i=1:fix(n/3)
            asj((i*len-len+1):i*len)=sj(1:len);
        end
    case 2  
        for i=fix(n/3):fix(n/3*2)
            asj((i*len-len+1):i*len)=sj(1:len);
        end
    case 3  
        for i=fix(n/3*2):n
            asj((i*len-len+1):i*len)=sj(1:len);
        end
    otherwise
        for i=fix(n/2):n
            asj((i*len-len+1):i*len)=sj(1:len);
        end
end
out=I+asj;
