function [modes,its]=iceemdan(x,Nstd,NR,MaxIter,SNRFlag)

%modes=iceemdan(x,Nstd,NR,MaxIter,SNRFlag)
%[modes its]=iceemdan(x,Nstd,NR,MaxIter,SNRFlag)

x=x(:)';
desvio_x=std(x);
x=x/desvio_x;

aux=zeros(size(x));
iter=zeros(NR,round(log2(length(x))+5));

white_noise = cell(1,NR);
for i=1:NR
    white_noise{i}=randn(size(x));%creates the noise realizations
end

modes_white_noise = cell(1,NR);
for i=1:NR
    modes_white_noise{i}=emd(white_noise{i});%calculates the modes of white gaussian noise
end

for i=1:NR %calculates the first mode
    xi=x+Nstd*modes_white_noise{i}(1,:)/std(modes_white_noise{i}(1,:));
    [temp, ~, it]=emd(xi,'MAXMODES',1,'MAXITERATIONS',MaxIter);
    temp=temp(1,:);
    aux=aux+(xi-temp)/NR;
    iter(i,1)=it;
end

modes= x-aux; %saves the first mode
medias = aux;
k=1;
aux=zeros(size(x));
es_imf = min(size(emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter)));

while es_imf>1 %calculates the rest of the modes
    for i=1:NR
        tamanio=size(modes_white_noise{i});
        if tamanio(1)>=k+1
            noise=modes_white_noise{i}(k+1,:);
            if SNRFlag == 2
                noise=noise/std(noise); %adjust the std of the noise
            end
            noise=Nstd*noise;
            try
                [temp,~,it]=emd(medias(end,:)+std(medias(end,:))*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
            catch    
                it=0;
                %disp('catch 1 '); disp(num2str(k))
                temp=emd(medias(end,:)+std(medias(end,:))*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
            end
            temp=temp(end,:);
        else
            try
                [temp, ~, it]=emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter);
            catch
                temp=emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter);
                it=0;
                %disp('catch 2 sin ruido')
            end
            temp=temp(end,:);
        end
        aux=aux+temp/NR;
    iter(i,k+1)=it;    
    end
    modes=[modes;medias(end,:)-aux];
    medias = [medias;aux];
    aux=zeros(size(x));
    k=k+1;
    es_imf = min(size(emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter)));
end
modes = [modes;medias(end,:)];
modes=modes*desvio_x;
its=iter;
