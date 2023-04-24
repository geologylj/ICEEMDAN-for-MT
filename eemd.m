function [modos its]=eemd(x,Nstd,NR,MaxIter)

%   modos=eemd(x,Nstd,NR,MaxIter)
%  [modos its]=eemd(x,Nstd,NR,MaxIter)


desvio_estandar=std(x);
x=x/desvio_estandar;
xconruido=x+Nstd*randn(size(x));
[modos, o, it]=emd(xconruido,'MAXITERATIONS',MaxIter);
modos=modos/NR;
iter=it;
if NR>=2
    for i=2:NR
        xconruido=x+Nstd*randn(size(x));
        [temp, ort, it]=emd(xconruido,'MAXITERATIONS',MaxIter);
        temp=temp/NR;
        lit=length(it);
        [p liter]=size(iter);
        if lit<liter
            it=[it zeros(1,liter-lit)];
        end;
        if liter<lit
            iter=[iter zeros(p,lit-liter)];
        end;
        
        iter=[iter;it];
        
        [filas columnas]=size(temp);
        [alto ancho]=size(modos);
        diferencia=alto-filas;
        if filas>alto
            modos=[modos; zeros(abs(diferencia),ancho)];
        end;
        if alto>filas
            temp=[temp;zeros(abs(diferencia),ancho)];
        end;
        
        modos=modos+temp;
    end;
end;
its=iter;
modos=modos*desvio_estandar;