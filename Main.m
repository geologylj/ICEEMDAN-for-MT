clc
clear all
close all
%%  Simulated charge discharge triangular wave
x1=100*randn(1,5000);
x=addsj(x1,5,8,200,4);
xx=[-x(1501:3000) x(1001:1700) -x(1201:2000) x(2001:4000)];
x1=[-x1(1501:3000) x1(1001:1700) -x1(1201:2000) x1(2001:4000)];
x2=xx-x1;
xx=load('OriginalSignal.dat');
xx=xx-mean(xx);
fids1=fopen('OriginalSignal.dat','wt');
fprintf(fids1,'%10.0f%10.0f',xx);
fclose(fids1);
[S,F,T,P] = spectrogram(xx,256,250,256,150);
figure(1)
subplot(221);
sf= surf(T,F,10*log10(P))
sf.EdgeColor = 'none';
% view(0,90);
% colorbar
axis tight
xlabel('时间 (s)')
ylabel('频率(Hz)')
figure(2)
subplot(221)
fs=150;t=(0:length(xx)-1)/fs;
plot(t,xx,'r');ylabel('幅值');legend('原始信号');
%% Simulated square wave signal
%% square wave
x3=load('moni_fangbo.dat');x3=x3*0.5;x2=x2*2;
%% Pulses and harmonics
x4=load('moni_xiebo.dat');
xn=[x2(1:1000),-3500,x4(1:1000),3500,x2(1:1000),x3(1001:2000),-x2(1:998)]
xx=xx+xn*0.5;

[S,F,T,P] = spectrogram(xx,256,250,256,150);
figure(1)
subplot(222);
sf= surf(T,F,10*log10(P))
sf.EdgeColor = 'none';
% view(0,90);
% colorbar
axis tight
xlabel('time (s)')
ylabel('frequency(Hz)')
figure(2)
subplot(222)
fs=150;t=(0:length(xx)-1)/fs;
plot(t,xx,'r');ylabel('Amplitude');legend('Simulated noise data');
%% mrsvd denoising
data_length=200;
L=fix(length(xx)/data_length);
for j=1:L
  signal=xx((j-1)*data_length+1:j*data_length);
  sig=signal;
for Q=1:5
 H=[sig(1:end-4);sig(2:end-3);sig(3:end-2);sig(4:end-1);sig(5:end)]; 
   [u,s,v]=svd(H,'econ');
%% Matrix reconstruction signal reconstruction
s(5,5)=0;s(4,4)=0;s(2,2)=0;s(3,3)=0;
H1=u*s*v';
b = rot90(H1);
r = [];
for i = 1 : sum(size(b))-1
k = i - size(b,1);
Diag = diag(b,k);
r = [r;mean(Diag)];
end   
 r=r';  sig=r;
end
da(1,(j-1)*data_length+1:j*data_length)=signal-r;
end
[S,F,T,P] = spectrogram(da,256,250,256,150);
figure(1)
subplot(223);
sf= surf(T,F,10*log10(P))
sf.EdgeColor = 'none';
% view(0,90);
% colorbar
axis tight
xlabel('时间 (s)')
ylabel('频率(Hz)')
figure(2)
subplot(223)
fs=150;t=(0:length(xx)-1)/fs;
plot(t,da,'r');ylabel('幅值');legend('原始信号');
%% iceemdan denoising
data_length=1000;
L=fix(length(da)/data_length);
for C=1:L
     signal=da((C-1)*data_length+1:C*data_length);
  sig=signal;
modes=iceemdan(signal,0.2,100,1000,1);
[a,b]=size(modes);
N=length(modes);
K=zeros(1,N);
for i=1:a
    for j=1:b
    Sj(i)=median(abs(modes(i,:)-median(abs(modes(i,:)))))/0.6745;
    Tj(i)=Sj(i)*sqrt(2*log10(N));
    w(i,j)=exp(-exp(Tj(i)*abs(modes(i,j)/Sj(i))-Tj(i)));
        if abs(modes(i,j))>Tj(i)
            mode(i,j)=modes(i,j)*w(i,j);
        else
            mode(i,j)=modes(i,j);
        end
     end
    K=K+mode(i,:);
end
y(1,(C-1)*data_length+1:C*data_length)=K;

end
 [S,F,T,P] = spectrogram(y,256,250,256,150);
figure(1)
subplot(224);
sf= surf(T,F,10*log10(P));
sf.EdgeColor = 'none';
axis tight
xlabel('time (s)')
ylabel('frequency(Hz)')
figure(2)
subplot(224)
fs=150;t=(0:length(xx)-1)/fs;
plot(t,y,'r');ylabel('amplitude');legend('Original signal');
