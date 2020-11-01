function y=process(x,h,N)
L=length(h);Lx1=length(x);
H=fft(h,N);
M=N-L+1;
P=Lx1/M;
P=ceil(P);
diff=P*M-Lx1;
x=[x zeros(1,diff)];
Lx2=length(x);
xb=zeros(1,M);
yb=zeros(1,M+L-1);
y=zeros(1,Lx2+L-1);
for i=0:P-1
    xb=x(i*M+1:(i+1)*M);
    yb=cc(xb,H);
    yb=yb(1:M+L-1);
    if i==0
        y(1:length(yb))=yb;
    else
        y(1+i*M : length(yb)+i*M)=y(1+i*M : length(yb)+i*M)+ yb;
    end
end
y=y(1:Lx1+L-1);
end