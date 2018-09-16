close all, clear all, clc
%system
A=[1 -0.8];
B=[0.4 0.6];
delta = [1 -1];
d=0; N=3; Nu=3;
%Initial condition
zi = filter(B,A,0);
%Horizons
N1 = d +1;
N2 = d + N;
At=conv(A,delta);
nA = length(A);
nAt = length(At);
nB = length(B);
%Formulation of F
f=[1 zeros(1, nAt-2)];
for j=1:N,
for i=1:nAt-2
    f(j+1,i)=f(j,i+1)-f(j,1)*At(i+1);
end
f(j+1,nAt-1) = -f(j)*At(nAt);
end
F=f(1+N1:1+N2,:);
%Formulation of E
E=zeros(N2);
e(1)=1;E(1,1)=e(1);
for j=2:N2,
 e(j)=f(j,1);
 E(j,1:j)=e;
end,
E=E(N1:N2,:);
%Predition horizon
for j=1:N
Gl(j,:) = conv(E(j,:),B);
end
Gl;
ngb=length(Gl);
G = zeros(N,Nu);
%Formulation of G - future
for i = 1:N,
 Gp(i,:)= Gl(i,i+1:ngb-N+i); % Gp - past
end,
G=toeplitz(Gl(N,1:N), zeros(1,Nu)); %G - future
%Control Formulation
lambda=0.8;
H=2*(G'*G+lambda*eye(length(G'*G)));
K_=inv(G'*G+lambda*eye(length(G'*G)))*G';
K=K_(1,:);
T=25; %time s
Ts = 1; %sample time
nts= T/Ts; %number of sample

%With RST configuration
Gp_=K*Gp; F_=K*F;
R_ = [1 Gp_];
S_ = K*F;
T_ = sum(K); % T=Kr
Bz1=conv(B,[0 1]);
num=conv(Bz1,T_);
num_u=conv(A,T_);
den=conv(R_,At)+conv(Bz1,S_);
yrst=filter(num,den,[0 ones(1,nts+1)]);
urst=filter(num_u,den,[0 ones(1,nts+1)]);
if nB > 1
 dup = zeros(nB-1,1);
end
u0 = 0;
y=0; yp = zeros(nA, 1);
for k = 1:nts+1
 yp = [y(k); yp(1:end-1)];
 if nB == 1
 fr = F*yp;
 else
 fr = F*yp + Gp*dup;
 end
 w = ones(N,1);
 du(k) = K*(w-fr);
 disp(du(k));
 u(k) = du(k) + u0;
 disp(u(k));
 [y(k+1), zf] = filter(B,A,u(k),zi);
 zi = zf;
 if nB>1
 dup = du(k); dup(1:end-1);
 end
 u0 = u(k);
end
ygpc=y; ugpc=u;
t = 0:Ts:T;
plot(t,yrst(2:end),'b+',t,ygpc(1:end-1),'b'), hold on
plot(t,urst(2:end),'g+',t,ugpc,'g'), grid on
legend('GPC_{RST}','GPC')
