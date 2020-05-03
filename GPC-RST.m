close all, clear all, clc
A=[1 -0.8];
B=[0.4 0.6];
delta = [1 -1];
d=0; N=3; Nu=3;

% Condicao inicial
zi = filter(B,A,0);

% horizontes
N1 = d +1;
N2 = d + N;
At=conv(A,delta);
nA = length(A);
nAt = length(At);
nB = length(B);

% Formulação de F
f=[1 zeros(1, nAt-2)];
for j=1:N,
for i=1:nAt-2
    f(j+1,i)=f(j,i+1)-f(j,1)*At(i+1);
end
f(j+1,nAt-1) = -f(j)*At(nAt);
end
F=f(1+N1:1+N2,:);

% Formulação de E
E=zeros(N2);
e(1)=1;E(1,1)=e(1);
for j=2:N2,
 e(j)=f(j,1);
 E(j,1:j)=e;
end,
E=E(N1:N2,:);


% HORIZONTE DE PREDIÇÃO
% Formulação de Gl
for j=1:N
Gl(j,:) = conv(E(j,:),B);
end
Gl;
ngb=length(Gl);
G = zeros(N,Nu);

% Formulação do G - futuro
for i = 1:N,
 Gp(i,:)= Gl(i,i+1:ngb-N+i); % Gp - passado
end,
G=toeplitz(Gl(N,1:N), zeros(1,Nu)); %G - futuro

% Formulação do controlador
lambda=0.8;
H=2*(G'*G+lambda*eye(length(G'*G)));
K_=inv(G'*G+lambda*eye(length(G'*G)))*G';
K=K_(1,:);
% Resposta no tempo
T=25; %tempo em segundos
Ts = 1; %tempo de amostragem
nts= T/Ts; % numero de amostras

% Utilizando RST
Gp_=K*Gp; F_=K*F;
R_ = [1 Gp_];
S_ = K*F;
T_ = sum(K); % T=Kr--> soma de todos os elementos de K

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
    % Saidas passadas
 yp = [y(k); yp(1:end-1)];

 if nB == 1
 fr = F*yp;
 else
 fr = F*yp + Gp*dup;
 end

 % Variacao de controle
 w = ones(N,1);
 du(k) = K*(w-fr);
 disp(du(k));
 u(k) = du(k) + u0;

 % Resolve a equacao do filtro
 disp(u(k));
 [y(k+1), zf] = filter(B,A,u(k),zi);

 % Atualiza a condicao inicial
 zi = zf;

 if nB>1
 % Atualiza \Delta u(t-1)
 dup = du(k); dup(1:end-1);
 end

 % Atualiza a inicializacao do sinal de controle u(t1)
 u0 = u(k);
end
ygpc=y; ugpc=u;
%Gera figura
t = 0:Ts:T;
plot(t,yrst(2:end),'b+',t,ygpc(1:end-1),'b'), hold on
plot(t,urst(2:end),'g+',t,ugpc,'g'), grid on
legend('GPC_{RST}','GPC')

% plot(t,yrst(2:end),'b'), hold on
% plot(t,urst(2:end),'g'), grid on
% legend('GPC_{RST}')
