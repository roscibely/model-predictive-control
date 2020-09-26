clear; clc;
%% Model parameters
A1 = [0.751 -0.0014; 0.986 0.063];
A2 = [0.751 -0.0014; 9.864 0.189];
A3 = [-0.136 -0.014; 9.864 0.189];
A4 = [-0.136 -0.014; 98.644 1.451];
B = [0.15 0; 0 -0.912];
N = 10;
umax = [0.5; 1];
x1set = [1.4e-1, 9e-2, 7e-2, 6e-2, 5e-2, 3e-2, 2e-2, 1.7e-2, 0.01, 0.001];
ymax = 1;
tsim = 11;
%% LMIs
%Parametros Qc Rc
Qc = eye(2); Rc = 0.2*eye(2);
n=size(B,1); m=size(B,2);
Q=sdpvar(n,n,'symmetric');
X=sdpvar(m,m,'symmetric');
Y=sdpvar(m,n,'full');
gama=sdpvar(1);
xk = sdpvar(2,1);

m111=[Q (A1*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m112=[(A1*Q+B*Y) Q zeros(n) zeros(n,m)];
m211=[Q (A2*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m212=[(A2*Q+B*Y) Q zeros(n) zeros(n,m)];
m311=[Q (A3*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m312=[(A3*Q+B*Y) Q zeros(n) zeros(n,m)];
m411=[Q (A4*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m412=[(A4*Q+B*Y) Q zeros(n) zeros(n,m)];
m3=[(sqrt(Qc)*Q) zeros(n) gama*eye(n) zeros(n,m)];
m4=[(sqrt(Rc)*Y) zeros(n,m)' zeros(n,m)' gama*eye(m)];
M1=[m111;m112;m3;m4]; %M1
M2=[m211;m212;m3;m4]; %M2
M3=[m311;m312;m3;m4]; %M3
M4=[m411;m412;m3;m4]; %M4
LMIs=[Q >= 0, M1 >= 0,  M2 >= 0,  M3 >= 0,  M4 >= 0, X>=0, gama >=0]; %restrições LMIs
LMIs=[LMIs, [1 xk';xk Q]>=0]; %restrição no estado
LMIs=[LMIs, [X Y;Y' Q]>=0, X(1,1) <= umax(1,1)^2, X(2,2)<=umax(2,1)^2 ]; %restrição no sinal de controle
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
% Obtendo o conjunto de valores para montagem da lookuptable
F11=[]; F12=[];
%% Processo iterativo off-line
for k=1:length(x1set)
model = optimizer(LMIs, gama,ops,xk,{X,Y,Q});
x=[x1set(k);0];
QY = model{x};
disp(QY{1});
QQ(:,:,k)=QY{3};
YY(:,:,k)=QY{2};
F(:,:,k) = YY(:,:,k)*inv(QQ(:,:,k));
F11=[F11; F(1,1,k) ];
F12=[F12; F(1,2,k) ];
LMIs=[LMIs, Q<=QQ(:,:,k)]; %acrescimo da restrição do algoritmo offline
end
%%
figure(1);
for k=1:length(x1set)
[xx(k,:), yy(k,:)] = elipse_matrix(QQ(:,:,k),40); 
hold on, plot(xx(k,:),yy(k,:),'k'), hold on,
xlabel('x_1'), ylabel('x_2'), title('Elipsoides de Estabilidade');
end
%% Condições iniciais
x = [0.1; 2];
u=[-0.1;0.95];
alpha = 1.1;
bet =1.1;
A = [0.85-0.0986*alpha -0.0014*alpha; 0.9864*alpha*bet 0.0487+0.01403*alpha*bet];
% tempo de simulação
t = [0:0.1:((tsim-1)/10)];
for k=1:tsim
x(:,k+1) = A*x(:,k) + B*u(:,k);
u(:, k+1) = F(:,:,5)*(x(:,k+1));
end
%%
figure(2)
subplot(121), plot(t,x(1,1:tsim)/0.1,'b',t,x(2,1:tsim)/2,'k');
title('Estados'); xlabel('Tempo (seg)'); ylabel('x');
legend('C_A','T')
subplot(122),
plot(t,u(:, 1:tsim),'k'); xlabel('Tempo (seg)'); ylabel('u(volts)');
title('Sinal de controle')
figure(3)
subplot(121), plot(x1set,F11); title('F(1,1)')
subplot(122), plot(x1set,F12); title('F(1,2)')





