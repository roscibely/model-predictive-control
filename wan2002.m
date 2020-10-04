% Robust output feedback model predictive control using
%  off-line linear matrix inequalities Wan and Kothare 2002
clear all, close all, clc
%% Exemple 1
% Parameters
alpha=[1 5];
beta=[1 5];
%Vertix
A1 = [0.85-0.0986*alpha(1) -0.0014*alpha(1);
    0.9864*alpha(1)*beta(1) 0.0487+0.01403*alpha(1)*beta(1)];
A2 = [0.85-0.0986*alpha(1) -0.0014*alpha(1);
    0.9864*alpha(1)*beta(2) 0.0487+0.01403*alpha(1)*beta(2)];
A3 = [0.85-0.0986*alpha(2) -0.0014*alpha(2);
    0.9864*alpha(2)*beta(1) 0.0487+0.01403*alpha(2)*beta(1)];
A4 = [0.85-0.0986*alpha(2) -0.0014*alpha(2);
    0.9864*alpha(2)*beta(2) 0.0487+0.01403*alpha(2)*beta(2)];
B = [0; -0.912];
C=[0 1];
% x-set
x1set=[1, 0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.035, 0.01];
tsim = 30;
% Control signal
umax=0.5;
%Weight matrix
Qc = [0 0; 0 1]; Rc = 2e-5;
%% LMIS
n=size(B,1); m=size(B,2);
Q=sdpvar(n,n,'symmetric');
X=sdpvar(m,m,'symmetric');
Y=sdpvar(m,n,'full');
gama=sdpvar(1);
xk = sdpvar(2,1);
% Matrix M1
m11=[Q (A1*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m12=[(A1*Q+B*Y) Q zeros(n) zeros(n,m)];
m3=[(sqrt(Qc)*Q) zeros(n) gama*eye(n) zeros(n,m)];
m4=[(sqrt(Rc)*Y) zeros(n,m)' zeros(n,m)' gama*eye(m)];
% Matrix M2
m21=[Q (A2*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m22=[(A2*Q+B*Y) Q zeros(n) zeros(n,m)];
%Matrix M3
m31=[Q (A3*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m32=[(A3*Q+B*Y) Q zeros(n) zeros(n,m)];
%Matrix M4
m41=[Q (A4*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m42=[(A4*Q+B*Y) Q zeros(n) zeros(n,m)];
M1=[m11;m12;m3;m4]; %M1
M2=[m21;m22;m3;m4]; %M2
M3=[m31;m32;m3;m4]; %M3
M4=[m41;m42;m3;m4]; %M4
LMIs=[Q >= 0, M1 >= 0, M2>=0, M3>=0, M4>=0, X>=0, gama >=0];
LMIs=[LMIs, [1 xk';xk Q]>=0]; %state constraints
LMIs=[LMIs, [X Y;Y' Q]>=0, X<=umax*umax]; %control signal constraints
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
%lookuptable
F11=[]; F12=[];
%% Off-line robust observer design
p = 0.1;
Pe = sdpvar(2,2, 'symmetric');
Ye = sdpvar(2,1);
Le =eye(2);
Lmi= [Pe>=0, [p^2*Pe-Le (Pe*A1-Ye*C)'; Pe*A1-Ye*C Pe]>=0];
Lmi = [Lmi, [p^2*Pe-Le (Pe*A2-Ye*C)'; Pe*A2-Ye*C Pe]>=0];
Lmi = [Lmi, [p^2*Pe-Le (Pe*A3-Ye*C)'; Pe*A3-Ye*C Pe]>=0];
Lmi = [Lmi, [p^2*Pe-Le (Pe*A4-Ye*C)'; Pe*A4-Ye*C Pe]>=0];
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
optimize(Lmi,-trace(Pe),ops);
Lp = inv(value(Pe))*value(Ye);

%% MPC off-line
for k=1:length(x1set)
    model = optimizer(LMIs, gama,ops,xk,{X,Y,Q});
    x=[x1set(k);0];
    QY = model{x};
    QQ(:,:,k)=QY{3};
    YY(:,:,k)=QY{2};
    F(:,:,k) = YY(:,:,k)*inv(QQ(:,:,k));
    F11=[F11; F(1,1,k) ];
    F12=[F12; F(1,2,k) ];
    LMIs=[LMIs, Q<=QQ(:,:,k)];
end
% figure(1)
for k=1:length(x1set)
    [xx(k,:), yy(k,:)] = elipse_matrix(QQ(:,:,k),40);
    hold on, plot(xx(k,:),yy(k,:),'k'), hold on,
    xlabel('x_1'), ylabel('x_2'), title('Elipsoides de Estabilidade');
end
%% Condições iniciais
N=5;
x = [0.1; 2]; u = 0;
alpha=linspace(1,5,tsim);
beta=linspace(1,5,tsim);
xhat=[0; 0];
X=[x; xhat];
for k=1:tsim
    A =  [0.85-0.0986*alpha(k) -0.0014*alpha(k);
        0.9864*alpha(k)*beta(k) 0.0487+0.01403*alpha(k)*beta(k)];
    Apoly = [A B*F(:,:,N); Lp*C A1+B*F(:,:,N)-Lp*C];
    X(:,:,k+1) = Apoly*X(:,:,k);
    z= X(:,:,k+1);
    u(k+1) = F(:,:,N)*(x(:,k));
    xhat(:,k+1)= [z(1,1); z(2,1)];
    x(:,k+1)=[z(3,1); z(4,1)];
    y(k)=C*x(:,k);
end
% time
t = [0:0.1:((tsim-1)/10)];
figure(2)
subplot(311), plot(t,x(1,1:tsim),'b-+',t,xhat(1,1:tsim),'k-+'); legend('x_1','\hat{x_1}');
subplot(312), plot(t,x(2,1:tsim),'b-+',t,xhat(2,1:tsim),'k-+'); legend('x_2','\hat{x_2}');
subplot(313), plot(t,u(1:end-1),'k-+'); xlabel('Tempo (seg)'); ylabel('u');  title('Sinal de controle')
figure(3)
subplot(121), plot(x1set,F11); title('F(1,1)')
subplot(122), plot(x1set,F12); title('F(1,2)')