
clc; clear all; close all;

% parameters
xa = -15;
ya = 18;
d = 0.90;
m = 85;
I = 1.7;
rho = 20;
Fmax = 180;
Mmax = 10;
Ts = 0.05;
n = 10;


% continuous system

A=[0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 0 0 -rho/m 0 0;
    0 0 0 0 -rho/m 0;
    0 0 0 0 0 -0.01/I];

B=[0 0 0;
    0 0 0;
    0 0 0;
    1/m 0 0;
    0 1/m 0;
    0 0 1/I];

C=[eye(3) zeros(3)];

D=zeros(3);

sys_c = ss(A,B,C,D);


% discrete system

Ad=[1 0 0 Ts 0 0;
    0 1 0 0 Ts 0;
    0 0 1 0 0 Ts;
    0 0 0 1-rho*Ts/m 0 0;
    0 0 0 0 1-rho*Ts/m 0;
    0 0 0 0 0 1-0.01*Ts/I];

Bd=[0 0 0;
    0 0 0;
    0 0 0;
    Ts/m 0 0;
    0 Ts/m 0;
    0 0 Ts/I];

Cd=[eye(3) zeros(3)];

Dd=zeros(3);

sys_d = ss(Ad,Bd,Cd,Dd,0.15);

poles = pole(sys_d)

contr = rank(ctrb(sys_d));

observ = rank(obsv(sys_d));

%  stabilizable
[v,d]=eig(Ad)
for i = 1:size(d)
    rk = [Ad-d(i,i)*eye(size(d)) Bd]
    r_stab(i) = rank(rk)
end
%

%  detectable
[v,d]=eig(Ad)
for i = 1:size(d)
    rk = [Ad-d(i,i)*eye(size(d)); Cd]
    r_det(i) = rank(rk)
end


z = tzero(sys_d)


%% discretization rules
sys = sys_c;
sys_tustin = c2d(sys_c,Ts,'tustin');
sys_zoh = c2d(sys_c,Ts,'zoh');

% backward rect
A_br = inv(Ad);
B_br = A_br * Bd;
C_br = C * A_br;
D_br = D + C_br * Bd;

sys_br = ss(A_br,B_br,C_br,D_br, Ts)

[y,t]=step(sys, 1);
[y1,t1]=step(sys_d, 1);
[y2,t2]=step(sys_tustin, 1);
[y3,t3]=step(sys_zoh, 1);
[y4,t4]=step(sys_br, 1);

figure;

%subplot(3,1,1);
hold on;
plot(t,y(:,:,1));
plot(t1,y1(:,:,1));
% plot(t2,y2(:,:,1));
% plot(t3,y3(:,:,1));
% plot(t4,y4(:,:,1));
% legend('Continuous', 'Euler', 'Tustin', 'ZOH','Back. Rect');
hold off;
xlabel('Time'); ylabel('x'); grid on

figure;
%subplot(3,1,2);
hold on;
plot(t,y(:,:,2));
plot(t1,y1(:,:,2));
% plot(t2,y2(:,:,2));
% plot(t3,y3(:,:,2));
% plot(t4,y4(:,:,2));
hold off;
xlabel('Time'); ylabel('y'); grid on
%legend('Continuous', 'Euler', 'Tustin', 'ZOH','Back. Rect');

%subplot(3,1,3);
figure;
hold on;
plot(t,y(:,:,3));
plot(t1,y1(:,:,3));
% plot(t2,y2(:,:,3));
% plot(t3,y3(:,:,3));
% plot(t4,y4(:,:,3));
hold off;
xlabel('Time'); ylabel('\theta'); grid on
%legend('Continuous', 'Euler', 'Tustin', 'ZOH','Back. Rect');


%% setpoints

load('trajectory_gk.mat');

figure(2)
subplot(3,1,1)
plot(y(1,:));ylabel('x [m]');xlabel('Time [sec]'); grid on

subplot(3,1,2)
plot(y(2,:));ylabel('y [m]');xlabel('Time [sec]');  grid on

subplot(3,1,3)
plot(y(3,:));ylabel('\theta [rad]');xlabel('Time [sec]');  grid on

%% get reference input

AA= sys_tustin.A
BB= sys_tustin.B
CC= sys_tustin.C
DD= sys_tustin.D

n=size(y,2);
W = zeros(3*n,1);

for i = 1:3:3*n-2
    W(i) = y(1,fix(i/3)+1);
    W(i+1) = y(2,fix(i/3)+1);
    W(i+2) = y(3,fix(i/3)+1);
end


Hn1 = kron(eye(n),DD);
for i = 2:n+1
    temp = diag(ones(n-i+1,1), -i+1);
    Markov = CC*AA^(i-2)*BB;
    Hn1 = Hn1 + kron(temp,Markov);
end
Hpr = Hn1'*Hn1;
fpr = -W'*Hn1;
umax =[Fmax Fmax 10]'
b = repmat(umax,n,1);
uref = quadprog(Hpr,fpr',[],[],[],[],-b,b);

uref1 = uref(1:3:3*n-2);
uref2 = uref(2:3:3*n-1);
uref3 = uref(3:3:3*n);

figure;
subplot(3,1,1)
plot(uref1); xlabel('Time [sec]'); ylabel('F_x'); grid on

subplot(3,1,2)
plot(uref2); xlabel('Time [sec]'); ylabel('F_y');  grid on

subplot(3,1,3)
plot(uref3); xlabel('Time [sec]'); ylabel('M'); grid on





%% calcualate reference states

xref = zeros(6,n);
x_prev = zeros(6,1);
xref(:,1) = x_prev;

for i=2:n;
    xref(:,i) = AA*xref(:,i-1)+BB*[uref1(i-1);uref2(i-1);uref3(i-1)];
end

xref3 = xref'

figure;
subplot(3,2,1)
plot(xref3(:,1)); xlabel('Time [sec]'); ylabel('x [m]');  grid on

subplot(3,2,3)
plot(xref3(:,2)); xlabel('Time [sec]'); ylabel('y [m]'); grid on


subplot(3,2,5)
plot(xref3(:,3)); xlabel('Time [sec]'); ylabel('\theta [rad]'); grid on


subplot(3,2,2)
plot(xref3(:,4)); xlabel('Time [sec]'); ylabel('v_x [m\s]');  grid on

subplot(3,2,4)
plot(xref3(:,5)); xlabel('Time [sec]'); ylabel('v_y [m\s]'); grid on


subplot(3,2,6)
plot(xref3(:,6)); xlabel('Time [sec]'); ylabel('\omega'); grid on




figure;
subplot(3,1,1)
hold on;
plot(xref3(:,1)); plot(y(1,:));
xlabel('Time [sec]'); ylabel('x'); grid on
hold off;

subplot(3,1,2)
hold on;
plot(xref3(:,2)); plot(y(2,:));
xlabel('Time [sec]'); ylabel('y'); grid on
hold off;

subplot(3,1,3)
hold on;
plot(xref3(:,3)); plot(y(3,:));
xlabel('Time [sec]'); ylabel('\theta'); grid on
hold off;

xref2 = zeros(6*n,1);

for i = 1:6:6*n-5
    xref2(i) = xref(1,fix(i/6)+1);
    xref2(i+1) = xref(2,fix(i/6)+1);
    xref2(i+2) = xref(3,fix(i/6)+1);
    xref2(i+3) = xref(4,fix(i/6)+1);
    xref2(i+4) = xref(5,fix(i/6)+1);
    xref2(i+5) = xref(6,fix(i/6)+1);
end



%% LQR

time = [0:1:171];
refs = zeros(3,n);
refs(1,:) = uref1';
refs(2,:) = uref2';
refs(3,:) = uref3';

yref_sim = timeseries(y',time);
xref_sim = timeseries(xref',time);
ref_sim = timeseries(refs',time);

R=1e-4*eye(3);
Q1=eye(size(AA));
Q2=[eye(3) zeros(3);
    zeros(3) zeros(3);]

K1=dlqr(AA,BB,Q1,R);
K2=dlqr(AA,BB,Q2,R);




%% LQR cost

xcost = zeros(6*size(y,2),1);
for i = 1:6:6*n-5
    xcost(i) = states.data(fix(i/6)+1,1);
    xcost(i+1) = states.data(fix(i/6)+1,2);
    xcost(i+2) = states.data(fix(i/6)+1,3);
    xcost(i+3) = states.data(fix(i/6)+1,4);
    xcost(i+4) = states.data(fix(i/6)+1,5);
    xcost(i+5) = states.data(fix(i/6)+1,6);
end


ucost = zeros(3*size(y,2),1);
for i = 1:3:3*n-2
    ucost(i) = u_contr.data(fix(i/3)+1,1);
    ucost(i+1) = u_contr.data(fix(i/3)+1,2);
    ucost(i+2) = u_contr.data(fix(i/3)+1,3);
end
temp = repmat(diag(Q1),n,1);
temp = [temp; repmat(diag(R),n,1)];
H = diag(temp);
J = ([xcost;ucost]-[xref2;uref])'*H*([xcost;ucost]-[xref2;uref])




%% MPC

n = size(y,2)
umpc = zeros(3,n-1);
xmpc = zeros(6,n-1);

N = 5*10;
nx = 6;
nu = 3;

Q2 = [eye(3) zeros(3);
    zeros(3) zeros(3);];
Q = Q2;
R = 10^(-4)*eye(3);

% 
%  for iter = 1:20,
%      tic;
temp = repmat(diag(Q),N,1);
temp = [temp; repmat(diag(R),N,1)];
H = 2*diag(temp);

Ai = [zeros(nu*N,nx*N)  eye(nu*N);
    zeros(nu*N,nx*N) -eye(nu*N);];

bi = [repmat([Fmax Fmax Mmax],1,2*N)];
temp1 = diag(ones(N,1));
Ae1 = kron(temp1,eye(nx));
temp2 = diag(ones(N-1,1), -1);
Ae1 = Ae1 + kron(temp2,-AA);
Ae2 = kron(temp1,-BB);
Ae = [Ae1 Ae2];



for k = 1:n-N
    f = -H*[xref2(k*6+1:(k+N)*6); uref((k-1)*3+1:(N+k-1)*3)];
    
    if (k == 1)
        temp = zeros(6,1);
    else
        temp= AA*xmpc(:,k-1) + BB*umpc(:,k-1);
    end
    be=zeros(N*nx,1);
    be(1:6,1) = temp;
    
    xtilde = quadprog(H,f,Ai,bi,Ae,be);
    if (k == n-N)
        xmpc(1,k+1:n) = xtilde(1:6:6*N-5);
        xmpc(2,k+1:n) = xtilde(2:6:6*N-4);
        xmpc(3,k+1:n) = xtilde(3:6:6*N-3);
        xmpc(4,k+1:n) = xtilde(4:6:6*N-2);
        xmpc(5,k+1:n) = xtilde(5:6:6*N-1);
        xmpc(6,k+1:n) = xtilde(6:6:6*N);
        umpc(1,k:n-1) = xtilde(6*N+1:3:9*N-2,1);
        umpc(2,k:n-1) = xtilde(6*N+2:3:9*N-1,1);
        umpc(3,k:n-1) = xtilde(6*N+3:3:9*N,1);
    else
        xmpc(:,k+1) = xtilde(1:6,1);
        umpc(:,k) = xtilde(6*N+1:6*N+3,1);
    end
end

% time = toc;
% t_iter = time/((n-N))
% t_sum = t_sum + t_iter;
% end
% t_avg = t_sum/100;

figure;
clf
subplot(3,1,1)
hold on
plot(uref1); plot(umpc(1,:));
ylabel('Fx'); xlabel('Time');
hold off
subplot(3,1,2)
hold on
plot(uref2); plot(umpc(2,:));
ylabel('Fy'); xlabel('Time');
hold off
subplot(3,1,3)
hold on
plot(uref3); plot(umpc(3,:));
ylabel('M'); xlabel('Time');
legend('ref','MPC');
hold off


Umpc = zeros(3*n,1)
for i = 1:3:3*n-3
    Umpc(i) = umpc(1,fix(i/3)+1);
    Umpc(i+1) = umpc(2,fix(i/3)+1);
    Umpc(i+2) = umpc(3,fix(i/3)+1);
end
outputMPC = Hn1*Umpc;
outMPC1 = outputMPC(1:3:3*n-2);
outMPC2 = outputMPC(2:3:3*n-1);
outMPC3 = outputMPC(3:3:3*n);

figure; %output
subplot(3,1,1)
hold on
plot(y(1,:)); plot(outMPC1);
ylabel('x');xlabel('Time');
hold off
subplot(3,1,2)
hold on
plot(y(2,:)); plot(outMPC2);
ylabel('y');xlabel('Time');
hold off
subplot(3,1,3)
hold on
plot(y(3,:)); plot(outMPC3);
ylabel('\theta');xlabel('Time');
legend('ref', 'MPC')
hold off

figure; % trajectory
hold on
plot(y(1,:),y(2,:));
plot(outMPC1,outMPC2);
legend('ref','MPC')
xlabel('x');
ylabel('y');
hold off


% cost

xcost = zeros(6*n,1);
for i = 1:6:6*n-5
    xcost(i) = xmpc(1,fix(i/6)+1);
    xcost(i+1) = xmpc(2,fix(i/6)+1);
    xcost(i+2) = xmpc(3,fix(i/6)+1);
    xcost(i+3) = xmpc(4,fix(i/6)+1);
    xcost(i+4) = xmpc(5,fix(i/6)+1);
    xcost(i+5) = xmpc(6,fix(i/6)+1);
end

temp = repmat(diag(Q),n,1);
temp = [temp; repmat(diag(R),n,1)];
H = diag(temp);
J = ([xcost;Umpc]-[xref2;uref])'*H*([xcost;Umpc]-[xref2;uref])



%% MPC constraints

Fmax = 10*180

n = size(y,2)
umpc = zeros(3,n-1);
xmpc = zeros(6,n-1);

N = 10/2;
nx = 6;
nu = 3;

Q2 = [eye(3) zeros(3);
    zeros(3) zeros(3);];
Q = Q2;
R = 10^(-4)*eye(3);

xmax = 20.16;
ymax = 16.5;
thetamax = pi/2;
vxmax = 4;
vymax = 4;
omegamax = 3;

temp = repmat(diag(Q),N,1);
temp = [temp; repmat(diag(R),N,1)];
H = 2*diag(temp);


Ai = [ eye(nx*N) zeros(nx*N,nu*N);
    -eye(nx*N) zeros(nx*N,nu*N);
    zeros(nu*N,nx*N)  eye(nu*N);
    zeros(nu*N,nx*N) -eye(nu*N);];

bi = [ repmat([xmax ymax thetamax vxmax vymax omegamax],1,2*N) repmat([Fmax Fmax Mmax],1,2*N) ];
temp1 = diag(ones(N,1));
Ae1 = kron(temp1,eye(nx));
temp2 = diag(ones(N-1,1), -1);
Ae1 = Ae1 + kron(temp2,-AA);
Ae2 = kron(temp1,-BB);
Ae = [Ae1 Ae2];

for k = 1:n-N
    f = -H*[xref2(k*6+1:(k+N)*6); uref((k-1)*3+1:(N+k-1)*3)];
    
    if (k == 1)
        temp = zeros(6,1);
    else
        temp= AA*xmpc(:,k-1) + BB*umpc(:,k-1);
    end
    be=zeros(N*nx,1);
    be(1:6,1) = temp;
    
    xtilde = quadprog(H,f,Ai,bi,Ae,be);
    if (k == n-N)
        xmpc(1,k+1:n) = xtilde(1:6:6*N-5);
        xmpc(2,k+1:n) = xtilde(2:6:6*N-4);
        xmpc(3,k+1:n) = xtilde(3:6:6*N-3);
        xmpc(4,k+1:n) = xtilde(4:6:6*N-2);
        xmpc(5,k+1:n) = xtilde(5:6:6*N-1);
        xmpc(6,k+1:n) = xtilde(6:6:6*N);
        umpc(1,k:n-1) = xtilde(6*N+1:3:9*N-2,1);
        umpc(2,k:n-1) = xtilde(6*N+2:3:9*N-1,1);
        umpc(3,k:n-1) = xtilde(6*N+3:3:9*N,1);
    else
        xmpc(:,k+1) = xtilde(1:6,1);
        umpc(:,k) = xtilde(6*N+1:6*N+3,1);
    end
end

figure;
clf
subplot(3,1,1)
hold on
plot(uref1); plot(umpc(1,:));
ylabel('Fx'); xlabel('Time');
hold off
subplot(3,1,2)
hold on
plot(uref2); plot(umpc(2,:));
ylabel('Fy'); xlabel('Time');
hold off
subplot(3,1,3)
hold on
plot(uref3); plot(umpc(3,:));
ylabel('M'); xlabel('Time');
legend('ref','MPC Constraints');
hold off


Umpc = zeros(3*n,1)
for i = 1:3:3*n-3
    Umpc(i) = umpc(1,fix(i/3)+1);
    Umpc(i+1) = umpc(2,fix(i/3)+1);
    Umpc(i+2) = umpc(3,fix(i/3)+1);
end
outputMPC = Hn1*Umpc;
outMPC1 = outputMPC(1:3:3*n-2);
outMPC2 = outputMPC(2:3:3*n-1);
outMPC3 = outputMPC(3:3:3*n);

figure; %output
subplot(3,1,1)
hold on
plot(y(1,:)); plot(outMPC1);
ylabel('x');xlabel('Time');
hold off
subplot(3,1,2)
hold on
plot(y(2,:)); plot(outMPC2);
ylabel('y');xlabel('Time');
hold off
subplot(3,1,3)
hold on
plot(y(3,:)); plot(outMPC3);
ylabel('\theta');xlabel('Time');
legend('ref', 'MPC Constraints')
hold off

figure; % trajectory
hold on
plot(y(1,:),y(2,:));
plot(outMPC1,outMPC2);
legend('ref','MPC Constraints')
xlabel('x');
ylabel('y');
hold off


% cost

xcost = zeros(6*n,1);
for i = 1:6:6*n-5
    xcost(i) = xmpc(1,fix(i/6)+1);
    xcost(i+1) = xmpc(2,fix(i/6)+1);
    xcost(i+2) = xmpc(3,fix(i/6)+1);
    xcost(i+3) = xmpc(4,fix(i/6)+1);
    xcost(i+4) = xmpc(5,fix(i/6)+1);
    xcost(i+5) = xmpc(6,fix(i/6)+1);
end

temp = repmat(diag(Q),n,1);
temp = [temp; repmat(diag(R),n,1)];
H = diag(temp);
J = ([xcost;Umpc]-[xref2;uref])'*H*([xcost;Umpc]-[xref2;uref])
