%----------------------------%
%Task One 
%----------------------------%
%% Part 0: Plot system responses and input
clear
close all
load u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0); %%% find index where pulse occurs
load u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;
%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));
%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);
ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets
t = [0:N-1]*ts - 1;
% figure(1);
% subplot(311)
% plot(t,u1,'b*','LineWidth',2)
% ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 1.1])
% subplot(312)
% plot(t,y11,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 0.1])
% subplot(313)
% plot(t,y21,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 0.1])
% figure(2);
% subplot(311)
% plot(t,u2,'b*','LineWidth',2)
% ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 1.1])
% subplot(312)
% plot(t,y12,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 0.1])
% subplot(313)
% plot(t,y22,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% axis([-0.2 2 -0.1 0.1])

%% Part 1: Model identification
H = zeros(200,200);
n = 1;
for k=1:100
    j = 1;
    for i = k:k+99
        H(n,j) = y11(41+i);
        H(n+1,j) = y12(41+i);
        H(n,j+1) = y21(41+i);
        H(n+1,j+1) = y22(41+i);
        j = j+2;    
    end
    n = n+2;
end

[U,D,V] = svd(H);
hsv100 = zeros(40,1);
for i = 1:40;
    hsv100(i) = (D(i,i));
end
i = 1:40;
% semilogy(i,hsv100','o');
hold on;
%% ns Problem
ns = 6;
for i = 1:ns
    Si(i,i) = hsv100(i);
end
U = U(:,1:ns);
V = V(:,1:ns);
O = U;
C = Si*V';

H_hat = zeros(200,200);
n = 1;
for k=1:100
    j = 1;
    for i = k:k+99
        H_hat(n,j) = y11(42+i);
        H_hat(n+1,j) = y12(42+i);
        H_hat(n,j+1) = y21(42+i);
        H_hat(n+1,j+1) = y22(42+i);
        j = j+2;    
    end
    n = n+2;
end

A = U'*H_hat*V/Si;
Aem = max(abs(eig(A)));
B = C(1:ns,1:2);
C = O(1:2,1:ns);

h = zeros(2,100);
    for k = 1:100
        a = C*A^(k-1)*B;
        h(1,2*k-1) = a(1,1);
        h(2,2*k-1) = a(2,1);
        h(1,2*k) = a(1,2);
        h(2,2*k) = a(2,2);

        qa(k+41)=a(1,1);
        qb(k+41)=a(2,1);
        qc(k+41)=a(1,2);
        qd(k+41)=a(2,2);
    end

%Plot Response   
figure(1);
subplot(321)
plot(t,u1,'b*','LineWidth',2);  ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');    grid on;    axis([-0.2 2 -0.1 1.1])
subplot(323)
plot(t,y11,'r*','LineWidth',2)
hold on

plot(t(1:141),qa,'cx','LineWidth',1.5);
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex'); grid on
axis([-0.2 2 -0.1 0.1])

subplot(325)
plot(t,y21,'r*','LineWidth',2)
hold on
plot(t(1:141),qc,'cx','LineWidth',1.5);
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
axis([-0.2 2 -0.1 0.1])

subplot(322)
plot(t,u2,'b*','LineWidth',2)
ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])
subplot(324)
plot(t,y12,'r*','LineWidth',2)
hold on
plot(t(1:141),qb,'cx','LineWidth',1.5);
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])
subplot(326)
plot(t,y22,'r*','LineWidth',2)
hold on
plot(t(1:141),qd,'cx','LineWidth',1.5);
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
axis([-0.2 2 -0.1 0.1])
%% task3&4
clear i;

ts = 1/40;
I = eye(ns);
fr_s = [];

for w = 0:0.1*2*pi:20*2*pi 
    fr_s = [fr_s (C*inv(exp(w*ts*i)*I-A)*B)];
end

y_f11 = [];
y_f12 = [];
for j=1:2:401
    y_f11 = [y_f11 fr_s(1,j)];
    y_f12 = [y_f12 fr_s(2,j)];
end

y_f21 = [];
y_f22 = [];
for j=2:2:402
    y_f21 = [y_f21 fr_s(1,j)];
    y_f22 = [y_f22 fr_s(2,j)];
end

y11f = fft(y11)./fft(u1);
N = length(y11f);
om1 = [0:N-1]/(ts*N);

y21f = fft(y21)./fft(u1);
N = length(y21f);
om2 = [0:N-1]/(ts*N);

y12f = fft(y12)./fft(u2);
N = length(y12f);
om3 = [0:N-1]/(ts*N);

y22f = fft(y22)./fft(u2);
N = length(y22f);
om4 = [0:N-1]/(ts*N);

%%
tf = om1*2*pi;
w = 0:0.1*2*pi:20*2*pi;

subplot(421)
loglog(w,abs(y_f11),'LineWidth',1.5);
grid on
hold on
loglog(tf,abs(y11f));xlim([0,125]);
grid on
hold on
title('y_{11} magnitude')

subplot(423)
loglog(w,abs(y_f12),'LineWidth',1.5);
grid on
hold on
loglog(tf,abs(y12f));xlim([0,125]);
grid on
hold on
title('y_{12} magnitude')

subplot(425)
loglog(w,abs(y_f21),'LineWidth',1.5);
grid on
hold on
loglog(tf,abs(y21f));xlim([0,125]);
grid on
hold on
title('y_{21} magnitude')

subplot(427)
loglog(w,abs(y_f22),'LineWidth',1.5);
grid on
hold on
loglog(tf,abs(y22f));xlim([0,125]);
grid on
hold on
title('y_{22} magnitude')

subplot(422)
semilogx(w,180/pi*phase(y_f11),'LineWidth',1.5);
grid on
hold on
semilogx(tf,180/pi*phase(y11f),'LineWidth',1.5);xlim([0,125]);
grid on
hold on
title('y_{11} Phase Angel')

subplot(424)
semilogx(w,180/pi*phase(y_f12),'LineWidth',1.5);
grid on
hold on
semilogx(tf,180/pi*phase(y12f),'LineWidth',1.5);xlim([0,125]);
grid on
hold on
title('y_{12} Phase Angel')

subplot(426)
semilogx(w,180/pi*phase(y_f21));
grid on
hold on
semilogx(tf,180/pi*phase(y21f),'LineWidth',1.5);xlim([0,125]);
grid on
hold on
title('y_{21} Phase Angel')

subplot(428)
semilogx(w,180/pi*phase(y_f22));
grid on
hold on
semilogx(tf,180/pi*phase(y22f),'LineWidth',1.5);xlim([0,125]);
grid on
hold on
title('y_{22} Phase Angel')
