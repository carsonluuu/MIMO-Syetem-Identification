%----------------------------%
%Task Two
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
figure(1);
subplot(311)
plot(t,u1,'b*','LineWidth',2)
ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])
subplot(312)
plot(t,y11,'r*','LineWidth',2)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])
subplot(313)
plot(t,y21,'r*','LineWidth',2)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])
figure(2);
subplot(311)
plot(t,u2,'b*','LineWidth',2)
ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])
subplot(312)
plot(t,y12,'r*','LineWidth',2)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])
subplot(313)
plot(t,y22,'r*','LineWidth',2)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])

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
scatter(i,hsv100');
hold on;
%% ns Problem
ns = 7;
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
% tzero(A,B,C,[]);
%% 
D = 0;
iopzmap(idss(A,B,C,D,'Ts',1/40))
% hold on
% plot(real(poles),imag(poles),'ro');
% hold on;
% plot(real(zeros),imag(zeros),'rx');
% axis square;
% axis([-3 1.5 -2.5 2.5]);
% axis([-1.5 1.5 -3 1.5 ],'square')
% command line:
% 40*log(-0.214 + 0.889i);

% -0.214 + 0.889i
% -0.214 - 0.889i
% -0.188 + 0.892i
% -0.188 - 0.892i
% 0.77
% 0.826
% 0.869

% -2.6310
% 0.9988
% -0.6078 + 0.6740i
% -0.6078 - 0.6740i
% -0.3241
% D = 0;
% 
% S2=[ A    B; 
%     -C   -D];
% I = [eye(size(A))  zeros(7,2);
%      zeros(2,7) zeros(2,2)];
% [v2 d2] = eig(S2,I);
% 
% plot(eig(A) ,'x');
% hold on
% plot(diag(d2),'o')
% xlim([-3 2])
% 
% xlabel('real')
% ylabel('complex')
% axis equal
% grid on
%  alpha=0:pi/20:2*pi;    %] 
% 
%  R=1;                   %
%  x=R*cos(alpha); 
% 
%  y=R*sin(alpha); 

%  plot(x,y,'m-.') 

%% H for H1,2,3,4
% H_1_100 = zeros(1,100);
% H_1_n   = zeros(1,100);
% H_2_100 = zeros(1,100);
% H_2_n   = zeros(1,100);
% H_3_100 = zeros(1,100);
% H_3_n   = zeros(1,100);
% H_4_100 = zeros(1,100);
% H_4_n   = zeros(1,100);
% 
% for k=1:100
%     H_1_100(k) = y11(42+k);
%     H_1_n(k)      = y11(141+k);
%     H_2_100(k) = y21(42+k);
%     H_2_n(k)      = y21(141+k);
%     H_3_100(k) = y12(42+k);
%     H_3_n(k)      = y12(141+k);
%     H_4_100(k) = y22(42+k);
%     H_4_n(k)      = y22(141+k);
% end
% H1 = hankel(H_1_100, H_1_n);
% H2 = hankel(H_2_100, H_2_n);
% H3 = hankel(H_3_100, H_3_n);
% H4 = hankel(H_4_100, H_4_n);
% 
% [~,S1,~] = svd(H1);
% [~,S2,~] = svd(H2);
% [~,S3,~] = svd(H3);
% [~,S4,~] = svd(H4);
% hsv100 = zeros(100,1);
% for i = 1:100;
%     hsv100_1(i) = (S1(i,i));
%     hsv100_2(i) = (S2(i,i));
%     hsv100_3(i) = (S3(i,i));
%     hsv100_4(i) = (S4(i,i));
% end
% i = 1:100;
% figure(7)
% subplot(221)
% semilogy(i,hsv100_1','o')
% xlabel('$Singular Value Index$','FontSize',14,'Interpreter','Latex');
% ylabel('$Hankel Singular Value H(1,1)$','FontSize',14,'Interpreter','Latex');
% subplot(222)
% semilogy(i,hsv100_2','o')
% xlabel('$Singular Value Index$','FontSize',14,'Interpreter','Latex');
% ylabel('$Hankel Singular Value H(2,1)$','FontSize',14,'Interpreter','Latex');
% subplot(223)
% semilogy(i,hsv100_3','o')
% xlabel('$Singular Value Index$','FontSize',14,'Interpreter','Latex');
% ylabel('$Hankel Singular Value H(1,2)$','FontSize',14,'Interpreter','Latex');
% subplot(224)
% semilogy(i,hsv100_4','o')
% xlabel('$Singular Value Index$','FontSize',14,'Interpreter','Latex');
% ylabel('$Hankel Singular Value H(2,2)$','FontSize',14,'Interpreter','Latex');
% hold on;