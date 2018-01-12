%----------------------------%
%Task Six 
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

% Define frequency response function
%------function continues-Time------

% function out = hinf61(sys,ttol,~)
% 
% H = 0;
% I = eye(length(D));
% ww = 0; %number of lo with good
% q = 0;
% N = 1;
% Nmax = 10000;
% zer = 1e-5;
% while N < Nmax
%     c = (ll+up)/2
%     Dl = (c)^2*I-D'*D;
%     a = A+B*inv(Dl)*D'*C;
%     b = -B*inv(Dl)*B';
%     c = C'*C+C'*D*inv(Dl)*D'*C;
%     d = -A'-C'*D*inv(Dl)*B';
%     Aclp = [a b; c d];
%     reale = real(eig(Aclp));
%     image = imag(eig(Aclp));
% 
%     for ii = (1: length(real(eig(Aclp))))
%         if abs(reale(ii)) < zer && image(ii) > zer || (b-a)/2 < tolerance
%             display(gam)
%             return 
%         end 
%     end
%     N = N+1;
% 
%     ll = c;
% end 
% Ga = unique(H);
% length(Ga)
% P       =C*inv((1j*om*I)-A)*B+D;
% step    = 0.1 ; %step size

% for lo      = [ll:step:ul]
% Dl = lo^2*I-D'D;
% a  = A+B*inv(Dl)*D'*C;
% b  = -B*inv(Dl)*B';
% c  = C'*C+C'*D*inv(Dl)*D'*C;
% d  = -A'-C'*D*inv(Dl)*B';
% 
% Aclp    = [a b; c d];
% 
% end 
% Ac = -inv(I + A)*(I-A)
% Bc = sqrt(2)*inv(I+A)*B
% Cc = sqrt(2)*inv(I+A)
% Dc = D - C*inv(I+A)*B


% narginchk(1,3);
% ni = nargin;
% 
% if (ni<2 || isempty(ttol))
%     ttol = 1e-3;
% end
% 
% 
% % SYSTEM matrices
% [a,b,c,d] = unpck(sys);
% sys=ss(a,b,c,d);
%     if max(real(eig(a)))< 0 % stable:compute gain
%         [f,w] = getPeakGain(sys,ttol);
%         out = [f f*(1+ttol) w];
%     end
% 


%------function Discrete-Time------

% function out = hinf62(sys,ttol,~)
% 
% narginchk(1,3);
% ni = nargin;
% 
% if (ni<2 || isempty(ttol))
%     ttol = 1e-3;
% end
% 
% 
% % SYSTEM matrices
% [a,b,c,d] = unpck(sys);
% sys=ss(a,b,c,d);
%     if max(real(eig(a)))< 0 % stable:compute gain
%         [f,w] = getPeakGain(sys,ttol);
%         out = [f f*(1+ttol) w];
%     end
% 

%% 

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

Hsv = zeros(2,201);
for n = 1:2:401
    Hsv(:,0.5*(n-1)+1) = svd(fr_s(:,n:n+1));
end
%% 
w = 0:0.1*2*pi:20*2*pi;
semilogx(w,20*log10(Hsv(1,:)),w,20*log10(Hsv(2,:)))
% hold on
% sigma(ss(A,B,C,0,1/40))