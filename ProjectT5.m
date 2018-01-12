%----------------------------%
%Task Five
%----------------------------%
load u_rand.mat
y1 = u_rand.Y(3).Data/2;
y2 = u_rand.Y(4).Data/2;
u1 = u_rand.Y(1).Data/2;
u2 = u_rand.Y(2).Data/2;

ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;

p = 12000;
u = [u1; u2];
y = [y1; y2];

a = zeros(1,401);b = zeros(1,401);c = zeros(1,401);d = zeros(1,401);

%%Problem 1
Ruu = [0 0 ; 0 0];
for q = 1:24001
    Ruu = Ruu + y(:,q)*(y(:,q))';
end 
y_rms = sqrt(diag((1/(2*p))*Ruu));
RMS = norm(y_rms)

%%problem 2
sys =ss(A,B,C,D,ts);
Gc = gram(sys,'c');
Go = gram(sys,'o');
PH2_1 = norm(sqrt(diag(B'*Go*B)))
PH2_2 = norm(sqrt(diag(C*Gc*C')))

% problem 3
 %experimental norm
 PH2_3  = 0;
for i = 1:401
    n = norm([y11(i) y12(i);y21(i) y22(i)],'fro');
    PH2_3 = PH2_3+n^2;
end 
PH2_3   = sqrt(PH2_3);
PH2_3   = norm(PH2_3);