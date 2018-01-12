%----------------------------%
%Task Four
%----------------------------%
load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;
ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;


p = 12000;
u = [u1; u2];
y = [y1; y2];
mean(u1)
mean(u2)

a = zeros(1,401);b = zeros(1,401);c = zeros(1,401);d = zeros(1,401);

for k = 1:401
    sm = zeros(2);
    for q = 202:2*p-202
          sm = sm + u(:,q+k-201)*(u(:,q))' ;   
    end
    matr = 1/(2*p)*sm;
    a(k) = matr(1,1);b(k) = matr(1,2);c(k) = matr(2,1);d(k) = matr(2,2);
end
%% 

figure(2222)
x = linspace(-5,5,401);
subplot(221)
plot (x,a,'g','LineWidth',0.5)
ylabel('$R_{uu}(1,1)$','FontSize',10,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
subplot(222)
plot (x,b,'g','LineWidth',0.5)
ylabel('$R_{uu}(1,2)$','FontSize',10,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
subplot(223)
plot (x,c,'g','LineWidth',0.5)
ylabel('$R_{uu}(2,1)$','FontSize',10,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
subplot(224)
plot (x,d,'g','LineWidth',0.5)
ylabel('$R_{uu}(2,2)$','FontSize',10,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');


%% problem 4
a = zeros(1,89);b = zeros(1,89);c = zeros(1,89);d = zeros(1,89);

for k = 1:89
    sm = zeros(2);
    for q = 89:2*p-89
          sm = sm + y(:,q+k-8)*(u(:,q))' ;   
    end
    matr = 1/(2*p)*sm;

    a(k) = matr(1,1)/var(u1);b(k) = matr(1,2)/var(u2);c(k) = matr(2,1)/var(u1);d(k) = matr(2,2)/var(u2);
end
%% 

x = linspace(-0.2,2,89);
figure(4444)
subplot(411)
plot (x+0.025,a,'go',t,y11,'r*','LineWidth',1)
ylabel('$y1(u1)$','FontSize',12,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
h41 = legend ('Cross-Correlation','Actual');
set(h41,'FontSize',9);
xlim([-0.2 2])

subplot(413)
plot (x+0.025,c,'go',t,y21,'r*','LineWidth',1)
ylabel('$y2(u1)$','FontSize',12,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
h42 = legend ('Cross-Correlation','Actual');
set(h42,'FontSize',9);
xlim([-0.2 2])

subplot(412)
plot (x+0.025,b,'go',t,y12,'r*','LineWidth',1)
ylabel('$y1(u2)$','FontSize',12,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
h43 = legend ('Cross-Correlation','Actual');
set(h43,'FontSize',9);
xlim([-0.2 2])

subplot(414)
plot (x+0.025,d,'go',t,y22,'r*','LineWidth',1)
ylabel('$y2(u2)$','FontSize',12,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
h44 = legend ('Cross-Correlation','Actual');
set(h44,'FontSize',9);
xlim([-0.2 2])