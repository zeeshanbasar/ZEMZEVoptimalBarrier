clc;clear;close all

%% LANDER IC %%
r0 = [-0.6 0.3 2.5]'*1000';
v0 = [-10 12 -75]';
m_wet = 1905;
t0 = 0;
tf = 100;

%% DE solve %%
N = 10000;
tSpan = linspace(t0,tf,N);
init = [r0' v0' m_wet']';

[t,x] = ode45(@odefun, tSpan, init);

%% PLOTS %%
N = length(x);

% TRAJECTORY %
figure(1)
plot3(x(:,1),x(:,2),x(:,3), 'LineWidth', 1.5);
hold on
plot3(x(1,1),x(1,2),x(1,3),'ks')
hold on
plot3(x(N,1),x(N,2),x(N,3),'ro')
hold off
grid on

% POSITION %
figure(2)
subplot(3,1,1)
plot(t,x(:,1))
grid on

subplot(3,1,2)
plot(t,x(:,2))
grid on

subplot(3,1,3)
plot(t,x(:,3))
grid on

% VELOCITY %
figure(3)
subplot(3,1,1)
plot(t,x(:,4))
grid on

subplot(3,1,2)
plot(t,x(:,5))
grid on

subplot(3,1,3)
plot(t,x(:,6))
grid on

% MASS %
figure(4)
plot(t,x(:,7))
grid on

% % THRUST %
% A = [0 0 1]';
% del = 1;
% phi = del^2/3;
% c = 500;
% tgo = tf - t;
% r = [x(:,1) x(:,2) x(:,3)];
% v = [x(:,4) x(:,5) x(:,6)];
% ZEM = rd' - (r + v.*tgo + 0.5*g'.*tgo.^2);
% ZEV = vd' - (v + g'.*tgo);
% a_av = c*(x(:,3).^2 - phi).*(tgo.^2)./(24*(x(:,3).^2 + phi).^2);
% T = ((6*ZEM./tgo.^2) - (2*ZEV./tgo) + a_av).*x(:,7);
% 
% figure(5)
% subplot(3,1,1)
% plot(t,T(:,1))
% grid on
% 
% subplot(3,1,2)
% plot(t,T(:,2))
% grid on
% 
% subplot(3,1,3)
% plot(t,T(:,3))
% grid on

k1 = [0.17;0.17;0.17]*1000;
k2 = [0;0;0];
k3 = [6;6;6];

bx = 0:0.1:650;
by = 0:0.1:650;
bzx = (bx/k1(1)).^(k3(1)) - k2(1);
bzy = (by/k1(2)).^(k3(2)) - k2(2);

figure(5)
xlabel('x (m)', 'FontSize', 14)
ylabel('z (m)', 'FontSize', 14)

plot(x(:,1),x(:,3), 'LineWidth', 1)
hold on
plot(x(1,1),x(1,3),'ks')
hold on
plot(x(N,1),x(N,3),'ro')
hold on
plot(-bx, bzx, '--');
hold on

grid on


figure(6)
xlabel('y (m)', 'FontSize', 14)
ylabel('z (m)', 'FontSize', 14)

plot(x(:,2),x(:,3), 'LineWidth', 1)
hold on
plot(x(1,2),x(1,3),'ks')
hold on
plot(x(N,2),x(N,3),'ro')
hold on
plot(by, bzy, '--');
hold on

grid on
