%%this is a test line
clc;clear;close all

%% INIT %%
IC = load('init-test_small.mat').IC;
% IC = [-483.767045759555,-584.583619137290,2413.04233229495,-19.8710395039958,-24.4932506265223,-82.2778907669142,1849.31049678031];
t0 = 0;
tf = 100;

%% DE solve %%
% run < parpool(8) > in cmd window before using < parfor >
% run < delete(gcp('nocreate')) > to close the parallel pool, must be done
% to ensure RAM space is released

N = 10000;
tSpan = linspace(t0,tf,N);

del_m = zeros(1, length(IC));

bf = 5;
% 1: smooth
% 2: glideslope
% 3: 1-step
% 4: flat top, 1-step
% 5: flat top, 2-step

tic
for cl = 2:3
% 1: classical ZEM-ZEV
% 2: self adjusting ZEM-ZEV
% 3: new ogl
parfor i = 1:length(IC) % set this to < parfor > when needed
    init = IC(i,:);

    [t,x] = ode45(@(t,x) odemain(t,x,bf,cl), tSpan, init);

    %% PLOTS %%
    N = length(x);

% %     % TRAJECTORY %
% %     figure(1)
% %     plot3(x(:,1),x(:,2),x(:,3), 'LineWidth', 1.5);
% %     hold on
% %     plot3(x(1,1),x(1,2),x(1,3),'ks')
% %     hold on
% %     plot3(x(N,1),x(N,2),x(N,3),'ro')
% %     hold on
% %     grid on
% % 
    % POSITION %
    figure(2)
    subplot(3,1,1)
    plot(t,x(:,1))
    hold on
    grid on

    subplot(3,1,2)
    plot(t,x(:,2))
    hold on
    grid on

    subplot(3,1,3)
    plot(t,x(:,3))
    hold on
    grid on

    % VELOCITY %
    figure(3)
    subplot(3,1,1)
    plot(t,x(:,4))
    hold on
    grid on

    subplot(3,1,2)
    plot(t,x(:,5))
    hold on
    grid on

    subplot(3,1,3)
    plot(t,x(:,6))
    hold on
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


    
    if bf == 1
        bx = 0:0.1:650;
        by = 0:0.1:650;
        bzx = zeros(1,6501);
        bzy = zeros(1,6501);
        
        k0 = [0;0;0];
        k1 = [0.17;0.17;0.17]*1000;
        k2 = [0;0;1];
        k3 = [6;6;6];
        
        bzx = ((bx-k0(1))/k1(1)).^(k3(1)) - k2(1);
        bzy = ((by-k0(2))/k1(2)).^(k3(2)) - k2(2);
    elseif bf == 2
        bx = 0:0.1:1000;
        by = 0:0.1:1000;
        bzx = zeros(1,1001);
        bzy = zeros(1,1001);
        
        k0 = [0;0;0];
        k1 = [1.5;1.5;1.5];
        k2 = [0;0;1];
        k3 = [1;1;1];
        
        bzx = ((bx-k0(1))/k1(1)).^(k3(1)) - k2(1);
        bzy = ((by-k0(2))/k1(2)).^(k3(2)) - k2(2);
    elseif bf == 3
        bx = 0:0.1:1000;
        by = 0:0.1:1000;
        bzx = zeros(1,10001);
        bzy = zeros(1,10001);
        
        for j = 1:width(bx)
            if bx(j) > 300
                k0 = [300;300;300];
                k1 = [178;178;178];
                k2 = [-200;-200;1];
                k3 = [6;6;6];
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            else
                k0 = [0;0;0];
                k1 = [231;231;231];
                k2 = [0;0;1];
                k3 = [20;20;20];
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            end
        end
        for k = 1:width(by)
            if by(k) > 300
                k0 = [300;300;300];
                k1 = [178;178;178];
                k2 = [-200;-200;1];
                k3 = [6;6;6];
                bzy(k) = ((by(k)-k0(2))/k1(2)).^(k3(2)) - k2(2);
            else
                k0 = [0;0;0];
                k1 = [231;231;231];
                k2 = [0;0;1];
                k3 = [20;20;20];
                bzy(k) = ((by(k)-k0(2))/k1(2)).^(k3(2)) - k2(2);
            end
        end
    elseif bf == 4
        bx = 0:0.1:1000;
        by = 0:0.1:1000;
        bzx = zeros(1,10001);
        bzy = zeros(1,10001);
        
        for j = 1:width(bx)
            if bx(j) > 300
                k0 = 0*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-200;-200;201];
                k3 = 1*ones(3,1);
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            else
                k0 = [0;0;0];
                k1 = [230;230;230];
                k2 = [0;0;1];
                k3 = [20;20;20]; 
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            end
        end
        for k = 1:width(by)
            if by(k) > 300
                k0 = 0*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-200;-200;201];
                k3 = 1*ones(3,1);
                bzy(k) = ((by(k)-k0(2))/k1(2)).^(k3(2)) - k2(2);
            else
                k0 = [0;0;0];
                k1 = [230;230;230];
                k2 = [0;0;1];
                k3 = [20;20;20]; 
                bzy(k) = ((by(k)-k0(2))/k1(2)).^(k3(2)) - k2(2);
            end
        end

    elseif bf == 5
        bx = 0:0.1:2000;
        by = 0:0.1:2000;
        bzx = zeros(1,20001);
        bzy = zeros(1,20001);        
        for j = 1:width(bx)
            if bx(j) > 1000
                k0 = -1000*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-1500;-1500;1501];
                k3 = ones(3,1);
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            elseif ((bx(j) <= 1000) && (bx(j) > 600))
                k0 = 600*ones(3,1);
                k1 = 142*ones(3,1);
                k2 = [-1000;-1000;1001];
                k3 = 6*ones(3,1);
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            else
                k0 = zeros(3,1);
                k1 = 426*ones(3,1);
                k2 = zeros(3,1);
                k3 = 20*ones(3,1);
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            end
        end
        
        for j = 1:width(by)
            if by(j) > 1000
                k0 = -1000*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-1500;-1500;1501];
                k3 = ones(3,1);
                bzy(j) = ((by(j)-k0(2))/k1(2)).^(k3(2)) - k2(2);
            elseif ((by(j) <= 1000) && (by(j) > 600))
                k0 = 600*ones(3,1);
                k1 = 142*ones(3,1);
                k2 = [-1000;-1000;1001];
                k3 = 6*ones(3,1);
                bzy(j) = ((by(j)-k0(2))/k1(2)).^(k3(2)) - k2(2);
            else
                k0 = zeros(3,1);
                k1 = 426*ones(3,1);
                k2 = zeros(3,1);
                k3 = 20*ones(3,1);
                bzy(j) = ((by(j)-k0(2))/k1(2)).^(k3(2)) - k2(2);
            end
        end 
    end

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
    plot(bx, bzx, '--');
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
    plot(-by, bzy, '--');
    hold on

    grid on

    del_m(i) = x(1,7) - x(N-1,7);

end
figure(7)
xlabel('Test Case #')
ylabel('\Deltam')

scatter(1:length(IC),del_m, 'LineWidth', 1);
hold on
grid on

mu = mean(del_m);
disp(['mu = ', num2str(mu)]);

sd = std(del_m);
disp(['sd = ', num2str(sd)]);

end
toc


