clc;clear;close all

%% INIT %%
IC = load('init_small-fin.mat').IC;

t0 = 0;
tf = 100;

%% DE solve %%
% run < parpool(8) > in cmd window before using < parfor >
% run < delete(gcp('nocreate')) > to close the parallel pool, must be done
% to ensure RAM space is released

% warning('off','all')

N = tf*100;
tSpan = linspace(t0,tf,N);
options = odeset('Events',@terminate);

del_m = zeros(length(IC),1);

bf = 5;
% 1: smooth
% 2: glideslope
% 3: 1-step
% 4: flat top, 1-step
% 5: flat top, 2-step

l1 = 6;

tic
for cl = 1
% 3: classical ZEM-ZEV
% 2: self adjusting ZEM-ZEV
% 1: new ogl
for i = 1:length(IC) % set this to < parfor > when needed
%     i = 1;
    init = (IC(i,:));
% for l1 = 1:5:40

    [t,x] = ode45(@(t,x) odefin(t,x,tf,bf,cl,init), tSpan, init, options);
    [~,T] = cellfun(@(t,x) odefin(t,x,tf,bf,cl,init), num2cell(t), num2cell(x,2),'uni',0);
    %% PLOTS %%
    N = length(x);

% %    % TRAJECTORY %
% %     figure(1)
% %     plot3(x(:,1),x(:,2),x(:,3), 'LineWidth', 1.5);
% %     hold on 
% %     plot3(x(1,1),x(1,2),x(1,3),'ks')
% %     hold on
% %     plot3(x(N,1),x(N,2),x(N,3),'ro')
% %     hold on
% %     grid on
% %
% %    % POSITION %
% %     figure(2)
% %     subplot(3,1,1)
% %     plot(t(1:100:end),x(1:100:end,1))
% %     hold on
% %     grid on
% % 
% %     subplot(3,1,2)
% %     plot(t(1:100:end),x(1:100:end,2))
% %     hold on
% %     grid on
% % 
% %     subplot(3,1,3)
% %     plot(t(1:100:end),x(1:100:end,3))
% %     hold on
% %     grid on
% % 
   % VELOCITY %
    figure(3)
    subplot(3,1,1)
    plot(t(1:end),x(1:end,4))
    hold on
    grid on

    subplot(3,1,2)
    plot(t(1:end),x(1:end,5))
    hold on
    grid on

    subplot(3,1,3)
    plot(t(1:end),x(1:end,6))
    hold on
    grid on

   % PLANAR TRAJECTORY %

    figure(5)
    plot(x(1:end,1),x(1:end,3), 'LineWidth', 1)
    hold on
    plot(x(1,1),x(1,3),'ks')
    hold on
    plot(x(N,1),x(N,3),'ro')
    hold on
    grid on


    figure(6)
    plot(x(1:end,2),x(1:end,3), 'LineWidth', 1)
    hold on
    plot(x(1,2),x(1,3),'ks')
    hold on
    plot(x(N,2),x(N,3),'ro')
    hold on
    grid on

    figure(12)
    plot(x(N,1),x(N,2),'ks')
    hold on
    grid on
    axis tight

    figure(13)
    scatter(i,x(end,6))
    hold on
    grid on

    % ACCELERATION %
    T = [T{:,1:end}]';

    figure(9)
    plot(t(1:length(T)),T(1:end,1))
    hold on
    grid on

    figure(10)
    plot(t(1:length(T)),T(1:end,2))
    hold on
    grid on

    figure(11)
    plot(t(1:length(T)),T(1:end,3))
    hold on
    grid on


% %     del_m(i,cl) = x(1,7) - x(N-1,7);
% %     del_m_ap(i,cl) = x(1,7) - x(N-1,7);

%% SAVE TRAJ AND THRUST FOR PLOTS %%

% % filename_x = append('[[final]]\mat files\x',num2str(i),'.mat');
% % filename_T = append('[[final]]\mat files\T',num2str(i),'.mat');
% % 
% % var_x = append('x');
% % var_T = append('T');
% % 
% % save(filename_x,var_x);
% % save(filename_T,var_T);
end

% FUEL STATS %

% % filename_del = append('[[final]]\mat files\del_m',num2str(cl),'_ap.mat');
% % 
% % var_del = append('del_m_ap');
% % 
% % save(filename_del,var_del);
% % 
% % figure(7)
% % scatter(1:length(IC),del_m(:,cl), 'LineWidth', 1);
% % hold on
% % grid on
% % 
% % mu(cl) = mean(del_m(:,cl));
% % disp(['mu = ', num2str(mu)]);
% % 
% % sd(cl) = std(del_m(:,cl));
% % disp(['sd = ', num2str(sd)]);

end
toc

% % figure(8)
% % boxchart(del_m(:,1:cl));
% % grid on

%% PLOT POSTPROCESSING %%
if bf == 1
    bx = 0:0.1:650;
    by = 0:0.1:650;

    k0 = [0;0;0];
    k1 = [0.17;0.17;0.17]*1000;
    k2 = [0;0;1];
    k3 = [6;6;6];

    bzx = ((bx-k0(1))/k1(1)).^(k3(1)) - k2(1);
    bzy = ((by-k0(2))/k1(2)).^(k3(2)) - k2(2);
elseif bf == 2
    bx = 0:0.1:1000;
    by = 0:0.1:1000;

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
        h = [500;1000];
        w = [600;1000];
        delta = 10;

        if bx(j) > w(2)
            alpha = deg2rad(0.05);

            k3 = ones(3,1);
            k0 = w(2)*ones(3,1);
            k1 = tan((pi/2) - alpha)*ones(3,1);
            k2 = [-h(2);-h(2);h(2) + delta];

            bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);

        elseif ((bx(j) <= w(2)) && (bx(j) > w(1)))

            k3 = 6*ones(3,1);
            k0 = w(1)*ones(3,1);
            k1 = ((w(2) - w(1))/(h(2) - h(1))^(1/k3(1)))*ones(3,1);
            k2 = [-h(1);-h(1);h(1) + delta];

            bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);

        else

            k3 = 20*ones(3,1);
            k0 = zeros(3,1);
            k1 = (w(1)/h(1)^(1/k3(1)))*ones(3,1);
            k2 = zeros(3,1);

            bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);

        end
    end
    for j = 1:width(by)
        h = [500;1000];
        w = [600;1000];
        delta = 10;

        if by(j) > w(2)
            alpha = deg2rad(0.05);

            k3 = ones(3,1);
            k0 = w(2)*ones(3,1);
            k1 = tan((pi/2) - alpha)*ones(3,1);
            k2 = [-h(2);-h(2);h(2) + delta];

            bzy(j) = ((by(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);

        elseif ((by(j) <= w(2)) && (by(j) > w(1)))

            k3 = 6*ones(3,1);
            k0 = w(1)*ones(3,1);
            k1 = ((w(2) - w(1))/(h(2) - h(1))^(1/k3(1)))*ones(3,1);
            k2 = [-h(1);-h(1);h(1) + delta];

            bzy(j) = ((by(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
        else

            k3 = 20*ones(3,1);
            k0 = zeros(3,1);
            k1 = (w(1)/h(1)^(1/k3(1)))*ones(3,1);
            k2 = zeros(3,1);

            bzy(j) = ((by(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);

        end
    end
end
% % 
% % % POSITION %
% % figure(2)
% % subplot(3,1,1)
% % xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
% % ylabel('$x$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
% %
% % subplot(3,1,2)
% % xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
% % ylabel('$y$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
% %
% % subplot(3,1,3)
% % xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
% % ylabel('$z$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
% %
% VELOCITY %
figure(3)
subplot(3,1,1)
xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
ylabel('$v_x$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])

subplot(3,1,2)
xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
ylabel('$v_y$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])

subplot(3,1,3)
xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
ylabel('$v_z$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])


% PLANAR TRAJECTORY %
figure(5)
xlabel('x (m)', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('z (m)', 'FontSize', 14, 'Interpreter', 'latex')
plot(-bx, bzx, '--', 'Color', 'black', 'LineWidth', 2);
hold on
plot(bx, bzx, '--', 'Color', 'black', 'LineWidth', 2);
hold on

figure(6)
xlabel('y (m)', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('z (m)', 'FontSize', 14, 'Interpreter', 'latex')
plot(by, bzy, '--', 'Color', 'black', 'LineWidth', 2);
hold on
plot(-by, bzy, '--', 'Color', 'black', 'LineWidth', 2);
hold on


% % % FUEL STATS %
% % figure(7)
% % xlabel('Test Case' ,'Interpreter', 'tex')
% % ylabel('$\Delta m$' ,'Interpreter', 'latex')
% % legend({'1: Terrain avoidance guidance', '2: Zhang et al. (2017)', '3: Ebrahimi et al. (2008)'}, 'Location', 'best')
% % 
% % figure(8)
% % xlabel('Algorithm', 'Interpreter', 'latex')
% % ylabel('$\Delta m$', 'Interpreter', 'latex')


% THRUST %

figure(9)
xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
ylabel('$T_x$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])

figure(10)
xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
ylabel('$T_y$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])

figure(11)
xlabel('$t$ (s)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
ylabel('$T_z$ (m)', 'FontSize', 14, 'Interpreter', 'latex', 'Color', [0 0 0])
% % 
% % 
% % 
% % 
