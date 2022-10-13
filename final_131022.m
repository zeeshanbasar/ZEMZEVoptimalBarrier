clc;clear;close all

%% INIT %%
IC = load('init-meh.mat').IC;
% IC = [351.965576688550,532.722905849107,2108.39674231030,25.5906010528496,27.1238598366170,-85.7789813376845,2094.60970946621];
t0 = 0;
tf = 75;

%% DE solve %%
% run < parpool(8) > in cmd window before using < parfor >
% run < delete(gcp('nocreate')) > to close the parallel pool, must be done
% to ensure RAM space is released

warning('off','all')

N = tf*1000;
tSpan = linspace(t0,tf,N);
options = odeset('Events',@terminate);

del_m = zeros(length(IC),1);

bf = 5;
% 1: smooth
% 2: glideslope
% 3: 1-step
% 4: flat top, 1-step
% 5: flat top, 2-step

tic
for cl = 1
% 3: classical ZEM-ZEV
% 2: self adjusting ZEM-ZEV
% 1: new ogl
for i = 1:length(IC) % set this to < parfor > when needed
    
    init = IC(i,:);

    [t,x] = ode45(@(t,x) odefin(t,x,tf,bf,cl), tSpan, init);

    %% PLOTS %%
    N = length(x);

    % TRAJECTORY %
    figure(1)
    plot3(x(:,1),x(:,2),x(:,3), 'LineWidth', 1.5);
    hold on
    plot3(x(1,1),x(1,2),x(1,3),'ks')
    hold on
    plot3(x(N,1),x(N,2),x(N,3),'ro')
    hold on
    grid on
% % 
% %     % POSITION %
% %     figure(2)
% %     subplot(3,1,1)
% %     plot(t,x(:,1))
% %     hold on
% %     grid on
% % 
% %     subplot(3,1,2)
% %     plot(t,x(:,2))
% %     hold on
% %     grid on
% % 
% %     subplot(3,1,3)
% %     plot(t,x(:,3))
% %     hold on
% %     grid on
% % 
% %     % VELOCITY %
% %     figure(3)
% %     subplot(3,1,1)
% %     plot(t,x(:,4))
% %     hold on
% %     grid on
% % 
% %     subplot(3,1,2)
% %     plot(t,x(:,5))
% %     hold on
% %     grid on
% % 
% %     subplot(3,1,3)
% %     plot(t,x(:,6))
% %     hold on
% %     grid on
% % 
% % % % 
% % % %     % % THRUST %
% % % %     % A = [0 0 1]';
% % % %     % del = 1;
% % % %     % phi = del^2/3;
% % % %     % c = 500;
% % % %     % tgo = tf - t;
% % % %     % r = [x(:,1) x(:,2) x(:,3)];
% % % %     % v = [x(:,4) x(:,5) x(:,6)];
% % % %     % ZEM = rd' - (r + v.*tgo + 0.5*g'.*tgo.^2);
% % % %     % ZEV = vd' - (v + g'.*tgo);
% % % %     % a_av = c*(x(:,3).^2 - phi).*(tgo.^2)./(24*(x(:,3).^2 + phi).^2);
% % % %     % T = ((6*ZEM./tgo.^2) - (2*ZEV./tgo) + a_av).*x(:,7);
% % % %     % 
% % % %     % figure(5)
% % % %     % subplot(3,1,1)
% % % %     % plot(t,T(:,1))
% % % %     % grid on
% % % %     % 
% % % %     % subplot(3,1,2)
% % % %     % plot(t,T(:,2))
% % % %     % grid on
% % % %     % 
% % % %     % subplot(3,1,3)
% % % %     % plot(t,T(:,3))
% % % %     % grid on
% % % % 
% % % % 
% %     
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
            if bx(j) > 1000
                k0 = -1000*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-1000;-1000;1001];
                k3 = ones(3,1);
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            elseif ((bx(j) <= 1000) && (bx(j) > 600))
                k0 = 600*ones(3,1);
                k1 = 142*ones(3,1);
                k2 = [-500;-500;501];
                k3 = 6*ones(3,1);
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            else
                k0 = zeros(3,1);
                k1 = 441*ones(3,1);
                k2 = zeros(3,1);
                k3 = 20*ones(3,1);
                bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            end
% %             if bx(j) > 1000
% %                 k0 = -1000*ones(3,1);
% %                 k1 = 300*ones(3,1);
% %                 k2 = [-1500;-1500;1501];
% %                 k3 = ones(3,1);
% %                 bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
% %             elseif ((bx(j) <= 1000) && (bx(j) > 600))
% %                 k0 = 600*ones(3,1);
% %                 k1 = 142*ones(3,1);
% %                 k2 = [-1000;-1000;1001];
% %                 k3 = 6*ones(3,1);
% %                 bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
% %             else
% %                 k0 = zeros(3,1);
% %                 k1 = 426*ones(3,1);
% %                 k2 = zeros(3,1);
% %                 k3 = 20*ones(3,1);
% %                 bzx(j) = ((bx(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
% %             end
        end
        
        for j = 1:width(by)
            if by(j) > 1000
                k0 = -1000*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-1000;-1000;1001];
                k3 = ones(3,1);
                bzy(j) = ((by(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            elseif ((by(j) <= 1000) && (by(j) > 600))
                k0 = 600*ones(3,1);
                k1 = 142*ones(3,1);
                k2 = [-500;-500;501];
                k3 = 6*ones(3,1);
                bzy(j) = ((by(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            else
                k0 = zeros(3,1);
                k1 = 441*ones(3,1);
                k2 = zeros(3,1);
                k3 = 20*ones(3,1);
                bzy(j) = ((by(j)-k0(1))/k1(1)).^(k3(1)) - k2(1);
            end
% %             if by(j) > 1000
% %                 k0 = -1000*ones(3,1);
% %                 k1 = 300*ones(3,1);
% %                 k2 = [-1500;-1500;1501];
% %                 k3 = ones(3,1);
% %                 bzy(j) = ((by(j)-k0(2))/k1(2)).^(k3(2)) - k2(2);
% %             elseif ((by(j) <= 1000) && (by(j) > 600))
% %                 k0 = 600*ones(3,1);
% %                 k1 = 142*ones(3,1);
% %                 k2 = [-1000;-1000;1001];
% %                 k3 = 6*ones(3,1);
% %                 bzy(j) = ((by(j)-k0(2))/k1(2)).^(k3(2)) - k2(2);
% %             else
% %                 k0 = zeros(3,1);
% %                 k1 = 426*ones(3,1);
% %                 k2 = zeros(3,1);
% %                 k3 = 20*ones(3,1);
% %                 bzy(j) = ((by(j)-k0(2))/k1(2)).^(k3(2)) - k2(2);
% %             end
        end 
    end

    figure(5)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('z (m)', 'FontSize', 14)

    plot(x(:,1),x(:,3), 'LineWidth', 1, 'DisplayName',num2str(i))
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

    plot(x(:,2),x(:,3), 'LineWidth', 1, 'DisplayName',num2str(i))
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

    del_m(i,cl) = x(1,7) - x(N-1,7);

end
figure(7)
xlabel('Test Case #')
ylabel('\Deltam')

scatter(1:length(IC),del_m(:,cl), 'LineWidth', 1);
hold on
grid on

mu(cl) = mean(del_m(:,cl));
disp(['mu = ', num2str(mu)]);

sd(cl) = std(del_m(:,cl));
disp(['sd = ', num2str(sd)]);

end
toc

figure(8)
xlabel('Algorithm')
ylabel('\Deltam stats')

boxplot(del_m(:,1:cl));
grid on

% % [h,p,ci,stats] = ttest(del_m(:,1)-mu(1))
% % 
% % [h2,p2,ci2,stats2] = ttest(del_m(:,2)-mu(2))






