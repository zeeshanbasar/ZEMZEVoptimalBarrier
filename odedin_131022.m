function dx = odefin(t,x,tf,bf,cl)
dx = zeros(7,1);
    
    %%% PLANET PARAMS %%%
    g = [0 0 -3.7114]';
    ge = 9.807;
    
    %%% LANDER TC %%%
    rd = [0 0 0]';
    vd = [0 0 0]';
%     tf = 200;
    tgo = tf - t;
    
    %%% THRUST PARAMS %%%
    Isp = 225;
    
    %%% ZEM/ZEV %%%
    r = [x(1) x(2) x(3)]';
    v = [x(4) x(5) x(6)]';
    
    ZEM = rd - (r + v*tgo + 0.5*g*tgo^2);
    ZEV = vd - (v + g*tgo);
    
    %%% ACCELERATION COMMAND %%%
    
    if cl == 3
        a = ((6*ZEM/tgo^2) - (2*ZEV/tgo));
        
    elseif cl == 2
        A = [0 0 1]';
        del = 1;
        phi = del^2/3;
        c = 500;
        a_av = A*c*(r(3)^2 - phi)*(tgo^2)/(24*(r(3)^2 + phi)^2);
        a = ((6*ZEM/tgo^2) - (2*ZEV/tgo) + a_av);
        
    elseif cl == 1
        %%% BARRIER AND AUG TERMS CONSTS %%%

        if bf == 1
            k0 = [0;0;0];
            k1 = [0.17;0.17;0.17]*1000;
            k2 = [0;0;1];
            k3 = [6;6;6];
            rho3 = k2(3);

        elseif bf == 2
            k0 = [0;0;0];
            k1 = [1.5;1.5;1.5];
            k2 = [0;0;1];
            k3 = [1;1;1];
            rho3 = k2(3);

        elseif bf == 3
            if x(3) > 200
                k0 = [300;300;300];
                k1 = [178;178;178];
                k2 = [-200;-200;1];
                k3 = [6;6;6];
            else
                k0 = [0;0;0];
                k1 = [231;231;231];
                k2 = [0;0;1];
                k3 = [20;20;20];
            end
            rho3 = k2(3);
            
        elseif bf == 4
            if (x(3) > 200)
                k0 = 0*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-200;-200;201];
                k3 = 1*ones(3,1);
            else
                k0 = [0;0;0];
                k1 = [230;230;230];
                k2 = [0;0;1];
                k3 = [20;20;20];            
            end

            if (x(3) > 200) && (norm([x(1) x(2)],"inf") > 300)
                rho3 = 210;
            elseif (x(3) > 200) && (norm([x(1) x(2)],"inf") < 300)
                rho3 = 10;
            else
                rho3 = 1;
            end        

        elseif bf == 5
            if x(3) > 1000
                k0 = -1000*ones(3,1);
                k1 = 300*ones(3,1);
                k2 = [-1000;-1000;1001];
                k3 = ones(3,1);
            elseif ((x(3) <= 1000) && (x(3) > 500))
                k0 = 600*ones(3,1);
                k1 = 142*ones(3,1);
                k2 = [-500;-500;501];
                k3 = 6*ones(3,1);
            else
                k0 = zeros(3,1);
                k1 = 441*ones(3,1);
                k2 = zeros(3,1);
                k3 = 20*ones(3,1);
            end
 %%% fix this %%%
 %%% this apparently does not need any fix now, idek just keep it like this


            if (x(3) > 1000) && (norm([x(1) x(2)],"inf") >= 1000)
                rho3 = 1010;
            elseif ((x(3)<=1000) && (x(3) >= 500)) && ((norm([x(1) x(2)],"inf") >= 600) && (norm([x(1) x(2)],"inf") <= 1000))
                rho3 = 510;
%             elseif ((x(3) <= 1000)) 
%                 rho3 = 10;
            else
                rho3 = 10;
            end        

% %             if (x(3) > 1500) && (min(abs([x(1),x(2)])) >= 1000)
% %                 rho3 = 1510;
% %             elseif ((x(3)<=1500) && (x(3) >= 1000)) && ((min(abs([x(1),x(2)])) >= 600) && (min(abs([x(1),x(2)])) <= 1000))
% %                 rho3 = 1010;
% % %             elseif ((x(3) <= 1000)) 
% % %                 rho3 = 10;
% %             else
% %                 rho3 = 10;
% %             end    
% %         end
        end

        l1 = [1;1;1];
        l2 = [2;2;2];
        l4 = [300;300;300];

        %%% BARRIER DEFN %%%    
        if x(1) >=0 
            rho1 = ((k1(1)*((x(3) + k2(1))^(1/k3(1)))) + k0(1));
        else
            rho1 = -((k1(1)*((x(3) + k2(1))^(1/k3(1)))) + k0(1));
        end

        if x(2) >=0 
            rho2 = ((k1(2)*((x(3) + k2(2))^(1/k3(2)))) + k0(2));
        else
            rho2 = -((k1(2)*((x(3) + k2(2))^(1/k3(2)))) + k0(2));
        end


        %%% DIST FROM BARRIER %%%    
        s1 = x(1) - rho1;
        s2 = x(2) - rho2;
        s3 = x(3) - rho3;

        d1 = abs(s1);
        d2 = abs(s2);
        d3 = abs(s3);

        %%% AUGMENTATION %%%
        phi1 = l2(1)/(d1^2 + l1(1));
        phi2 = l2(2)/(d2^2 + l1(2));
        phi3 = l2(3)/(d3^2 + l1(3));

        b1 = exp(-phi1);
        b2 = exp(-phi2);
        b3 = exp(-phi3);



        %%% NEW OGL %%%
        p1 = (-2*sign(s1)*d1*l2(1)*l4(1)*b1)/(d1^2 + l1(1))^2;
        p2 = (-2*sign(s2)*d2*l2(2)*l4(2)*b2)/(d2^2 + l1(2))^2;
        p3 = (-2*sign(s3)*d3*l2(3)*l4(3)*b3)/(d3^2 + l1(3))^2;

        p = [p1; p2; p3];
        a = ((6*ZEM/tgo^2) - (2*ZEV/tgo) - (p/18)*(tgo^2)); 
        
    end
    


    %%% THRUST GENERATION %%%    
    T = a*x(7);
    
    % uncomment for thrust constraints
% %     Tmax = 18600;
% %     for i=1:3
% %         T(i) = min(max(T(i), -.8*Tmax), .8*Tmax);
% %     end

    
    %%% ATM PERTURB %%%
%     ap = 0.5*(T./x(7))*sin(pi*t/3);
    ap = [0;0;0];
    

    %%% ODE EQNS %%%
    dx(1) = x(4);%x
    dx(2) = x(5);%y
    dx(3) = x(6);%z
    dx(4) = g(1)  + T(1)/x(7) + ap(1);%vx
    dx(5) = g(2)  + T(2)/x(7) + ap(2);%vy
    dx(6) = g(3)  + T(3)/x(7) + ap(3);%vz
    dx(7) = -norm(T)/(Isp*norm(ge));%m
end

