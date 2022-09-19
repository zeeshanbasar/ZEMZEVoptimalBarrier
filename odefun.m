function dx = odefun(t,x)
dx = zeros(7,1);
    
    %%% PLANET PARAMS %%%
    g = [0 0 -3.7114]';
    ge = 9.807;
    
    %%% LANDER TC %%%
    rd = [0 0 0]';
    vd = [0 0 0]';
    tf = 100;
    tgo = tf - t;
    
    %%% THRUST PARAMS %%%
    Isp = 226;
    
    %%% BARRIER AND AUG TERMS CONSTS %%%
    k1 = [0.17;0.17;0.17]*1000;
    k2 = [0;0;1];
    k3 = [6;6;6];
    
    l1 = [1;1;1];
    l2 = [6;6;6];
    
    %%% BARRIER DEFN %%%    
    if x(1) >=0 
        rho1 = k1(1)*(x(3) + k2(1))^(1/k3(1));
    else
        rho1 = -k1(1)*(x(3) + k2(1))^(1/k3(1));
    end
    
    if x(2) >=0 
        rho2 = k1(2)*(x(3) + k2(2))^(1/k3(2));
    else
        rho2 = -k1(2)*(x(3) + k2(2))^(1/k3(2));
    end
    
    rho3 = k2(3);
    
    %%% DIST FROM BARRIER %%%    
    d1 = abs(x(1) - rho1);
    d2 = abs(x(2) - rho2);
    d3 = abs(x(3) - rho3);
    
    %%% AUGMENTATION %%%
    phi1 = l2(1)/(d1^2 + l1(1));
    phi2 = l2(2)/(d2^2 + l1(2));
    phi3 = l2(3)/(d3^2 + l1(3));
    
    b1 = exp(-phi1);
    b2 = exp(-phi2);
    b3 = exp(-phi3);
    

    
    %%% ZEM/ZEV %%%
    r = [x(1) x(2) x(3)]';
    v = [x(4) x(5) x(6)]';
    
    ZEM = rd - (r + v*tgo + 0.5*g*tgo^2);
    ZEV = vd - (v + g*tgo);
    
    %%% THRUST GENERATION %%%

    %%% NEW OGL %%%
%     p1 = (-2*d1*l2(1)*b1)/(d1^2 + l1(1))^2;
%     p2 = (-2*d2*l2(2)*b2)/(d2^2 + l1(2))^2;
%     p3 = (-2*d3*l2(3)*b3)/(d3^2 + l1(3))^2;
%     
%     p = [p1; p2; p3];
%     a = ((6*ZEM/tgo^2) - (2*ZEV/tgo) - (p/18)*(tgo^2));
%     
%     T = a*x(7);
    
    
    %%% ORIGINAL OGL %%%
    A = [0 0 1]';
    del = 1;
    phi = del^2/3;
    c = 500;
    a_av = A*c*(r(3)^2 - phi)*(tgo^2)/(24*(r(3)^2 + phi)^2);
    a = ((6*ZEM/tgo^2) - (2*ZEV/tgo) + a_av);
    
    T = a*x(7);
    
    %%% ATM PERTURB %%%
%     ap = 0.2*(T./x(7))*sin(pi*t/3);
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

