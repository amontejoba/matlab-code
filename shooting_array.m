function [Z] = shooting_array(parametros, estados)
    mu = parametros(1);
    R_f = parametros(2);
    Fmax = parametros(3);
    I_sp = parametros(4);
    g_0 = parametros(5);
    v_f = parametros(6);
    ro_s = parametros(7);
    w = parametros(8);
    m_s = parametros(9);

    K = 1;

    x = estados(end,1);
    y = estados(end,2);
    vx = estados(end,3);
    vy = estados(end,4);
    m = estados(end,5);
    px = estados(end,6);
    py  = estados(end,7);
    pvx = estados(end,8);
    pvy = estados(end,9);
    pm = estados(end,10);


    omega_s_x = x - (mu*(2*mu+2*x-2)/(2*((mu+x-1)^2+y^2)^(3/2))) - (m_s*(2*x-2*ro_s*cos(w)))/(2*((x-ro_s*cos(w))^2+(y-ro_s*sin(w))^2)^(3/2)) + ((2*mu+2*x)*(mu-1))/(2*((mu+x)^2+y^2)^(3/2)) - m_s*cos(w)/ro_s^2;
    omega_s_y = y - (m_s*(2*y-2*ro_s*sin(w)))/(2*((x-ro_s*cos(w))^2+(y-ro_s*sin(w))^2)^(3/2)) - (mu*y)/(((mu+x-1)^2+y^2)^(3/2)) - (m_s*sin(w))/(ro_s^2) + (y*(mu-1))/(((mu+x)^2+y^2)^(3/2));
    
    p_v_mod = sqrt(pvx^2+pvy^2);
    S = -pm-(I_sp*g_0*p_v_mod)/m;
    N = 0;
    if S>0 
        N = 0;
    elseif S<=0
        N = 1;
    end

    u_1 = -N*(pvx/p_v_mod);
    u_2 = -N*(pvy/p_v_mod);

    u_mod = sqrt(u_1^2+u_2^2);
    

    hf_1 = (x + mu -1)^2 + y^2 - R_f^2;
    hf_2 = (x + mu -1)*(vx - y) + y*(vy + x + mu - 1);
    hf_3 = (vx - y)^2 + (vy + x + mu - 1)^2 - v_f^2;
    tc = px*y - py*(x + mu - 1) + pvx*vy - pvy*vx;
    pm_t0 = pm + K;
    H = px*vx + py*vy + pvx*(2*vy + omega_s_x + (u_1*Fmax)/m) + pvy*(-2*vx + omega_s_y + (u_2*Fmax)/m) + pm*(-u_mod*(Fmax/(I_sp*g_0)));
    
    Z = [hf_1 hf_2 hf_3 tc pm_t0 H]; 
end

