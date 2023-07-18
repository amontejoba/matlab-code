function [estados_d] = motor_bicircular(t,estados, parametros)
    x = estados(1);
    y = estados(2);
    v_x = estados(3);
    v_y = estados(4);
    m = estados(5);
    p_x = estados(6);
    p_y = estados(7);
    p_vx = estados(8);
    p_vy = estados(9);
    p_m = estados(10);
    


    Isp = parametros(1);
    g0 = parametros(2);
    Fmax = parametros(3);
    ro_s = parametros(4);
    mu = parametros(5);
    w_0 = parametros(6);
    w_s = parametros(7);
    m_s = parametros(8);

    w = w_0 + w_s*t;
   
    omega_s_x = x - (mu*(2*mu+2*x-2)/(2*((mu+x-1)^2+y^2)^(3/2))) - (m_s*(2*x-2*ro_s*cos(w)))/(2*((x-ro_s*cos(w))^2+(y-ro_s*sin(w))^2)^(3/2)) + ((2*mu+2*x)*(mu-1))/(2*((mu+x)^2+y^2)^(3/2)) - m_s*cos(w)/ro_s^2;
    omega_s_y = y - (m_s*(2*y-2*ro_s*sin(w)))/(2*((x-ro_s*cos(w))^2+(y-ro_s*sin(w))^2)^(3/2)) - (mu*y)/(((mu+x-1)^2+y^2)^(3/2)) - (m_s*sin(w))/(ro_s^2) + (y*(mu-1))/(((mu+x)^2+y^2)^(3/2));

    sigma_1xx = (x-ro_s*cos(w))^2+(y-ro_s*sin(w))^2;
    sigma_2xx = (mu+x-1)^2+y^2;
    sigma_3xx = (mu+x)^2+y^2;

    omega_s_xx = (mu-1)/(sigma_3xx^(3/2)) - m_s/(sigma_1xx^(3/2)) - mu/(sigma_2xx^(3/2)) + (3*mu*(2*mu+2*x-2)^2)/(4*sigma_2xx^(5/2)) + (3*m_s*(2*x-2*ro_s*cos(w))^2)/(4*sigma_1xx^(5/2)) - (3*((2*mu+2*x)^2)*(mu-1))/(4*sigma_3xx^(5/2))+1;
    omega_s_yy = (mu-1)/(sigma_3xx^(3/2)) - m_s/(sigma_1xx^(3/2)) - mu/(sigma_2xx^(3/2)) - (3*(y^2)*(mu-1))/(sigma_3xx^(5/2)) + (3*m_s*(2*y-2*ro_s*sin(w))^2)/(4*(sigma_1xx^(5/2))) + (3*mu*(y^2))/(sigma_2xx^(5/2))+1;

    t1 = 2*x-ro_s*cos(w);
    t2 = 2*y-ro_s*sin(w); 
    t3 = (x-ro_s*cos(w))^2;
    t4 = (y-ro_s*sin(w))^2;
    t5 = 2*mu+2*x-2;
    t6 = (mu+x-1)^2;
    t7 = (t6+y^2)^(5/2);
    t8 = (2*mu+2*x);
    t9 = (mu-1);
    t10= (mu+x)^2;
    t11= (t10+y^2)^(5/2);

    t2 = 1. - mu;
    t3 = x + mu;
    t4 = t3 ^ 2;
    t5 = y ^ 2;
    t6 = t4 + t5;
    t7 = sqrt(t6);
    t13 = x + mu - 1.;
    t14 = t13 ^ 2;
    t15 = t14 + t5;
    t16 = sqrt(t15);
    t22 = cos(w);
    t24 = x - ro_s * t22;
    t25 = t24 ^ 2;
    t26 = sin(w);
    t28 = y - ro_s * t26;
    t29 = t28 ^ 2;
    t30 = t25 + t29;
    t31 = sqrt(t30);
    t66 = t6 ^ 2;
    t69 = t2 / t7 / t66;
    t73 = t15 ^ 2;
    t76 = mu / t16 / t73;
    t80 = t30 ^ 2;
    t83 = m_s / t31 / t80;
    t98 = 3. * t69 * y * t3 + 3. * t76 * y * t13 + 3. * t83 * t28 * t24;

    %omega_s_xy = (3*m_s*t1*t2)/(4*(t3+t4)^(5/2)) + (3*mu*y*t5)/(2*t7) - (3*y*t8*t9)/(2*t11)
    %omega_s_yx = omega_s_xy;
    omega_s_xy = t98;
    omega_s_yx = omega_s_xy;
    p_v_mod = sqrt(p_vx^2+p_vy^2);
    S = -p_m-(Isp*g0*p_v_mod)/m;
    N = 0;
    if S>0 
        N = 0;
    elseif S<=0
        N = 1;
    end

    u_1 = -N*(p_vx/p_v_mod);
    u_2 = -N*(p_vy/p_v_mod);

    u_mod = sqrt(u_1^2+u_2^2);

    estados_d(1) = v_x;
    estados_d(2) = v_y;
    estados_d(3) = 2*v_y + omega_s_x + u_1*Fmax/m;
    estados_d(4) = -2*v_x + omega_s_y + u_2*Fmax/m;
    estados_d(5) = -u_mod*Fmax/(Isp*g0);
    estados_d(6) = -p_vx*omega_s_xx - p_vy*omega_s_yx;
    estados_d(7) = -p_vx*omega_s_xy - p_vy*omega_s_yy;
    estados_d(8) = -p_x + 2*p_vy;
    estados_d(9) = -p_y - 2*p_vx;
    estados_d(10)= (Fmax*(p_vx*u_1 + p_vy*u_2))/m^2;
    estados_d = estados_d';
end