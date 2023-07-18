function [estados_d] = motor_bicircular_sin_coestados(t,estados, parametros)
    x = estados(1);
    y = estados(2);
    v_x = estados(3);
    v_y = estados(4);
    m = estados(5);

    Isp = parametros(1);
    g0 = parametros(2);
    Fmax = parametros(3);
    ro_s = parametros(4);
    mu = parametros(5);
    w_0 = parametros(6);
    w_s = parametros(7);
    m_s = parametros(8);
    u_1 = parametros(9);
    u_2 = parametros(10);

    w = w_0 + w_s*t;
    
    omega_s_x = x - (mu*(2*mu+2*x-2)/(2*((mu+x-1)^2+y^2)^(3/2))) - (m_s*(2*x-2*ro_s*cos(w)))/(2*((x-ro_s*cos(w))^2+(y-ro_s*sin(w))^2)^(3/2)) + ((2*mu+2*x)*(mu-1))/(2*((mu+x)^2+y^2)^(3/2)) - m_s*cos(w)/ro_s^2;
    omega_s_y = y - (m_s*(2*y-2*ro_s*sin(w)))/(2*((x-ro_s*cos(w))^2+(y-ro_s*sin(w))^2)^(3/2)) - (mu*y)/(((mu+x-1)^2+y^2)^(3/2)) - (m_s*sin(w))/(ro_s^2) + (y*(mu-1))/(((mu+x)^2+y^2)^(3/2));

    u_mod = sqrt(u_1^2+u_2^2);

    estados_d(1) = v_x;
    estados_d(2) = v_y;
    estados_d(3) = 2*v_y + omega_s_x + u_1*Fmax/m;
    estados_d(4) = -2*v_x + omega_s_y + u_2*Fmax/m;
    estados_d(5) = -u_mod*Fmax/(Isp*g0);
    estados_d = estados_d';
end

