function [Z, estados, t] = shooting_function(parametros, s_0)
    alfa = parametros(1);
    gamma_1 = parametros(2);
    gamma_2 = parametros(3);
    gamma_3 = parametros(4);
    p_m0 = parametros(5);
    tf = parametros(6);

    %Fmax = 0.194226;
    Fmax = 0.826103977267828;
    %Isp = 8.15499E-3;
    Isp = 0.007983105820217;
    %g0 = 3806.66;
    g0 = 3.598141767655429e+03;
    %m_s= 3.2891E5;
    m_s = 328900;
    %ro_s = 388.97;
    ro_s = 3.887800000000000e+02;
    %mu = 0.012150664267;
    mu = 0.012150664267000;
    w_0 = 0;
   % w_s = -0.92518;
    w_s = -0.925179999954674;

    m_0 = 1;
    R_e = 1.658E-2;
    d_v_0 = 3.021867; %3.032757; %2.5-3.1
    H_0 = 4.3418E-4;
    R_f = 0.1301;

    r_0 = R_e + H_0;
    [x_0, y_0, v_x_0, v_y_0 , p_x_0, p_y_0, p_vx_0, p_vy_0] = calcular_ci(alfa, mu, s_0, r_0, d_v_0, gamma_1, gamma_2, gamma_3);
    ci = [x_0 y_0 v_x_0 v_y_0 m_0 p_x_0 p_y_0 p_vx_0 p_vy_0 p_m0];
    
    options=odeset('RelTol',1e-8,'AbsTol',1e-8);
    tspan=0:.01:tf;
    parametros = [Isp g0 Fmax ro_s mu w_0 w_s m_s];
    [t,estados]=ode45(@(t,st)motor_bicircular(t,st,parametros),tspan,ci,options);

    v_f = sqrt(mu/R_f);
    w = w_0 + w_s*tf;
    parametros = [mu R_f Fmax Isp g0 v_f ro_s w m_s];
    Z = shooting_array(parametros, estados);
    
end

