function [x_0, y_0, v_x_0, v_y_0, p_x_0, p_y_0, p_vx_0, p_vy_0] = calcular_ci(alfa,mu, s_0, r_0, d_v_0, gamma_1, gamma_2, gamma_3)
    v_0 = sqrt((1-mu)/r_0);
    x_0 = -mu-r_0*sin(alfa);
    y_0 = r_0*cos(alfa);
    v_x_0 = (s_0*(v_0+d_v_0)+r_0)*cos(alfa);
    v_y_0 = (s_0*(v_0+d_v_0)+r_0)*sin(alfa);

    a_0 = (s_0*(v_0+d_v_0)+r_0)/sqrt((s_0*(v_0+d_v_0)+r_0)^2+r_0^2);
    b_0 = (-r_0)/sqrt((s_0*(v_0+d_v_0)+r_0)^2+r_0^2);

    w_0_1 = [-sin(alfa) cos(alfa) 0 0];
    w_0_2 = [0 0 cos(alfa) sin(alfa)];
    w_0_3 = [a_0*cos(alfa) a_0*sin(alfa) b_0*sin(alfa) -b_0*cos(alfa)];

    p_x_0 = gamma_1*w_0_1(1) + gamma_2*w_0_2(1) + gamma_3*w_0_3(1);
    p_y_0 = gamma_1*w_0_1(2) + gamma_2*w_0_2(2) + gamma_3*w_0_3(2);
    p_vx_0 = gamma_1*w_0_1(3) + gamma_2*w_0_2(3) + gamma_3*w_0_3(4);
    p_vy_0 = gamma_1*w_0_1(4) + gamma_2*w_0_2(4) + gamma_3*w_0_3(4);
end

