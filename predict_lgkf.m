function [g, P] = predict_lgkf(g, P, Q, dt)
    v = lie2vec_radar(g);
    Omega = [v(4)*dt; 0; v(5)*dt; 0; 0];
    g_next = g * expm(hat_map(Omega));
    F = Ad_G_radar(expm(-hat_map(Omega))); 
    G_mat = zeros(5,5); G_mat(1,4) = dt; G_mat(3,5) = dt;
    F = F + G_mat;
    g = g_next; P = F * P * F' + Q;
end

