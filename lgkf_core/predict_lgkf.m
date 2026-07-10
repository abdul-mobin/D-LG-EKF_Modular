function [g, P] = predict_lgkf(g, P, Q, dt)
    v = lie2vec_radar(g);
    Omega = [v(4)*dt; 0; v(5)*dt; 0; 0];
    
    % Compute the next state using the exponential map
    g_next = g * expm(hat_map(Omega));
    g = g_next;

    % Covariance update using the adjoint map
    Phi_Omega = phi(Omega);
    G_mat = zeros(5,5); G_mat(1,4) = dt; G_mat(3,5) = dt;
    F = Ad_G_radar(expm(-hat_map(Omega)))+ Phi_Omega * G_mat;
    P = F * P * F' + Phi_Omega * Q * Phi_Omega';
end

