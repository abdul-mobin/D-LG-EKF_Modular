function [g, P] = update_lgkf(g, P, z, az_std)
    R = diag([max(az_std^2, 1e-8), 10.0^2, 1.0^2]);
    
    % g: Predicted state (5x1), \mu_{k|k-1} of the paper
    % P: Predicted covariance (5x5), \P_{k|k-1} of the paper
    % Equation implementation sequence: 26 -> 25 -> 24 -> 23 -> 21 -> 22  
    
    v = lie2vec_radar(g);
    px=v(1); py=v(2); pth=v(3); pv=v(4); 
    rho_sq = px^2+py^2; alpha = atan2(py,px);
    
    H = zeros(3,5);
    H(1,1) = (px*sin(pth) - py*cos(pth)) / max(rho_sq, 1e-6); 
    H(1,2) = (px*cos(pth) + py*sin(pth)) / max(rho_sq, 1e-6);
    H(2,1) = (px*cos(pth) + py*sin(pth)) / max(sqrt(rho_sq), 1e-6); 
    H(2,2) = (-px*sin(pth) + py*cos(pth)) / max(sqrt(rho_sq), 1e-6);
    H(3,1) = -pv*sin(alpha-pth)*H(1,1); H(3,2) = -pv*sin(alpha-pth)*H(1,2);
    H(3,3) = pv*sin(alpha-pth); H(3,4) = cos(alpha-pth);
    
    S = H * P * H' + R; K = P * H' / S;
    innov = [z(1)-alpha; z(2)-sqrt(rho_sq); z(3)-pv*cos(alpha-pth)];
    innov(1) = atan2(sin(innov(1)), cos(innov(1)));
    
    m = K * innov;
    g = g * expm(hat_map(m));   
    P = phi(m)*(eye(5) - K * H) * P*phi(m)';
end