function X_hat = hat_map(v)
    X_hat = zeros(6);
    X_hat(1:2, 1:2) = [0, -v(3); v(3), 0]; 
    X_hat(1:2, 3) = [v(1); v(2)]; 
    X_hat(4:5, 6) = [v(4); v(5)]; 
end

