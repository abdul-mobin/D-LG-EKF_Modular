function true_state = generate_ground_truth(c, n)
    true_state = zeros(5, c.N);
    true_state(:,1) = [0; 0; 0; 200; 0];
    for k = 2:c.N
        v = true_state(4, k-1); w = true_state(5, k-1); th = true_state(3, k-1);
        if abs(w) > 1e-5
            dx = v * (sin(th + w*c.dt) - sin(th))/w;
            dy = v * (cos(th) - cos(th + w*c.dt))/w;
        else
            dx = v * cos(th) * c.dt; dy = v * sin(th) * c.dt;
        end
        true_state(:, k) = [true_state(1,k-1)+dx; true_state(2,k-1)+dy; th+w*c.dt; ...
                            v + n.std_v*randn; w + n.std_omega*randn];
    end
end
