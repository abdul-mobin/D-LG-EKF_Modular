function true_state = generate_ground_truth(c, n)
    % Generate a ground-truth trajectory for the CTRV motion model.
    % The state vector is [x; y; theta; v; w], where v is speed and w is turn rate.

    % Allocate storage for the full trajectory.
    true_state = zeros(5, c.N);

    % Initialize the vehicle at the origin with forward speed 200 and zero turn rate.
    true_state(:,1) = [0; 0; 0; 200; 0];

    % Propagate the state over time, adding process noise to the speed and turn rate.
    for k = 2:c.N
        v = true_state(4, k-1);
        w = true_state(5, k-1);
        th = true_state(3, k-1);

        % Update position using the CTRV kinematics.
        if abs(w) > 1e-5
            dx = v * (sin(th + w*c.dt) - sin(th))/w;
            dy = v * (cos(th) - cos(th + w*c.dt))/w;
        else
            % Small-turn-rate limit for the singular case.
            dx = v * cos(th) * c.dt;
            dy = v * sin(th) * c.dt;
        end

        % Update the state vector with the new pose and noisy motion parameters.
        true_state(:, k) = [true_state(1,k-1)+dx;
                            true_state(2,k-1)+dy;
                            th + w*c.dt;
                            v + n.std_v*randn;
                            w + n.std_omega*randn];
    end
end
