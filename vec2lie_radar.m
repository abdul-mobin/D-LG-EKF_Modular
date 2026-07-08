function G = vec2lie_radar(vec)
    G = eye(6);
    G(1:2, 1:2) = [cos(vec(3)), -sin(vec(3)); sin(vec(3)), cos(vec(3))]; 
    G(1:2, 3) = [vec(1); vec(2)]; 
    G(4:5, 6) = [vec(4); vec(5)]; 
end

