function Ad = Ad_G_radar(G)
    th = atan2(G(2,1), G(1,1));
    Ad = eye(5);
    Ad(1:2, 1:2) = [cos(th), -sin(th); sin(th), cos(th)]; 
    Ad(1:2, 3) = [G(2,3); -G(1,3)];
end