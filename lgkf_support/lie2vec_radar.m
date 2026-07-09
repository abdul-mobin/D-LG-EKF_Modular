function vec = lie2vec_radar(G)
    vec = [G(1,3); G(2,3); atan2(G(2,1), G(1,1)); G(4,6); G(5,6)];
end

