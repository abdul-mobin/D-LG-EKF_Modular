function Phi = phi(a)
%PHI   Left-Jacobian-like map for the radar Lie algebra.
%   Phi(a) = \sum_{m=0}^\infty (-1)^m/(m+1)! * ad_G(a)^m.
%   For the radar model, the only nonzero entries of ad_G(a) lie in the
%   3x3 upper-left block, so phi reduces to a closed-form expression.
%
%   Input:
%       a - 3- or 5-vector of Lie algebra coordinates [x; y; theta; v; omega]
%   Output:
%       Phi - 5x5 matrix

    a = a(:);
    if numel(a) == 3
        a = [a; 0; 0];
    elseif numel(a) ~= 5
        error('phi: input must be a 3- or 5-vector');
    end

    theta = a(3);
    A = zeros(5);
    A(1,2) = -theta;
    A(1,3) = a(2);
    A(2,1) = theta;
    A(2,3) = -a(1);

    if abs(theta) < 1e-8
        c = 1/2 - theta^2/24;
        d = 1/6 - theta^2/120;
    else
        c = (1 - cos(theta)) / theta^2;
        d = (theta - sin(theta)) / theta^3;
    end

    Phi = eye(5) - c * A + d * (A * A);
end
