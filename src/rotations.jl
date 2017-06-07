module Rotations

Rx(θ::Real) = begin
    s = sin(θ)
    c = cos(θ)
    [1 0 0;
     0 c -s;
     0 s c]
end

Ry(θ::Real) = begin
    s = sin(θ)
    c = cos(θ)
    [c 0 s;
     0 1 0;
     -s 0 c]
end

Rz(θ::Real) = begin
    s = sin(θ)
    c = cos(θ)
    [c -s 0;
     s c 0;
     0 0 1]
end

end
