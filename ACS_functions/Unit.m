
% Unit vector function

function unit = Unit(u)
magnitude = sqrt(u(1)^2 + u(2)^2 + u(3)^2);
unit = [u(1)/magnitude u(2)/magnitude u(3)/magnitude].';
end