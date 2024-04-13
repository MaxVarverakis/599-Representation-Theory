clear all; close all; clc

syms a b c d

X = [ a b ; c d ];

J = @(s) [ 2*s -1 ; -1 2*s ];
chk = @(s) [ conj(s.a) conj(s.c) ; conj(s.b) conj(s.d) ]*[ s.a s.b ; s.c s.d ];
%%

% clear sol;
clear A b1
clc

sol = solve(X'*X == J(sqrt(3)/2))
% sol = solve(X'*X == J(-1/2))
% sol = vpasolve(X'*X == J(-1/2),'Random',true)

% latex(sym(J(sqrt(3)/2)))

c = chk(sol)

% t = [-1.25 .4 ; .4 -1.25] % for J(sqrt(3)/2), but not quite exact answer for that s
% b = [ -7/6 3/7 ; 3/7 -7/6 ]

A = solve(a^2/4 + 1/a^2 == sqrt(3))

b1 = @(k) [ -A(k)/2 1/A(k) ; 1/A(k) -A(k)/2 ]
%% 

% test = sym([ exp(1i*-pi/4)/2 exp(1i*pi/4) ; exp(1i*pi/4) exp(1i*-pi/4)/2 ])
% test = @(a) sym([ (-(2^(1/2)*a)/2)/(2*(2^(1/2)*a)/2) ((2^(1/2)*a)/2)/((2^(1/2)*a)/2) ; ((2^(1/2)*a)/2)/((2^(1/2)*a)/2) (-(2^(1/2)*a)/2)/(2*(2^(1/2)*a)/2) ])
% test(1)'*test(1)

%%

clc

for ii = 1:4
    latex(b1(ii))
    % eval(b1(ii))
    latex(sym(simplify(b1(ii)'*b1(ii))))
end

%%

clear P p1 p2; clc

s = exp(1i*5*pi/6);

P = [ s^1 0 ; 0 s^2 ];
p1 = [ -s^2 0 ; 1 1 ];
p2 = [ 1 s^2 ; 0 -s^2 ];


% us1 = simplify((b1(1)*P) * p1 * inv(b1(1)*P))
% us2 = simplify((b1(1)*P) * p2 * inv(b1(1)*P))

us1 = eval((b1(1)*P) * p1 * (b1(1)*P)^-1);
us2 = eval((b1(1)*P) * p2' * (b1(1)*P)^-1);

us1^-1
us1'

% u1 = 1/2 * exp(-1i*pi/6) * [ sqrt(3)*exp(i*atan(1/sqrt(2))) 1 ; 1 -sqrt(3)*exp(-i*atan(1/sqrt(2))) ]

% diagonal entries have magnitude cos(pi/6) = sqrt(3)/2
% solve(a^2 + us1(1,2)^2 == conj(a)) % weird property from mat. mult us1*us1
% diagonals have 1/2^(1/2) +/- 1/2 and angle -(2*pi)/3 or pi/3 respectively for +/-
% solve(3/4 * (exp(i*a) + exp(i*b)) == trace(us1), a+b == -2*pi/3)
% from trace and det
% sum of diagonal angles is -2*pi/3
% diag angles: -pi/3+asec(sqrt(3)) and -pi/3 - asec(sqrt(3)) % NO

% d1 = sqrt(3)/2*exp(-1i*(pi/3-asec(sqrt(3))));
% d2 = sqrt(3)/2*exp(-1i*(pi/3+asec(sqrt(3))));

d1 = 1/4*(i+sqrt(2))*(-i+sqrt(3));
d2 = -1/4*(-i+sqrt(2))*(-i+sqrt(3));
% off-diagonal: 1/2 * exp(-i*pi/6)
od = 1/2 * exp(1i*pi/6);

% [ d1 od ; od d2 ] % inverse analytical approximation

% us2^-1
% us2'
% 
% [ d2 od ; od d1 ]

% u2 = 1/2 * exp(-1i*pi/6) * [ -sqrt(3)*exp(-i*atan(1/sqrt(2))) 1 ; 1 sqrt(3)*exp(i*atan(1/sqrt(2))) ]

% us1*us2*us1
% us2*us1*us2

% us1*us2 = (1/2 * exp(-1i*pi/6))^2 * [ -2 2*sqrt(3)*exp(i*atan(1/sqrt(2))) ; -2*sqrt(3)*exp(-i*atan(1/sqrt(2))) -2 ]
