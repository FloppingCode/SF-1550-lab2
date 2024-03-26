%% Interpolation

h = 0.25; 
[t, x, y, ~, ~] = kastbana(h);

%% Linjär interpolation

h = 0.25; 
[t, x, y, ~, ~] = kastbana(h);

[xroot, xmax, ymax, interx, intery] = linear_interpolation(t, x, y);

figure;
plot(x, y, 'g.');
hold on;
plot(x, y, 'k');
plot(interx, intery, 'b*');
title('Linjär interpolation');

fprintf("x: ")
disp(xmax)
fprintf("y: ")
disp(ymax)
%% Kvadratisk interpolation

% [t, x, y, ~, ~] = kastbana(h);
h = 0.25; 
[t, x, y, ~, ~] = kastbana(h);

[xroot, xmax, ymax, interx, intery] = quadratic_interpolation(t, x, y);

figure;
plot(x, y, 'k');
hold on;
plot(interx, intery, 'k', 'LineWidth', 1);
plot(xmax, ymax, 'b.', 'MarkerSize', 10); 
plot(xroot, 0, 'b*', 'MarkerSize', 10); 
title('Kvadratisk interpolation');

fprintf("x: ")
disp(xmax)
fprintf("y: ")
disp(ymax)

%% Funktioner
function [xroot, xmax, ymax, interx, intery] = linear_interpolation(times, x, y)
    xroot = 0;
    interx = [];
    intery = [];
    [ymax, i] = max(y);
    xmax = x(i);
    for t = 1:length(times) - 1
        if y(t+1) < 0 && y(t) > 0
            x_points = [x(t) x(t+1)];
            y_points = [y(t) y(t+1)];
            c = [[1; 1]  x_points(:)]\y_points(:);
            slope = c(2);
            intercept = c(1);
            f = [slope, intercept];
            xroot = roots(f);
        end
    end
    
    interx = [interx, xmax, xroot];
    intery = [intery, ymax, 0];
end

function [xroot, xmax, ymax, interx, intery] = quadratic_interpolation(t, x, y)
    xroot = 0;
    
    n = length(x);
    interx = [];
    intery = [];
    for i = 2:2:n-1
        i1 = i - 1; i2 = i; i3 = i+1;
        c0 = y(i1);
        c1 = (y(i2) - y(i1)) / (x(i2) - x(i1));
        c2 = (((y(i3) - y(i2))/(x(i3) - x(i2))) - ((y(i2) - y(i1))/(x(i2) - x(i1)))) / (x(i3) - x(i1));
        
        fm = @(xm) c2*(xm-x(i1)).*(xm-x(i2)) + c1*(xm-x(i1)) + c0; % Elementvis multiplikation med .*
        x_interp = linspace(x(i1), x(i3), 100);
        y_interp = fm(x_interp);
        
        interx = [interx, x_interp];
        intery = [intery, y_interp];
        
        derivative = @(xm) 2*c2*xm - c2*(x(i1) + x(i2)) + c1;
        pot_max = fzero(derivative, x(i2));
        if pot_max >= x(i1) && pot_max <= x(i3)
            xmax = pot_max;
            ymax = fm(xmax);
        end

        if (y(i1) > 0 && y(i2) < 0) || (y(i1) > 0 && y(i3) < 0) || (y(i2) > 0 && y(i3) < 0)
            xroot = fzero(fm, x(i2));
        end
    end
end

%% Kastbana

% Er kod här...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,x,y,vx,vy]=kastbana(h)

%KASTBANA(H) beräknar banan för ett kast med en liten boll.
%
%   Dynamiken ges av en ODE som inkluderar effekten av luftmotståndet,
%
%      r'' = -g*ez-sigma*r'*|r'|/m.
%
%   Funktionen beräknar bollens position och hastighet vid
%   tidpunkter separerade med en given steglängd. Bollen kastas från
%   (X,Y)=(0,2) med hastigheten 30 m/s i 45 graders vinkel uppåt.
%
%   Syntax:
%
%   [T,X,Y,VX,VY] = KASTBANA(H)
%
%   H       - Steglängden mellan tidpunkterna.
%   T       - Vektor med tidpunkter där bollens position och hastighet beräknats.
%   X, Y    - Vektorer med bollens x- och y-koordinater vid tidpunkterna.
%   VX, VY  - Vektorer med bollens hastigheter i x- och y-led vid tidpunkterna.

m = 56e-3;     % Massan (kg) = 56 gram
ra = 6.6e-2/2; % 6.6 cm in diameter

g=9.81;      % Tyngdaccelerationen (m/s^2)

rho=1.2;     % Luftens densitet (kg/m^3)
A=ra^2*pi;   % Kroppens tvärsnittsarea (m^2)
Cd=0.47;     % Luftmotståndskoefficient,
% "drag coefficient" (dimensionslös)
% Läs mer på http://en.wikipedia.org/wiki/Drag_coefficient

sigma = rho*A*Cd/2; % Totala luftmotståndet

T  = 5;      % Sluttid
v0 = 30;     % Utkasthastighet
al = pi/4;   % Utkastvinkel

% Begynnelsevärden

r0 = [0 2]';                   % Position
r1 = [v0*cos(al) v0*sin(al)]'; % Hastighet

% ODEns högerled

f = @(u) [u(3:4); -u(3:4)*norm(u(3:4),2)*sigma/m - [0;g]];  % RHS

u = [r0;r1];
U = u';
t = 0:h:T;

% Runge-Kutta 4

for tn=t(1:end-1)
    s1 = f(u);
    s2 = f(u + h/2*s1);
    s3 = f(u + h/2*s2);
    s4 = f(u + h*s3);
    u = u + h/6*(s1 + 2*s2 + 2*s3 + s4);
    U = [U; u'];
end

x  = U(:,1);
y  = U(:,2);
vx = U(:,3);
vy = U(:,4);

end

