%% 4a Trapetsregeln och Simpson

ns = [30, 60, 200];

for n = ns
    fprintf('Beräkning för n = %d:\n', n);
    
    V_trap = trapets(n);
    fprintf('Volym med trapetsmetoden för n=%d: %.24f\n', n, V_trap);
    
    V_simpson = simpson(n);
    fprintf('Volym med simpson för n=%d: %.24f\n', n, V_simpson);
end

%% 4b Konvergensstudier

% 4b.1
R = 3;
ref = simpson(1000);
err_trap = [];
err_simp = [];

x_max = 1000;
x_values = 2:2:x_max;

for n = x_values
    err_trap = [err_trap abs(ref - trapets(n))];
    err_simp = [err_simp abs(ref - simpson(n))];
end

loglog(R*x_values.^-1,err_trap, 'DisplayName', 'Trapets');
hold on;
loglog(R*x_values.^-1,err_simp, 'DisplayName', 'Simpson');
grid on

x_trap_fit = linspace(R/x_values(1),R/x_values(end),100);
y_trap_fit = x_trap_fit.^2 / x_trap_fit(1) * err_trap(1);
loglog(x_trap_fit, y_trap_fit, '--', 'DisplayName', 'p = 2');

y_simp_fit = x_trap_fit.^4 / x_trap_fit(1) * err_simp(1);
loglog(x_trap_fit, y_simp_fit, '--', 'DisplayName', 'p = 4');

legend('Location', 'best');

%% 4b.2 ändra n till h
R = 3;

series_of_n = 30:2:50;
ratio_trap = zeros(length(series_of_n), 1);

for i = 1:length(series_of_n)
    n = series_of_n(i);
    ratio = (trapets(n) - trapets(2*n)) / (trapets(2*n) - trapets(4*n));
    ratio_trap(i) = ratio;
end

ratio_simp = zeros(length(series_of_n), 1);
for i = 1:length(series_of_n)
    n = series_of_n(i);
    ratio = (simpson(n) - simpson(2*n)) / (simpson(2*n) - simpson(4*n));
    ratio_simp(i) = ratio;
end

T = table(series_of_n', ratio_trap, ratio_simp, 'VariableNames', {'n', 'Trapets, kvot', 'Simpson, kvot'});
disp(T);

%% 4c Trapetsregeln i 2D

trap_20 = trapets2d(20);
disp("Volym = " + trap_20)

series_of_n = 80:2:100;

ratio_trap = [];
for n = series_of_n 
    ratio = (trapets2d(n) - trapets2d(2*n)) / (trapets2d(2*n) - trapets2d(4*n));
    ratio_trap = [ratio_trap; ratio];
end

T = table(series_of_n', ratio_trap, 'VariableNames', {'n', 'Kvot'});
disp(T);
%% Funktioner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = trapets(n)
    a = 0; b = 3;
    h = (b-a)/n;
    x = linspace(a, b, n+1);
    
    g = @(r) ((3*r.^3).*exp(-r)./(1+1/3*sin(r*pi/2))).*r;
    V0 = ((3*b^3)*exp(-b)/(1+1/3*sin(b*pi/2)))*b^2*pi;

    V = V0 - h*2*pi*(sum(g(x)) - g(x(1))/2 - g(x(end))/2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = simpson(n)
    a = 0; b = 3;
    h = (b-a)/n;
    x = linspace(a, b, n+1);
    
    g = @(r) ((3*r.^3).*exp(-r)./(1+1/3*sin(r*pi/2))).*r;
    V0 = ((3*b^3)*exp(-b)/(1+1/3*sin(b*pi/2)))*b^2*pi;

    g_x0 = g(x(1));
    g_xn = g(x(end));
    
    odd_i = 2:2:n; even_i = 3:2:n-1;
    
    g_odd = g(x(odd_i));
    g_even = g(x(even_i));
    
    I = (h/3)*(g_x0 + 4*sum(g_odd) + 2*sum(g_even) + g_xn);
    
    V = V0 - 2*pi*I;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = trapets2d(n)
    R = 3; 
    L = 3 * sqrt(2);
    h = L / n;
    x = linspace(-L/2, L/2, n+1);
    y = linspace(-L/2, L/2, n+1);
    
    g_R = (3 * R^3) * exp(-R) / (1 + 1/3 * sin(R * pi / 2));
    g = @(r) g_R - (3 * r.^3) .* exp(-r) ./ (1 + 1/3 * sin(r * pi / 2));
    
    F_xj = zeros(1, length(x));
    
    for j = 1:length(x)
       r = sqrt(x(j).^2 + y.^2);
       f_xy = g(r);
       F_xj(j) = h * (sum(f_xy) - f_xy(1) / 2) - f_xy(end) / 2;
    end
    
    V = h * (sum(F_xj) - F_xj(1) / 2 - F_xj(end) / 2);
end
