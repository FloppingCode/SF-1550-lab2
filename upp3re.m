load dollarkurs.mat
X = USDSEK;
N = length(X);
tt=(1:N)';

%% 3a Linjär modell

A = zeros(N,2);
for i = 1:N
    A(i,1) = 1;
    A(i,2) = i;
end
size(A);
C = A\X;
f = @(t) C(1)+C(2)*t;
t = 1:N;
modely = arrayfun(f,t);
E = dot(X-modely',X-modely')/N; % MSE

%% 3b Linjär + periodisk modell
% skattar L från plot
L = 485;
B = zeros(N,4);
for i = 1:N
   B(i,1) = 1;
   B(i,2) = i;
   B(i,3) = sin(2*pi*i/L);
   B(i,4) = cos(2*pi*i/L);
end
d = B\X;
f_periodic = @(t) d(1) + d(2)*t + d(3)*sin(2*pi*t/L)+d(4)*cos(2*pi*t/L);
newModely = arrayfun(f_periodic,t);
E_periodic = dot(X-newModely',X-newModely')/N; % MSE

%% 3c Icke-linjär modell
% vi tar föregående som startgissning
x0 = [d(1), d(2), d(3), d(4), 486]';
MKMsol = gaussNewton(x0,0.2,X);
f_nonlinear = @(t) MKMsol(1) + MKMsol(2)*t + MKMsol(3)*sin(2*pi*t/MKMsol(5))+MKMsol(4)*cos(2*pi*t/MKMsol(5));
newModely_nonlinear = arrayfun(f_nonlinear,t);
E_nonlinear = dot(X-newModely_nonlinear',X-newModely_nonlinear')/N; % MSE

%% Plotting
figure;
plot(tt, X, '-k', 'LineWidth', 1.5, 'DisplayName', 'Dollarkurs');
hold on;
plot(t, modely, '-r', 'LineWidth', 1.5, 'DisplayName', 'Linjär modell');
plot(t, newModely, '-.g', 'LineWidth', 1.5, 'DisplayName', 'Linjär + periodisk modell');
plot(t, newModely_nonlinear, ':b', 'LineWidth', 1.5, 'DisplayName', 'Ickelinjär modell');

legend('Dollarkurs', 'Linjär modell', 'Linjär + periodisk modell', 'Ickelinjär modell');
xlabel('Time');
ylabel('USD/SEK');
%% Functions
function MKMsol = gaussNewton(x0, tol, X)
    L = length(X);
    x_new = x0;
    disp("first x_n")
    disp(x_new)
    err = inf;
    iter = 0;
    
    while err >= tol
        % linjärisera runt x_n och Låt x_n+1 vara lösningsvektor till MKM
        J = zeros(L, 5);
        F = zeros(L, 1);
        for k = 1:length(X)
           J(k,1) = 1;
           J(k,2) = k;
           J(k,3) = sin(2*pi*k/x_new(5));
           J(k,4) = cos(2*pi*k/x_new(5));
           J(k,5) = 2*pi*k*((x_new(4)*sin(2*pi*k/x_new(5))-x_new(3)*cos(2*pi*k/x_new(5))))/(x_new(5))^2;
           F(k) = (x_new(1) + x_new(2)*k + x_new(3)*sin(2*pi*k/x_new(5)) + x_new(4)*cos(2*pi*k/x_new(5)))-X(k);
        end
        err = dot(F,F)/L;
        if err < tol
            disp(["MKF"])
            dot(F,F)/L
            break
        end
        % update
        sol = J\(-F);
        disp(['sol'])
        disp(sol)
        disp(['x_new'])
        disp(x_new)
        x_new = x_new + sol;
        disp("x_new")
        disp(x_new)
        iter = iter + 1;
        disp(['Iteration: ', num2str(iter), char(10), 'Fel: ', num2str(err), char(10)]);
    end
    MKMsol = x_new;
end
