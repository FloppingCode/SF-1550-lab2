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
E = dot(X-modely',X-modely')/N % MSE
figure;
plot(t, modely)
hold on;
plot(t, X, 'DisplayName', 'Financial market data model');

%% 3b Linjär + periodisk modell
% skattar L från plot
L = 485;
B = zeros(N,2);
for i = 1:N
   B(i,1) = 1;
   B(i,2) = i;
   B(i,3) = sin(2*pi*i/L);
   B(i,4) = cos(2*pi*i/L);
end
d = B\X
t = 1:N;
f = @(t) d(1) + d(2)*t + d(3)*sin(2*pi*t/L)+d(4)*cos(2*pi*t/L);
newModely = arrayfun(f,t);
E = dot(X-newModely',X-newModely')/N % MSE
figure; 
plot(t, newModely);
hold on;
plot(t, X, 'DisplayName', 'Financial market data model');
%% 3c Icke-linjär modell
% vi tar föregående som startgissning
x0 = [d(1), d(2), d(3), d(4), 486]'
%x0 = [1,1,1,1,1]'
disp(x0)
MKMsol = gaussNewton(x0,0.2,X)
f = @(t) MKMsol(1) + MKMsol(2)*t + MKMsol(3)*sin(2*pi*t/MKMsol(5))+MKMsol(4)*cos(2*pi*t/MKMsol(5));
figure;
newModely = arrayfun(f,t);
plot(t, newModely);
hold on;
plot(t, X, 'DisplayName', 'Financial market data model');

function MKMsol = gaussNewton(x0, tol, X)
    L = length(X);
    x_new = x0;
    disp("first x_n")
    disp(x_new)
    err = inf;
    iter = 0;
    
    while err >= tol
    %while iter < tol
        % linjärisera runt x_n och Låt x_n+1 vara lösningsvektor till MKM
        J = zeros(L, 5);
        F = zeros(L, 1);
        for k = 1:length(X)
           J(k,1) = 1;
           J(k,2) = k;
           J(k,3) = sin(2*pi*k/x_new(5));
           J(k,4) = cos(2*pi*k/x_new(5));
           J(k,5) = 2*pi*k*((x_new(4)*sin(2*pi*k/x_new(5))-x_new(3)*cos(2*pi*k/x_new(5))))/(x_new(5))^2;
           %delta_X(k) = (x_new(1) + x_new(2)*k + x_new(3)*sin(2*pi*k./x_new(5)) + x_new(4)*cos(2*pi*k./x_new(5)));
           F(k) = (x_new(1) + x_new(2)*k + x_new(3)*sin(2*pi*k/x_new(5)) + x_new(4)*cos(2*pi*k/x_new(5)))-X(k);
        end
        err = dot(F,F)/L
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
        %sol = sol'
        x_new = x_new + sol;

        disp("x_new")
        disp(x_new)
        iter = iter + 1;
        disp(['Iteration: ', num2str(iter), char(10), 'Fel: ', num2str(err), char(10)]);
    end
    MKMsol = x_new;
end