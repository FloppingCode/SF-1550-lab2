Iexact = 6.231467927023725;  % Ett noggrannt värde för I

%% 5a Trapetsregeln i 10 dimensioner

% Er kod här...
tic
format long
I = trapets10d(7)
toc
Trapfelet = abs(I-Iexact)
disp("Felet blev:")
disp(Trapfelet)

%% 5b Monte-Carlo
format long

Hypervolym = 1.2^10;
N = 10^6;
runs = 5

all_integ = zeros(N, runs);
all_errors = zeros(N, runs);

for run = 1:runs
    dp = 0;
    Is2 = zeros(1, N);
    for k = 1:N
        samp = (1.2) * rand(1, 10);
        f = @(x) exp(prod(x));
        I_k = (Hypervolym * f(samp) + dp * (k - 1)) / k;
        dp = I_k;
        Is2(k) = I_k;
    end
    all_integ(:, run) = Is2';
    all_errors(:, run) = abs(Iexact - Is2)';
end

figure;
hold on;
for run = 1:runs
    plot(1:length(all_integ), all_integ(:, run), 'DisplayName', ['Run ', num2str(run)]);
end

xlabel('Antalet punkter');
ylabel('Integralens värde');
title('Monte Carlo - integral');
legend('Location', 'Best');
grid on;

figure;
hold on;
for run = 1:runs
    loglog(1:length(all_errors), all_errors(:, run), 'DisplayName', ['Run ', num2str(run)]);
end
arr = ones(N, 1) * Trapfelet;
plot(1:N, arr,'DisplayName', "trapfelet");
axis([0 N -10^-3 10^-3]);
xlabel('Antalet punkter')
ylabel('Fel');
title('Monte Carlo - fel');
legend('Location', 'Best');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = trapets10d(n)

%  Indata:
%
%  n  - antal delintervall i varje koordinatriktning (skalär)
%
%  Utdata:
%
%  I - integralvärdet (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=1.2;

h = L/n;

x = 0:h:L;

I1 = zeros(n+1,1);
I2 = zeros(n+1,1);
I3 = zeros(n+1,1);
I4 = zeros(n+1,1);
I5 = zeros(n+1,1);
I6 = zeros(n+1,1);
I7 = zeros(n+1,1);
I8 = zeros(n+1,1);
I9 = zeros(n+1,1);
I10 = zeros(n+1,1);

for j1=0:n
    for j2=0:n
        for j3=0:n
            for j4=0:n
                for j5=0:n
                    for j6=0:n
                        for j7=0:n
                            for j8=0:n
                                for j9=0:n
                                    for j10=0:n
                                        I10(j10+1) = exp(j1*j2*j3*j4*j5*j6*j7*j8*j9*j10*h^10);
                                    end
                                    I9(j9+1) = h*(sum(I10) - I10(1)/2 - I10(end)/2);
                                end
                                I8(j8+1) = h*(sum(I9) - I9(1)/2 - I9(end)/2);
                            end
                            I7(j7+1) = h*(sum(I8) - I8(1)/2 - I8(end)/2);
                        end
                        I6(j6+1) = h*(sum(I7) - I7(1)/2 - I7(end)/2);
                    end
                    I5(j5+1) = h*(sum(I6) - I6(1)/2 - I6(end)/2);
                end
                I4(j4+1) = h*(sum(I5) - I5(1)/2 - I5(end)/2);
            end
            I3(j3+1) = h*(sum(I4) - I4(1)/2 - I4(end)/2);
        end
        I2(j2+1) = h*(sum(I3) - I3(1)/2 - I3(end)/2);
    end
    I1(j1+1) = h*(sum(I2) - I2(1)/2 - I2(end)/2);
end
I = h*(sum(I1) - I1(1)/2 - I1(end)/2);

end
