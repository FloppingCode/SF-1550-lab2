%% 1a) bifogas som pdf


%% 1b -- Visualisering
load eiffel1.mat

[V D] = eig(A);
% V är egenvektorerna.
% D är egenvärdena. 

eigvalArray = zeros(1,size(D,2));
for i=1:size(D,2)
    eigvalArray(i) =  D(i,i);
end
[eigvalSorted, index] = sort(eigvalArray);

eigvec = [];
for i=1:size(D,2)
    eigvec = cat(2,eigvec,[V(:,index(i))]);
end

trussplot(xnod, ynod, bars); %ingen förskjutning
for i=1:4
    y = eigvec(:,i);
    hold on;
    if i == 1
       trussplot(xnod+y(1:2:end), ynod+y(2:2:end), bars, 'r');
    elseif i == 2
       trussplot(xnod+y(1:2:end), ynod+y(2:2:end), bars, 'g');
    elseif i == 3
       trussplot(xnod+y(1:2:end), ynod+y(2:2:end), bars, 'b');
    else
       trussplot(xnod+y(1:2:end), ynod+y(2:2:end), bars, 'y');
    end
end
% does the last one
trussanim(xnod, ynod, bars, y);

% skriv upp deras frekvens. Det är bara en funktion av egenvärdet tror jag
%f = @(k) 1/k
%frekvenser = arrayfun()
%% 1c -- Beräkning av största och minsta egenvärdena

for i= 1:4
    format long
    if i == 1
        load eiffel1.mat
    elseif i == 2
        load eiffel2.mat
    elseif i == 3
        load eiffel3.mat
    else
        load eiffel4.mat
    end
    disp(["round: ", i])
    disp(["time eig"])
    tic;
    [V, D] = eig(A);
    toc;
    eigvalArray = zeros(1,size(D,2));
    for j=1:size(D,2)
        eigvalArray(j) =  D(j,j);
    end
    [eigvalSort, index] = sort(eigvalArray);

    disp(["time potens"])
    tic;
    [mu, iter] = potens(A,10^(-10));
    toc;
    disp(["iter: ", iter])
    disp(["approx eigenvalue: ", mu])
    disp(["actual value: ", eigvalSort(end)])
    disp(["ratio lambda(2)/lambda(1): ", eigvalSort(end-1)/eigvalSort(end)])
    if i == 4
        disp(["ratio lambda(3)/lambda(2): ", eigvalSort(end-2)/eigvalSort(end-1)])
    end
    disp(["eig 2 och 1", eigvalSort(end-1), eigvalSort(end)])
    disp(["time invpotens"])
    tic;
    [mu, iter] = inverspotens(A,10^(-10));
    toc;
    disp(["iter: ", iter])  
    disp(["approx eigenvalue: ", mu])
    %finds min eigenval
    disp(["actual value: ", eigvalSort(1)])
    disp(["ratio lambda(n)/lambda(n-1): ", eigvalSort(1)/eigvalSort(2)])
    
end
%% 1d -- Beräkning av andra egenvärden
load eiffel1.mat

[V, D] = eig(A);
eigvalArray = zeros(1,size(D,2));
for i=1:size(D,2)
    eigvalArray(i) =  D(i,i);
end
[eigvalSort, index] = sort(eigvalArray);
% noterar att dessa egenvärden ligger nära 67 och påverkar konvergensen.
% ty två största egenvärden till (A-67*I)^(-1) blir nära varandra. 
disp(["egenvärden nära 67"])
eigvalSort(69)
eigvalSort(70)
% två nära 10
disp(["egenvärden nära 10"])
eigvalSort(18)
eigvalSort(19)
% två nära 50
disp(["egenvärden nära 50"])
eigvalSort(53)
eigvalSort(54)

[mu iter] = inverspotens(A,10^(-10),10);
disp(["eig for 10: ", mu, " iter: ", iter])
%disp(["ratio lambda(2)/lambda(1): ", eigvalSort(end-1)/eigvalSort(end)])
[mu iter] = inverspotens(A,10^(-10),50);
disp(["eig for 50: ", mu, " iter: ", iter])
%[mu iter] = inverspotens(A,10^(-10),66.75);
[mu iter] = inverspotens(A,10^(-10),67);
disp(["eig for 67: ", mu, " iter: ", iter])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = potens(A,tau)

%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (skälär)
%
%  Utdata:
%
%  mu - största egenvärdet till A (skalär)
%  iter - antal iterationer som använts (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = sparse(A);
y_prev = ones(size(A,1),1);
fel = 10^9;
iter = 0;
    while fel >= tau
        v_new = A*y_prev;
        mu = dot(v_new,y_prev);
        if iter > 0
            fel = abs(mu-eigval_prev);
        end
        eigval_prev = mu;
        y_new = (v_new)/norm(v_new);
        y_prev = y_new;
        iter = iter + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = inverspotens(A,tau,shift)

%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (skälär)
%
%  Utdata:
%
%  mu - minsta egenvärdet till A (skalär)
%  iter - antal iterationer som använts (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %A = lu(A);
    if (~exist('shift', 'var'))
       [mu, iter] = potens(A^(-1),tau);
        mu = 1/mu;
    else
        I = eye(size(A,1));
        [mu, iter] = potens((A-shift*I)^(-1),tau);
        mu = 1/(mu) + shift;
    end

end
