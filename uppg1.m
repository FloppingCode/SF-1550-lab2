%% 1a) bifogas som pdf


%% 1b -- Visualisering
load eiffel1.mat

[V D] = eig(A);
% V är egenvektorerna.
% D är egenvärdena. 

eigvalArray = ones(1,size(D,2));
%eigvalArray = ones(1,5);
for i=1:size(D,2)
    eigvalArray(i) =  D(i,i);
end
disp("sort");
[eigvalSorted, index] = sort(eigvalArray);

%eigvalSorted
eigvec = [];
for i=1:size(D,2)
    eigvec = cat(2,eigvec,[V(:,index(i))]);
end
y = eigvec;
trussplot(xnod+y(1:2:end), ynod+y(2:2:end), bars);
trussanim(xnod, ynod, bars, y);
%% 1c -- Beräkning av största och minsta egenvärdena

for i= 1:4
    if i == 1
        load eiffel1.mat
    elseif i == 2
        load eiffel2.mat
    elseif i == 3
        load eiffel3.mat
    else
        load eiffel4.mat
    end
    disp("round")
    disp(i)
    tic;
    [mu, iter] = potens(A,10^(-10));
    toc;
    disp([mu, iter])
    mu1 = 0;
    [V, D] = eig(A);
    for j = 1:size(D,2)
        if D(j,j) > mu1
            mu1 = D(j,j);
        end
    end
    disp(mu1)

    
    tic;
    [mu, iter] = inverspotens(A,10^(-10));
    toc;
    disp([mu, iter])
    mu1 = 10^12;
    [V, D] = eig(A);
    for j = 1:size(D,2)
        if D(j,j) < mu1
            mu1 = D(j,j);
        end
    end
    disp(mu1)
end
%% 1d -- Beräkning av andra egenvärden

% Er kod här...
%med shift

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
    if (~exist('shift', 'var'))
       [mu, iter] = potens(A^(-1),tau);
        mu = 1/mu;
    else
        I = eye(size(A,1))
        [mu, iter] = potens((A-shift*I)^(-1),tau);
        mu = 1/(mu) + shift;
    end

end
