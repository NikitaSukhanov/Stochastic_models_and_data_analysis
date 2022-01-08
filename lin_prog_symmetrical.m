function [x, w] = lin_prog_symmetrical(K, inf_b, sup_b)
    m = size(K, 1);
    n = size(K, 2);
    rad_b = (sup_b - inf_b) / 2;
    mid_b = (sup_b + inf_b) / 2;
    
    % LP canonic form prepare
    A = zeros(2*m, n + m);
    for i = 1:m
        A(i, i + n) = -rad_b(i);
        A(i + m, i + n) = -rad_b(i); 
        for j = 1:n
            A(i, j) = K(i, j);
            A(i + m, j) = -K(i, j);
        end
    end
    
    f = cat(1, zeros(n, 1), ones(m, 1));
    c = cat(1, mid_b, -mid_b);
    
    % Lp solve
    linprog_solution = linprog(f, A, c, [], [], cat(1, zeros(n, 1), ones(m, 1)), []);
    x = linprog_solution(1:n);
    w = linprog_solution(n + 1: n + m);
    # display(flag);
    # display(f'*linprog_solution);
end
