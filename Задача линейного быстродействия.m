%%
clear, clc;
 
[A, B, f, x_0, k, p, r, x_1, t_0] = odeInput();
[step, coeff, t_max, n] = progParamsInput();

%На каждой итерации проверяем лежит ли стартовая точка внутри X_0 (решаем с конца)
if (abs(x_1(1) - x_0(1)) <= k/2) && (abs(x_1(2) - x_0(2)) <= k/2)

    disp(['Задача разрешима, T = ', num2str(0)]);
    disp(['Погрешность из УТ_0 = ', num2str(0)]);
else

    %Создаем равномерную сетку направлений
    l(1 : n, 1 : 2) = zeros(n, 2);
    for i = 1 : n 
        l(i, 1) = cos(2*pi*i/n);
        l(i, 2) = sin(2*pi*i/n);
    end

    [t_1_best, x_memory, s_memory, u_memory, i_best, t] = trajectoriesBuilding(A, B, f, x_0, k, p, r, x_1, t_0, step, coeff, t_max, n, l);

    %Результаты
    if i_best

        s_best = s_memory(end, 2 * i_best - 1 : 2 * i_best);

        disp(t_1_best)
        disp(['Задача разрешима, T = ', num2str(t_1_best - t_0)]);
        x_best_end = x_memory(end, 2 * i_best - 1 : 2 * i_best);
        disp(['Погрешность из УТ_0 = ', num2str(UT_0(s_best, x_0, k, x_best_end))])
    else
    
        s_best = [];
        disp('Задача не разрешима');
    end

    %Визуализация
    t_best = [transpose([t_0 : 1/coeff : t(1)]); t(2 : end)];
    draw_graphics(x_memory, s_memory, u_memory, x_0, k, p, r, x_1, t_0, t_1_best, i_best, s_best, t_best);    
end
  
%%
function [A, B, f, x_0, k, p, r, x_1, t_0] = odeInput()
    A = input('Введите матрицу A\n');
    while ~check_2x2(A)
        A = input('Введите матрицу A\n');
    end
    
    B = input('Введите матрицу B\n');
    while ~check_2x2(B)
        B = input('Введите матрицу B\n');
    end
    
    f = input('Введите матрицу f\n');
    while ~check_2x1(f)
        f = input('Введите матрицу f\n');
    end
    
    disp('Введите параметры начального множества X_0:');
    x_0 = input('1) координаты центра квадрата\n');
    while ~check_1x2(x_0)
        x_0 = input('1) координаты центра квадрата\n');
    end
    
    k = input('2) длина стороны квадрата\n');
    while (k <= 0)
        k = input('2) длина стороны квадрата\n');
    end
    
    disp('Введите параметры множества допустимых управлений P:');
    p = input('1) координаты центра P \n');
    while ~check_1x2(p)
        p = input('1) координаты центра P \n');
    end
    p = transpose(p);
    
    r = input('2) радиус круга \n');
    while (r < 0)
        r = input('2) радиус круга \n');
    end
    
    disp('Введите параметры конечного множества X_1:');
    x_1 = input('координаты конечной точки\n');
    while ~check_1x2(x_1)
        x_1 = input('координаты конечной точки\n');
    end
    
    disp('Введите начальный момент времени t_0:');
    t_0 = input('');
end

function q = check_2x2(A)
    if (size(A, 1) == 2 && size(A, 2) == 2)
        q = 1;
    else 
        q = 0;
    end
end

function q = check_2x1(f)
    if (size(f, 1) == 2 && size(f, 2) == 1)
        q = 1;
    else 
        q = 0;
    end
end

function q = check_1x2(x)
    if (size(x, 1) == 1 && size(x, 2) == 2)
        q = 1;
    else 
        q = 0;
    end
end

function [step, coeff, t_max, n] = progParamsInput()
    step = input('Введите шаг численного метода\n');
    coeff = input('Введите необходимое количество отсчётов за шаг\n');
    t_max = input('Введите ограничение сверху на время\n');
    n = input('Введите число элементов в сетке по начальным значениям сопряженной переменной\n');
end

function [t_1_best, x_memory, s_memory, u_memory, i_best, t] = trajectoriesBuilding(A, B, f, x_0, k, p, r, x_1, t_0, step, coeff, t_max, n, l)
   
    t_1_best = t_max;
    s_memory = zeros(t_max * coeff + 1, 2);
    u_memory = zeros(t_max * coeff + 1, 2);
    x_memory = zeros(t_max * coeff + 1, 2);
    i_best = 0;
    for i = 1 : n

        t_left = t_0;
        t_right = t_0 + step;
        s(1, :) = l(i, :);
        x(1, :) = x_1;
        i_left = 1;
        u = [];
        while(t_right <= t_1_best)       

            [t, s_curr, s] = s_search(A, t_left, t_right, s, i_left, coeff, []);
            [u, u_curr] = u_search(B, t_left, t_0, r, p, s_curr, u);        
            [t, x] = x_search(t, x, i_left, u_curr, A, B, f, x_0, k);

            if (t(end) < t_right)
                t_1_best = t(end);
                i_best = i;
            end

            t_left = t_right;
            t_right = t_right + step;
            i_left = size(s, 1);

        end
        
        %Так как последний отрезок интегрирования целиком не лежит в рассматриваемом временном 
        %промежутке, обрезаем его до t_1_best, а далее интегрируем как и раньше        
        if (t_right > t_1_best) && (t_right < t_1_best + step)      

            t = t_left : 1/coeff : t_right;
            t = t(t < t_1_best);
            t(end + 1) = t_1_best;
            
            [t, s_curr, s] = s_search(A, t_left, t_1_best, s, i_left, coeff, t);
            [u, u_curr] = u_search(B, t_left, t_0, r, p, s_curr, u);        
            [t, x] = x_search(t, x, i_left, u_curr, A, B, f, x_0, k); 

            if (t(end) < t_1_best)
                t_1_best = t(end);
                i_best = i;
            end    

        end 

        %Сохранение и очистка
        s_memory = saving(s_memory, s);
        u_memory = saving(u_memory, u);
        x_memory = saving(x_memory, x);  
        s = [];
        u = [];
        x = [];   
    end
end

function [t, s_curr, s] = s_search(A, t_left, t_right, s, i_left, coeff, t)

    dsdt = @(t, s) [-A(1)*s(1) - A(2)*s(2);
                    -A(3)*s(1) - A(4)*s(2)];
    if isempty(t)
        tspan = [t_left : 1 / coeff : t_right];
    else 
        %Соответствует последнему отрезку времени
        tspan = t;
    end
    s0 = s(i_left, :);
    [t, s_curr] = ode45(dsdt, tspan, s0); 
    
    s = [s; s_curr(2 : end, :)];
end

function [u, u_curr] = u_search(B, t_left, t_0, r, p, s_curr, u)
    
    eps = 0.001;
    %Случай "плохого" s
    if s_curr(1) == 0
        s_curr(1) = eps;
    else 
        if s_curr(2) == 0
            s_curr(2) = eps;
        end
    end    
    %Случай вырожденной B
    if (B(1)*B(4) - B(3)*B(2) == 0)
            B = B - [eps 0; 0 eps];
    end
    
    rho_u = transpose(B) * transpose(s_curr);       
    u_curr = sign(rho_u) * 0.5 + rho_u ./ sqrt(sum(rho_u.^2, 1)) .* r + p;
    u_curr = transpose(u_curr); 
    if (t_left == t_0)
        u = u_curr;
    else 
        u = [u; u_curr(2 : end, :)];
    end
end

function [t, x] = x_search(t, x, i_left, u_curr, A, B, f, x_0, k)
    
    tspan = t;
    t_old = t;
    x0 = x(i_left, :);
    [t, x_curr] = ode45(@(t, x) myode(t, x, t_old, u_curr, A, B, f), tspan, x0, odeset('Events', @(t, x) X_0_Event(t, x, x_0, k), 'RelTol',1e-2, 'AbsTol',1e-4)); 
    x = [x; x_curr(2 : end, :)];
end

function dxdt = myode(t, x, t_old, u_curr, A, B, f)
    u_curr = interp1(t_old, u_curr, t);   
    dxdt = A * x + B * transpose(u_curr) + f;
end

%Функция события: попали в начальное мн-во (решение ищется с конца)
function [value, isterminal, direction] = X_0_Event(t, x, x_0, k)
    value = (abs(x(1) - x_0(1)) > k/2) | (abs(x(2) - x_0(2)) > k/2);  
    isterminal = 1;   
    direction = 0;   
end

function arr_memory = saving(arr_memory, arr)
    arr_memory = arr_memory(1 : size(arr, 1), :);
    arr_memory(:, end - 1 : end) = arr;
    arr_memory(:, end + 1 : end + 2) = zeros(size(arr_memory, 1), 2);
end

%Погрешность условия трансверсальности на левом конце (так определяем точность найденного решения)
function res = UT_0(s_best, x_0, k, x_best_end)
    %разность между сдвигами многообразий
    s_best = s_best ./ norm(s_best);
    res = abs( (abs(s_best(1)) + abs(s_best(2))) * (k/2) + (-s_best) * transpose(x_0) - x_best_end * transpose(-s_best) );
end

function draw_graphics(x_memory, s_memory, u_memory, x_0, k, p, r, x_1, t_0, t_1_best, i_best, s_best, t_best)  
    if i_best
        
        disp('Введите номера осей, графики результатов в которых вы хотите увидеть. Вторым параметром после номера для каждой оси обозначьте полноту вывода (a - только оптимальную, b -все ).')
        disp('1) (x1, x2)')
        disp('2) (u1, u2)')
        disp('3) (t, x1)')
        disp('4) (t, x2)')
        disp('5) (t, s1)')
        disp('6) (t, s2)')
        disp('7) (t, u1)')
        disp('8) (t, u2)')
        disp('Пример ввода:')
        disp("[1 'a'; 1 'b'; 2 'a'; 2 'b'; 3 'a'; 3 'b'; 4 'a'; 4 'b'; 5 'a'; 5 'b'; 6 'a'; 6 'b'; 7 'a'; 7 'b'; 8 'a'; 8 'b';]")
        params = input('');
    else       
        
        params = [1 'b'; 2 'b'; 3 'b'; 4 'b'; 5 'b'; 6 'b'; 7 'b'; 8 'b'];
    end    
    
    a = params(:, 1);
    b = params(:, 2);
    F = figure;
    G = uitabgroup('Parent', F); 
    for i = 1 : size(params, 1)
        tab = uitab(G);
        ax = axes('parent', tab); 
        
        switch a(i)
            case 1
                name1 = 'x_1';
                name2 = 'x_2';
                hold(tab.Children(end), 'on');
               
                if (b(i) == 'a')
                   
                    plot(tab.Children(end), x_memory(:, 2*i_best-1), x_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    rectangle(tab.Children(end), 'Position', [x_0(1) - k/2, x_0(2) - k/2, k, k], 'FaceColor', [0, 0.5, 0.5]);
                    plot(tab.Children(end), x_1(1), x_1(2), '.', 'Color', 'black', 'MarkerSize', 5, 'DisplayName', 'X_1');
                    plot(tab.Children(end), x_0(1), x_0(2), '.', 'Color', [0, 0.5, 0.5], 'MarkerSize', 5, 'DisplayName', 'X_0');
                
                    s_best = -s_best;
                    c = -(s_best(1) * (x_0(1) + sign(s_best(1)) * k/2) + s_best(2) * (x_0(2) + sign(s_best(2)) * k/2));
                    x_tmp = [x_0(1) - k : 0.1 : x_0(1) + k];
                    if (s_best(2) ~= 0)
                        y = (s_best(1) * x_tmp + c) ./ (-s_best(2));
                    else
                        y = x_0(1) + sign(s_best(1)) * k/2 + 0 .* x_tmp;
                    end
                    s_best = -s_best;
                    
                    hold(tab.Children(end), 'off');                 
                else                   
                    plot(tab.Children(end), x_memory(:, 1 : 2 : end - 2), x_memory(:, 2 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), x_memory(:, 2*i_best-1), x_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    else
                        xlim([min(x_1(1),x_0(1) - k/2 ) - 1, max(x_1(1),x_0(1) + k/2 ) + 1]);
                        ylim([min(x_1(2),x_0(2) - k/2 ) - 1, max(x_1(2),x_0(2) + k/2 ) + 1]); 
                    end
                    rectangle(tab.Children(end), 'Position', [x_0(1) - k/2, x_0(2) - k/2, k, k], 'FaceColor', [0, 0.5, 0.5]);
                    plot(tab.Children(end), x_1(1), x_1(2), '.', 'Color', 'black', 'MarkerSize', 5, 'DisplayName', 'X_1');
                    plot(tab.Children(end), x_0(1), x_0(2), '.', 'Color', [0, 0.5, 0.5], 'MarkerSize', 5, 'DisplayName', 'X_0');
                    
                    if i_best
                        s_best = -s_best;
                        c = -(s_best(1) * (x_0(1) + sign(s_best(1)) * k/2) + s_best(2) * (x_0(2) + sign(s_best(2)) * k/2));
                        x_tmp = [x_0(1) - k : 0.1 : x_0(1) + k];
                        if (s_best(2) ~= 0)
                            y = (s_best(1) * x_tmp + c) ./ (-s_best(2));
                        else
                            y = x_0(1) + sign(s_best(1)) * k/2 + 0 .* x_tmp;
                        end
                        s_best = -s_best;
                    end
                    
                    hold(tab.Children(end), 'off');
                end
            case 2
                name1 = 'u_1';
                name2 = 'u_2';
                hold(tab.Children(end), 'on');
                xlim([p(1) - 0.5 - r - 1, p(1) + 0.5 + r + 1]);
                ylim([p(2) - 0.5 - r - 1, p(2) + 0.5 + r + 1]); 
                
                if (b(i) == 'a')
                    rectangle(tab.Children(end), 'Position', [p(1)-0.5-r, p(2)-0.5-r, 2*r+1, 2*r+1], 'Curvature', r/(r + 0.5), 'FaceColor', 'y');
                    plot(tab.Children(end), u_memory(:, 2*i_best-1), u_memory(:, 2*i_best), '.', 'Color', 'red', 'MarkerSize', 10, 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    plot(tab.Children(end), p(1), p(2), '.', 'Color', 'yellow', 'MarkerSize', 5, 'DisplayName', 'P');
                    hold(tab.Children(end), 'off');
                else
                    rectangle(tab.Children(end), 'Position', [p(1)-0.5-r, p(2)-0.5-r, 2*r+1, 2*r+1], 'Curvature', r/(r + 0.5), 'FaceColor', 'y');
                    plot(tab.Children(end), u_memory(:, 1 : 2 : end - 2), u_memory(:, 2 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), u_memory(:, 2*i_best-1), u_memory(:, 2*i_best), '.', 'Color', 'red', 'MarkerSize', 10, 'DisplayName', 'Opt');
                    end
                    plot(tab.Children(end), p(1), p(2), '.', 'Color', 'yellow', 'MarkerSize', 5, 'DisplayName', 'P');
                    hold(tab.Children(end), 'off');
                end
            case 3
                name1 = 't';
                name2 = 'x_1';
                hold(tab.Children(end), 'on');
                xlim([t_0 - 1, t_1_best + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), rot90(rot90(t_best)), x_memory(:, 2*i_best-1), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), rot90(rot90(t_best)), x_memory(:, 1 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), rot90(rot90(t_best)), x_memory(:, 2*i_best-1), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 4
                name1 = 't';
                name2 = 'x2';
                hold(tab.Children(end), 'on');
                xlim([t_0 - 1, t_1_best + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), rot90(rot90(t_best)), x_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), rot90(rot90(t_best)), x_memory(:, 2 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), rot90(rot90(t_best)), x_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 5
                name1 = 't';
                name2 = 's_1';
                hold(tab.Children(end), 'on');
                xlim([t_0 - 1, t_1_best + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), rot90(rot90(t_best)), s_memory(:, 2*i_best-1), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), rot90(rot90(t_best)), s_memory(:, 1 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), rot90(rot90(t_best)), s_memory(:, 2*i_best-1), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 6
                name1 = 't';
                name2 = 's_2';
                hold(tab.Children(end), 'on');
                xlim([t_0 - 1, t_1_best + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), rot90(rot90(t_best)), s_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), rot90(rot90(t_best)), s_memory(:, 2 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), rot90(rot90(t_best)), s_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 7
                name1 = 't';
                name2 = 'u_1';
                hold(tab.Children(end), 'on');
                xlim([t_0 - 1, t_1_best + 1]);
                ylim([p(2)-0.5-r-1, p(2)+0.5+r+1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), rot90(rot90(t_best)), u_memory(:, 2*i_best-1), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), rot90(rot90(t_best)), u_memory(:, 1 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), rot90(rot90(t_best)), u_memory(:, 2*i_best-1), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 8
                name1 = 't';
                name2 = 'u_2';
                hold(tab.Children(end), 'on');
                xlim([t_0 - 1, t_1_best + 1]);
                ylim([p(2)-0.5-r-1, p(2)+0.5+r+1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), rot90(rot90(t_best)), u_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), rot90(rot90(t_best)), u_memory(:, 2 : 2 : end - 2), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), rot90(rot90(t_best)), u_memory(:, 2*i_best), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
        end
        
        xlabel(tab.Children(end), name1);
        ylabel(tab.Children(end), name2); 
        
        if (b(i) == 'a')
            type = 'opt';
        else
            type = 'all';
        end        
        
        tab.Title = strcat(type, ' (', name1, ',', name2, ')');
    end
end

%Измельчение сетки в определённом угловом секторе и построение соответствующих траекторий
function [x_memory, s_memory, u_memory, s_best, i_best, t_best, t_1_best] = local_add(min_angle, max_angle, n, A, B, f, x_0, k, p, r, x_1, t_0, step, coeff, t_max, s_old, u_old, x_old, i_best_old, t_1_best_old)
    
    l(1 : n, 1 : 2) = zeros(n, 2);
    for i = 1 : n
        l(i, 1) = cos(min_angle + (max_angle - min_angle)*i/n);
        l(i, 2) = sin(min_angle + (max_angle - min_angle)*i/n);
    end

    [t_1_best, x_memory, s_memory, u_memory, i_best, t] = trajectoriesBuilding(A, B, f, x_0, k, p, r, x_1, t_0, step, coeff, t_max, n, l);
  
    [x_memory, s_memory, u_memory, s_best, i_best, t_best] = results(t_1_best, t_1_best_old, x_memory, s_memory, u_memory, i_best, t, i_best_old, x_old, s_old, u_old, coeff, x_0, k, p, r, x_1, t_0);
end

%Измельчение сетки в целое число раз глобально (по всем направлениям) и построение соответствующих траекторий
function [x_memory, s_memory, u_memory, s_best, i_best, t_best, n_new, t_1_best] = global_add(n_old, n_new, A, B, f, x_0, k, p, r, x_1, t_0, step, coeff, t_max, s_old, u_old, x_old, i_best_old, t_1_best_old)
    
    k = n_new / n_old;
    n = n_new - n_old;
    l(1 : n_new, 1 : 2) = zeros(n_new, 2);
    for i = 1 : n_new 
        l(i, 1) = cos(2*pi*i/n_new);
        l(i, 2) = sin(2*pi*i/n_new);
    end
    l(k : k : end, :) = [];
    
    [t_1_best, x_memory, s_memory, u_memory, i_best, t] = trajectoriesBuilding(A, B, f, x_0, k, p, r, x_1, t_0, step, coeff, t_max, n, l);
    
    [x_memory, s_memory, u_memory, s_best, i_best, t_best] = results(t_1_best, t_1_best_old, x_memory, s_memory, u_memory, i_best, t, i_best_old, x_old, s_old, u_old, coeff, x_0, k, p, r, x_1, t_0);
end

function [x_memory, s_memory, u_memory, s_best, i_best, t_best] = results(t_1_best, t_1_best_old, x_memory, s_memory, u_memory, i_best, t, i_best_old, x_old, s_old, u_old, coeff, x_0, k, p, r, x_1, t_0)

    if (i_best ~= 0 || i_best_old ~= 0)

        if t_1_best_old < t_1_best && t_1_best_old ~= 0
            i_best = i_best_old;
            t_1_best = t_1_best_old;
        else
            i_best = i_best + (size(x_old, 2) - 2) / 2;
        end
        m_min = min(size(x_old, 1), size(x_memory, 1));
        s_old = s_old(1 : m_min, 1 : end - 2);
        s_memory = s_memory(1 : m_min, :);
        u_old = u_old(1 : m_min, 1 : end - 2);
        u_memory = u_memory(1 : m_min, :);
        x_old = x_old(1 : m_min, 1 : end - 2);
        x_memory = x_memory(1 : m_min, :);
        s_memory = horzcat(s_old, s_memory);
        u_memory = horzcat(u_old, u_memory);
        x_memory = horzcat(x_old, x_memory);

        s_best = s_memory(end, 2 * i_best - 1 : 2 * i_best);

        disp(['Задача разрешима, T = ', num2str(t_1_best - t_0)]);

        x_best_end = x_memory(end, 2 * i_best - 1 : 2 * i_best);
        disp(['Погрешность из УТ_0 = ', num2str(UT_0(s_best, x_0, k, x_best_end))])
    else           

        s_best = [];
        disp('Задача не разрешима');
    end
    
    t_best = [transpose([t_0 : 1/coeff : t(1)]); t(2 : end)];    
    draw_graphics(x_memory, s_memory, u_memory, x_0, k, p, r, x_1, t_0, t_1_best, i_best, s_best, t_best); 
end
