    % 參數設定
tf = 0.28;     % m
Lf = 2;        % m
Lj = 0.5;      % m
Ls = 3;        % m
theta_fs = 50; % 度
alphas = [1, 1/2, 1/3];  % alpha 的值
theta_as = linspace(-30, 30, 10000);  % degrees
m = 1.40;
n = 3.49;
p = 1.31;
q = 4.14;


for i = 1:length(alphas)
    alpha = alphas(i);

    % 計算
    Nn_max = 3 * pi + 2 + (tf / Lf) * (alpha + (1 + alpha) / sqrt(2));
    Nt_max = 2 * alpha + 15 * (tf / Lf);
    Nm_max = (pi / 2) * (1 + (tf / Lf) ^ 2);
    c1 = sind(theta_fs + theta_as);
    c2 = cosd(theta_fs + theta_as);
    c3 = sind(theta_fs + theta_as) * ((Lj / Lf) + (Ls / Lf) * cosd(theta_fs) - (1 / 2)) - cosd(theta_fs + theta_as) * sind(theta_fs * (Ls / Lf));
    Nem = Nm_max ./ abs(c3);
    Nen = Nn_max ./ abs(c1);
    Net = Nt_max ./ abs(c2);

    % 找到每個 x 對應的 y 值的最小值
    Noninteractive = min([Nem; Nen; Net]);
    
end

Lf = 2;        % m
Ls = 3;        % m
Lj = 0.5;      % m
Af = 6;        % m^2
tf = 0.28;     % m
Wf = 3;        % m
dl = 0.073;    % m
theta_fs = 50; % 度
Su = 2;        % kPa
Su0 = 2;       % kPa
k = 0;         % kPa/m
alpha = 1;
Ncl = 12;
Nn_max = 11.98;
Nt_max = 4.39;
Nm_max = 1.645;
m = 1.56;
n = 4.19;
p = 1.57;
q = 4.43;
Za = 0;       % Zai+1 = Zai + deltaZi 
theta_a = 0.1;
theta_s = 0;
theta_as = 0; % theta_as = theta_a - theta_s 
En = 1;
theta_0 = 0;
theta_f = 50;
delta_t = 0.01;
last_Ne = 0;
last_theta_as = 0;
last_theta_a = 0;
last_theta_s = 0;
last_hat_z = 0;
delta_x = zeros(10000, 1);
delta_z = zeros(10000, 1);
x_positions = zeros(10000, 1);
z_positions = zeros(10000, 1);
i = 1;

% 初始化第一個位置
x_positions(1) = 0;
z_positions(1) = 0;

% 主循環
while theta_f ~= 0 && i <= 100
    fprintf('---------------------------------------------------------------------\n');
    fprintf('這是第 %d 次迴圈 \n', i);
    fprintf('theta_f = %.4f, theta_a = %.4f, Za = %.4f\n',  theta_f, theta_a, Za);
 
    % 計算 C1, C2, C3
    C1 = sind(theta_as + theta_fs);
    C2 = cosd(theta_as + theta_fs);
    C3 = (sind(theta_fs + theta_as) * (Lj / Lf + (Ls / Lf) * cosd(theta_fs) - 0.5) - cosd(theta_fs + theta_as) * sind(theta_fs) * Ls / Lf);
    Nm = (Ne * C3);
    Nt = (Ne * C2);
    Nn = (Ne * C1);
    
    fprintf('C1 = %f, C2 = %f, C3 = %f\n', C1, C2, C3);
    fprintf('Nm = %f, Nt = %f, Nn = %f\n', Nm, Nt, Nn);
    
    delta_n = ((((Nt_max/Nn_max) * (p * q / n)) / (abs(Nm / Nm_max)^m + abs(Nt / Nt_max)^n)^((1/p) - 1)) * ((abs(Nn / Nn_max)^(q-1)) / ((Nt / Nt_max)^(n-1)))) * delta_t;
    
    fprintf('delta_n = %f\n', delta_n);
    
    % 計算 delta_x 和 delta_z
    delta_x(i) = delta_t * cosd(theta_f) + delta_n * sind(theta_f);
    delta_z(i) = delta_t * sind(theta_f) - delta_n * cosd(theta_f);
    
    % 更新累積位置
    x_positions(i) = sum(delta_x(1:i));
    z_positions(i) = sum(delta_z(1:i));
    
    fprintf('delta_x(%d) = %f, delta_z(%d) = %f\n', i, delta_x(i), i, delta_z(i));
    
    % 計算 delta_beta
    delta_beta = ((C3 * m * Nt_max * (abs(Nm) / Nm_max)^(m-1)) / (abs(C3) * n * Nm_max * (abs(Nt) / Nt_max)^(n-1))) * (delta_t / Lf);
    fprintf('delta_beta = %f\n', delta_beta);
    
    % 計算其他參數
    hat_z = (z_positions(i)) / dl;
    hat_Ta = Ne * Af / (dl ^ 2);
    eta = dl * k / Su;
    hat_Z = Za / dl;
    delta_Ne = (last_Ne - Ne);
    delta_hat_Ta = delta_Ne * Af / (dl ^ 2);

    if i == 1
        delta_theta_a = ((En * Ncl / hat_Ta) - theta_a ^ 2 / (2 / (eta + hat_Z))) / (theta_a) * (hat_z - last_hat_z); 
    else
        delta_theta_a = -((En * Ncl / hat_Ta) - theta_a ^ 2 / (2 / (eta + hat_Z))) / (theta_a + ((delta_hat_Ta * (theta_a ^ 2)) / (hat_Ta * (theta_as - last_theta_as) * 2)) * (1 - (theta_s - last_theta_s) / (theta_a - last_theta_a))) * (hat_z - last_hat_z);
    end    
    
    
    fprintf('delta_theta_a = %f\n', delta_theta_a);
    
    last_Ne = Ne;
    last_theta_as = theta_as;
    last_theta_a = theta_a;
    last_theta_s = theta_s;
    last_hat_z = hat_z;

    % 更新參數
    theta_a = theta_a + delta_theta_a;
    theta_f = theta_f - delta_beta;
    theta_s = theta_s + delta_beta;
    theta_as = theta_a - theta_s;
    Za = Za + delta_z(i);      
    i = i + 1; 
end

% 輸出結果
if i > 10000
    disp('達到最大迭代次數');
elseif theta_f <= 0
    fprintf('theta_f = %f\n',theta_f);
    disp('計算完成：theta_f 達到或小於 0');
end

fprintf('迭代次數: %d\n', i-1);
if i > 1
    fprintf('最終 x 位置: %.4f\n', x_positions(i-1));
    fprintf('最終 z 位置: %.4f\n', z_positions(i-1));
else
    fprintf('沒有有效的計算結果\n');
end

% 繪圖
if i > 1
    figure;
    hold on;
    plot(x_positions(1:i-1), -z_positions(1:i-1), 'b-', 'LineWidth', 2);
    xlabel('拖曳距離 (Drag Distance)');
    ylabel('貫入深度 (Penetration Depth)');
    title('嵌入式拖錨 DEA 貫入軌跡示意圖');
    set(gca, 'XAxisLocation', 'top');
    grid on;
    max_x = max(x_positions(1:i-1));
    max_z = max(z_positions(1:i-1));
    xlim([0, max(max_x, 1)]);
    ylim([0, max(max_z, 1)]);
    hold off;
else
    disp('沒有足夠的數據來繪圖');
end