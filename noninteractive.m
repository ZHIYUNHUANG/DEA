% 參數設定
tf = 0.28;     % m
Lf = 2;        % m
Lj = 0.5;      % m
Ls = 3;        % m
theta_fs = 50; % 度
alphas = [1, 1/2, 1/3];  % alpha 的值
theta_as = linspace(-30, 30);  % degrees

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
    if alpha == 1
        figure;  % 創建一個新圖窗口
        % 繪製原始曲線
        plot(theta_as, Nem, 'b--', 'LineWidth', 2);  % 繪製第一條線，藍色
        hold on;  % 保持圖形
        plot(theta_as, Nen, 'r--', 'LineWidth', 2);  % 繪製第二條線，紅色虛線
        hold on;
        plot(theta_as, Net, '--', 'LineWidth', 2);  % 繪製第三條線，紫色虛線
    end
    % 找到每個 x 對應的 y 值的最小值
    min_values = min([Nem; Nen; Net]);

    % 將最小值點連線成黑色線
    plot(theta_as, min_values, '-', 'Color', [0 0 0], 'LineWidth', 2); 
end

xlabel('Force Angle');
ylabel('Bearing Factor');
title('Anchor capacity curve of noninteractive resistance model');
ylim([0 16]); 
legend('Nem', 'Nen', 'Net', 'Minimum Values');  
