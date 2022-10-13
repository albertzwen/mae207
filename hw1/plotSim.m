function fig = plotSim(pos, title_var, subtitle_var)
cla; clf; close all;
fig = figure;
hold on;
grid on;
view(-45, 30);
xlabel("x_1");
ylabel("x_2");
zlabel("x_3");
title(title_var);
subtitle(subtitle_var);
plot3(pos(:, 1), pos(:, 2), pos(:, 3), 'MarkerSize', 5);

