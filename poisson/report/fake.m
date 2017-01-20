f = @(x, y) 1 + sin(x.*y)
g = @(x, y) 1 + sin(x.*y) + 0.001.*(1 - sqrt(x.^2 + y.^2))

x = -2 : 0.08 : 2
y = -2 : 0.08 : 2

[X, Y] = meshgrid(x, y)

hold on
surf(X, Y, f(X, Y), 'FaceColor', 'red', 'FaceAlpha', 0.5)
surf(X, Y, f(X, Y), 'FaceColor', 'green', 'FaceAlpha', 0.5)
