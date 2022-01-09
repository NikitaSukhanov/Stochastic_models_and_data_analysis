clear all
pkg load optim
addpath(genpath(fullfile('octave-interval-examples', 'm')));
load data.mat

# plot initial dataset
figure
hold on
set(gca, 'fontsize', 22, "linewidth", 2)
plot(data_x, data_y, "linewidth", 2)
xlabel('n');
ylabel('cond_{med}(A_n)');

left = 1;
right = 20;
n = right - left + 1
x = data_x(left:right)';
y = data_y(left:right)';
X = [ones([n, 1]) x];

# plot chosen fragment
figure
hold on
set(gca, 'fontsize', 22, "linewidth", 2)
ylim ([0, 6]);
scatter(x, y, "filled")
xlabel('n');
ylabel('cond_{med}(A_n)');

rad_y = ones([n, 1]) * 0.1;
irp_initial = ir_problem(X, y, rad_y);

# plot LSM regression
lsm_betta = (X \ y)
x_min = min(x);
x_max = max(x);
figure
hold on
set(gca, 'fontsize', 22, "linewidth", 2)
ir_scatter(irp_initial, 'b.');
plot([x_min, x_max], [lsm_betta(1) + lsm_betta(2) * x_min, lsm_betta(1) + lsm_betta(2) * x_max], "linewidth", 2)
xlabel('n');
ylabel('cond_{med}(A_n)');
legend({"interval data", "LSM linear regression"},  "location", "northwest")

# expand intervals
[linprog_betta, w] = lin_prog_symmetrical(X, y - rad_y, y + rad_y);
display(linprog_betta);
rad_y_expanded = rad_y .* w * 1.1;
irp_expanded = ir_problem(X, y, rad_y_expanded);
centroid = mean(ir_beta2poly(irp_expanded))
inf_y = y - rad_y_expanded;
sup_y = y + rad_y_expanded;

# plot infset
figure
hold on
set(gca, 'fontsize', 22, "linewidth", 2)
ir_plotbeta(irp_expanded)
xlim ([0.95, 1.04]);
scatter(linprog_betta(1), linprog_betta(2), 60, "g", "filled")
scatter(lsm_betta(1), lsm_betta(2), 60, "r", "filled")
scatter(centroid(1), centroid(2), 60, "k", "filled")
xlabel('\beta_1')
ylabel('\beta_2')
legend({"inf. set", "inf. set boundary", "LP betta", "LSM betta", "centroid"}, "location", "east")

# find boundary points
predict = ir_predict(irp_expanded, irp_expanded.X);
predict_lower = predict(:, 1);
predict_upper = predict(:, 2);

eps = 1e-8;
indices = [];
for i = 1:n
  if abs(predict_lower(i) - inf_y(i)) < eps
    indices = [indices, i];
    continue
  endif
  if abs(predict_lower(i) - sup_y(i)) < eps
    indices = [indices, i];
    continue
  endif
  if abs(predict_upper(i) - inf_y(i)) < eps
    indices = [indices, i];
    continue
  endif
  if abs(predict_upper(i) - sup_y(i)) < eps
    indices = [indices, i];
    continue
  endif
end
display(indices)
other = setdiff(1:n, indices);
irp_boundary = ir_problem(X(indices, :), y(indices), rad_y_expanded(indices));
irp_other = ir_problem(X(other, :), y(other), rad_y_expanded(other));

# plot compatibility cone
x_min = 0;
x_max = 29;
figure
hold on
set(gca, 'fontsize', 22, "linewidth", 2)
ir_plotmodelset(irp_expanded, [x_min, x_max])
ir_scatter(irp_boundary, 'ro');
ir_scatter(irp_other, 'b.');
xlabel('n');
ylabel('cond_{med}(A_n)');
legend({"compatibility cone", "", "", "boundary points", "other points"}, "location", "northwest")

# plot compatibility cone close
figure
hold on
set(gca, 'fontsize', 22, "linewidth", 2)
xlim ([0, 2]);
ylim([0.8, 1.3]);
ir_plotmodelset(irp_expanded, [x_min, x_max])
plot([x_min, x_max], [linprog_betta(1) + linprog_betta(2) * x_min, linprog_betta(1) + linprog_betta(2) * x_max], "g", "linewidth", 2)
plot([x_min, x_max], [lsm_betta(1) + lsm_betta(2) * x_min, lsm_betta(1) + lsm_betta(2) * x_max], "r", "linewidth", 2)
plot([x_min, x_max], [centroid(1) + centroid(2) * x_min, centroid(1) + centroid(2) * x_max], "k--", "linewidth", 2)
ir_scatter(irp_boundary, 'ro');
ir_scatter(irp_other, 'b.');
xlabel('n');
ylabel('cond_{med}(A_n)');
legend({"compatibility cone", "", "", "LP betta", "LSM betta", "centroid", "(x_1, y_1)"}, "location", "northwest")


test_n = 10;
x_test = data_x(right + 1: right + test_n)';
y_test = data_y(right + 1: right + test_n)';
X_test = [ones([test_n, 1]) x_test];
y_predict = ir_predict(irp_expanded, X_test)
mid_y_predict = (y_predict(:, 1) + y_predict(:, 2)) / 2;
rad_y_predict = (y_predict(:, 2) - y_predict(:, 1)) / 2
irp_predict = ir_problem(X_test, mid_y_predict, rad_y_predict);

figure
hold on
set(gca, 'fontsize', 22, "linewidth", 2)
xlim ([27, 50]);
ylim([3, 9]);
ir_plotmodelset(irp_expanded, [27, 50])
ir_scatter(irp_predict, 'b.')
scatter(x_test, y_test, 'r', "filled")
xlabel('n');
ylabel('cond_{med}(A_n)');
legend({"compatibility cone", "", "", "predicted points", "data points"}, "location", "northwest")






