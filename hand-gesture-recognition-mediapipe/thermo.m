% Define the parameters
sigma = 10;
beta = 8/3;
rho = 28;

% Define the ordinary differential equations
ode = @(t, xyz) [sigma * (xyz(2) - xyz(1)); 
                 rho * xyz(1) - xyz(2) - xyz(1) * xyz(3);
                 -beta * xyz(3) + xyz(1) * xyz(2)];

% Define the time span
t_span = [0 100];

% Define the initial conditions
initial_conditions_a = [0.1, 0.1, 0.1];
initial_conditions_b = [0.100001, 0.1, 0.1];

% Solve the system of differential equations for initial conditions a
[t_a, xyz_a] = ode45(ode, t_span, initial_conditions_a);

% Solve the system of differential equations for initial conditions b
[t_b, xyz_b] = ode45(ode, t_span, initial_conditions_b);

% Plot the results
figure;
subplot(3,1,1);
plot(t_a, xyz_a(:,1), 'b', t_b, xyz_b(:,1), 'r');
xlabel('Time');
ylabel('x(t)');
legend('Initial conditions a', 'Initial conditions b');

subplot(3,1,2);
plot(t_a, xyz_a(:,2), 'b', t_b, xyz_b(:,2), 'r');
xlabel('Time');
ylabel('y(t)');
legend('Initial conditions a', 'Initial conditions b');

subplot(3,1,3);
plot(t_a, xyz_a(:,3), 'b', t_b, xyz_b(:,3), 'r');
xlabel('Time');
ylabel('z(t)');
legend('Initial conditions a', 'Initial conditions b');

sgtitle('Variations of x(t), y(t), and z(t) from t = 0 to t = 100');

% Find the index corresponding to t = 100
index_100_a = find(t_a == 100);
index_100_b = find(t_b == 100);

% Extract the final values of x(t), y(t), and z(t) for initial conditions a
final_values_a = xyz_a(index_100_a, :);

% Extract the final values of x(t), y(t), and z(t) for initial conditions b
final_values_b = xyz_b(index_100_b, :);