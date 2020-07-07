clear;

% Parameters
m = 0.5;
L = 0.25;
k = 3 * 10^(-6);
b = 1 * 10^(-7);
g = 9.81;

Ixx = 5 * 10^(-3);
Iyy = 5 * 10^(-3);
Izz = 1 * 10^(-2);
Tmax = 50;

% System Parameters
nx = 12;    % Nb of states
nu = 4;     % Nb of inputs
ny = 6;     % Nb of outputs

% Inputs at equilbrium point (all states 0)

u_eq = g*m/4;

%--------------------------------------------
% Linear state space model
%--------------------------------------------

% symbolic variables
syms x y z v_x v_y v_z phi theta psi w_x w_y w_z u1 u2 u3 u4

% state vector
state = [x; y; z; v_x; v_y; v_z; phi; theta; psi; w_x; w_y; w_z];

% input vector
input = [u1; u2; u3; u4];

% Rotation Matrix
R = [cos(psi)*cos(theta)-sin(phi)*sin(psi)*sin(theta)  -cos(phi)*sin(psi)  cos(psi)*sin(theta)+cos(theta)*sin(phi)*sin(psi);
     cos(theta)*sin(psi)+cos(psi)*sin(phi)*sin(theta)  cos(phi)*cos(psi)   sin(psi)*sin(theta)-cos(psi)*cos(theta)*sin(phi);
     -cos(phi)*sin(theta)                              sin(phi)            cos(phi)*cos(theta)];

% The functions f:

f1 = v_x;
f2 = v_y;
f3 = v_z;

lin_acc = (1/m)*([0; 0; -m*g] + R*[0; 0; u1]);

f4 = lin_acc(1);
f5 = lin_acc(2);
f6 = lin_acc(3);

R_pqr_ptp = [cos(theta)  0  -cos(phi)*sin(theta);
             0           1  sin(phi);
             sin(theta)  0  cos(phi)*cos(theta)];

eul_ang_vel = inv(R_pqr_ptp)*[w_x; w_y; w_z];         

f7 = eul_ang_vel(1);
f8 = eul_ang_vel(2);
f9 = eul_ang_vel(3);

I = [Ixx  0    0;
     0    Iyy  0;
     0    0    Izz];
 
ang_acc = inv(I)*([u2; u3; u4] -cross([w_x; w_y; w_z], I*[w_x; w_y; w_z])); 

f10 = ang_acc(1);
f11 = ang_acc(2);
f12 = ang_acc(3);

fun = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12];

% Deriving the functions in the state variables (Jacobian)
J = jacobian(fun, state);
 
% Evaluating the jacobian in the equilibrium values: the result is A
A = subs(J,[state; input],[zeros(nx,1); [u_eq; 0; 0; 0]]);
A = double(A)

% Deriving the functions in the input variables
J = jacobian(fun, input);

% Evaluating the derivatives in the equilibrium values: the result is B
B = subs(J, [state; input],[zeros(nx,1); [u_eq; 0; 0; 0]]);
B = double(B)

% The output consists of states 1 to 3 and 7 to 9, so C selects these and D
% is zero
C = [eye(3), zeros(3,9);
     zeros(3,6), eye(3), zeros(3,3)];
D = zeros(ny,nu);

% Creating the continuous time system
c_sys = ss(A,B,C,D);

%-----------------------------------------
% stability
%-----------------------------------------

% Output poles
disp('Poles:')
disp(eig(A))

% Plotting the location of the poles
% pzmap(c_sys);

% Step response
%figure;
%step(c_sys,Tmax)

% Impulse response
%figure;
%impulse(c_sys,Tmax)

% Rank of the matrix

disp('Rank of a matrix');

rank(A)

%-----------------------------------------
% controllablity
%-----------------------------------------

disp ('Controllability matrix');

CO = ctrb(A,B);

disp('Rank of the controllability matrix:');
rank(CO)

%-----------------------------------------
% observablity
%-----------------------------------------

disp ('Observability matrix');

OB = obsv(A,C);

disp('Rank of the observability matrix:');
rank(OB)



