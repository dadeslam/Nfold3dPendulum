clc;
clear all;
close all;

%% Notes about the code:
%The representation we use for the unkown vector and hence the order of the
%equations is (q_1,w_1,...,q_N,w_N), this just to reflect the structure of
%the manifold TS^2 x ... TS^2.

%The implemented numerical schemes are Lie Euler, Euler Heun and the
%Commutator Free RK4.

%
%% N connected 3D pendulums

prompt = 'How many 3D pendulums do you want to connect?\n\n';
P = input(prompt);

L = rand(P,1)+0.5; 
m = rand(P,1)+0.5;

t0 = 0;
T = 10;
N = 1000;
time = linspace(t0,T,N);
dt = time(2)-time(1);

getq = @(v) extractq(v); %(q_1,w_1,q_2,w_2,...,q_N,w_N)-->(q_1,q_2,...,q_N)
getw = @(v) extractw(v); %(q_1,w_1,q_2,w_2,...,q_N,w_N)-->(w_1,w_2,...,w_N)

f = @(v) fManiToAlgebra(getq(v),getw(v),L,m); %It computes a possible function from the manifold to the Lie algebra
%of the group acting transitively

action = @(B,input) actionSE3N(B,input); 
%Here in the action the matrix B has size 3*N x 4, because the
%representation of the element (A,a)€SE(3) we choose is the matrix [A a].
%While the vector input is of size 6*N x 1

vecField = @(sigma,p) dexpinvSE3N(sigma,f(action(exponentialSE3N(sigma),p))); %Vector field on the Lie algebra of
% (SE(3))^N

%% Solving the problem

[q0,w0,z0] = initializeSE3N(P); %It generates random admissible initial conditions
z = z0;
q = zeros(3*P,N); %matrix where we store the evolution of the state variable q
p = q; %here we store q with the lengths of the pendulums

Len = zeros(3*P,1);
for i = 1:P
    Len(3*i-2:3*i) = L(i)*ones(3,1);
end
Mat = diag(Len);
if P>1
    for i = 3:3:3*(P-1)
        Mat = Mat + diag(Len(1:3*P-i),-i);
    end
end
%The matrix Mat is a lower diagonal block matrix built to get Mat*q = p,
%just a transformation taking into account the lengths of the N pendulums

q(:,1) = q0;
p(:,i) = Mat*q0;

prompt = 'Do you want to see the convergence rate of the three methods? Write 1 for yes, 0 for no\n\n';
C = input(prompt);
if C==1
    checkConvergenceRate(f,action,vecField,z0,L,m) %In the function different time steps are tested and log-log
    %plots comparing them are shown. Here the tested methods are the ones
    %mentioned at the top of this script.
end

%% Time evolution of the solution

prompt = 'Do you want to see the Time Evolution of the solution? Write 1 for yes, 0 for no\n\n';
C1 = input(prompt);

if C1==1
    z = z0;
    for i = 1:N-1
        z = FreeRK4SE3N(f,action,dt,z);
        q(:,i+1) = extractq(z);
        p(:,i+1) = Mat*q(:,i+1);
    end
    
    figure('Units','normalized','Position',[0 0 1 1])
    for i = 1:2:N

        plot3([0,p(1,i)],[0,p(2,i)],[0,p(3,i)],'r-*',...
            [p(3*(1:P-1)-2,i),p(3*(1:P-1)+1,i)],[p(3*(1:P-1)-1,i),p(3*(1:P-1)+2,i)],...
            [p(3*(1:P-1),i),p(3*(1:P-1)+3,i)],'k-o','Markersize',4);

        str = "Time evolution of the pendulum, current time t = "+string(time(i));
        axis([-sum(L) sum(L) -sum(L) sum(L) -sum(L) sum(L)]);
        title(str)
        pause(0.00000000000001);
    end
end