function a = checkConvergenceRate(f,action,vecField,z0,L,m)

    T = 0.1;
    Nrange = 25:5:40; 
    errorLE = zeros(length(Nrange),1);
    errorRK4 = errorLE;
    errorEH = errorLE;
    count = 1;
    
    q0 = extractq(z0);
    w0 = extractw(z0);
    
    z0Ref = [q0;w0]; %We reorder the initial condition to make ode45 work withour reordering the equations
    %i.e. we pass from [q_1,w_1,...,q_N,w_N] to
    %[q_1,q_2,...,q_N,w_1,w_2,...,w_N]
    Nref = 1000;
    time = linspace(0,T,Nref);
    odeFunc = @(t,z) [FuncQ(z);FuncW(z,L,m)]; %FuncQ,FuncW just assemble the first N and last N equations
    [timeSol,zSol] = ode45(odeFunc, time, z0Ref);
    zSol = zSol';
    zRef = zSol(:,end);
    zRef = reorder(zRef); %we pass from [q_1,q_2,...,q_N,w_1,w_2,...,w_N] to
    %[q_1,w_1,...,q_N,w_N] so that we can compare it easily with the
    %solution coming from our numerica schemes
    
    for N = Nrange
        
        time = linspace(0,T,N);
        dt = time(2)-time(1);
        z1 = z0;
        z2 = z0;
        z3 = z0;
        
        for i = 1:N-1
            z1 = LieEulerSE3N(vecField,action,z1,dt);
            z2 = EulerHeunSE3N(vecField,action,z2,dt);
            z3 = FreeRK4SE3N(f,action,dt,z3);
        end
        errorLE(count) = norm(z1-zRef,inf);
        errorEH(count) = norm(z2-zRef,inf);
        errorRK4(count) = norm(z3-zRef,inf);
        count = count + 1;
    end
    
    %% Convergence rate of Lie Euler
    ord1 = (Nrange/Nrange(end)).^(-1)*errorLE(end); 
    ord2 = (Nrange/Nrange(end)).^(-2)*errorLE(end);
    ord3 = (Nrange/Nrange(end)).^(-3)*errorLE(end);
    ord4 = (Nrange/Nrange(end)).^(-4)*errorLE(end);

    figure;
    loglog(Nrange,errorLE,'r-*',Nrange,ord1,'k-',Nrange,ord2,'c-',Nrange,ord3,'r-',Nrange,ord4,'g-')
    title("Convergence rate Lie Euler against the reference ODE45 solution")
    legend('Error','ord1','ord2','ord3','ord4');

    %% Convergence rate of Euler Heun
    ord1 = (Nrange/Nrange(end)).^(-1)*errorEH(end); 
    ord2 = (Nrange/Nrange(end)).^(-2)*errorEH(end);
    ord3 = (Nrange/Nrange(end)).^(-3)*errorEH(end);
    ord4 = (Nrange/Nrange(end)).^(-4)*errorEH(end);

    figure;
    loglog(Nrange,errorEH,'r-*',Nrange,ord1,'k-',Nrange,ord2,'c-',Nrange,ord3,'r-',Nrange,ord4,'g-')
    title("Convergence rate EULER HEUN against the reference ODE45 solution")
    legend('Error','ord1','ord2','ord3','ord4');


    %% Convergence rate of RKMK4
    ord1 = (Nrange/Nrange(end)).^(-1)*errorRK4(end);
    ord2 = (Nrange/Nrange(end)).^(-2)*errorRK4(end);
    ord3 = (Nrange/Nrange(end)).^(-3)*errorRK4(end);
    ord4 = (Nrange/Nrange(end)).^(-4)*errorRK4(end);

    figure;
    loglog(Nrange,errorRK4,'r-*',Nrange,ord1,'k-',Nrange,ord2,'c-',Nrange,ord3,'r-',Nrange,ord4,'g-')
    title("Convergence rate RKMK4 against the reference ODE45 solution")
    legend('Error','ord1','ord2','ord3','ord4');

end