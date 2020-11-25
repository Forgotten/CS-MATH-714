% wave propagation

% Grid and initial data:
  N = 64; 
  x = cos(pi*(0:N)/N); 
  y = x';
  % final time
  T = 1;
  dt = 6/N^2;
  [xx,yy] = meshgrid(x,y);
  % number of iterations
  Nit = round(T/dt)
  
  plot_bool = true;
  plot_n    = 100;
  
  % initial condition
  ut = @(x,y)(exp(-400*(x).^2).*exp(-400*(y).^2));

  % evaluating the initial condition
  dvvdt = ut(xx,yy);
  
  % using the initial condition to compute u^1
  vv = dt*dvvdt + 1/6*dt^3*laplacian(dvvdt,x,y);
  % setting u^0 to zero
  vvold = zeros(N+1, N+1); 

% Time-stepping by leap frog formula:
  [ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
  
  for n = 0:Nit
    t = n*dt;
    % plot if necessary
    if plot_bool && (mod(n,plot_n) == 0 )
        figure(1);
        surf(xx, yy, vv);
        axis([-1 1  -1 1  -0.05 0.05]);
        shading interp
        drawnow;
        fprintf('iteration number : %i \n', n)
        
    end
    
    lapu = laplacian_opt(vv,x,y);
    % we could compute both simultaneously
    lap2u = laplacian_opt(lapu,x,y);
    vvnew = 2*vv - vvold + dt^2*(lapu) + 1/12*dt^4*lap2u; 
    vvold = vv; vv = vvnew;
  end
