% clear workspaces
clear
clc

% define variables
xmin =0 ;
xmax = 1;
N = 100;   % no. nodes -1
dt = 0.009;  % time step
t = 0;
tmax = 0.5;
v = 0.9; % velocity

% discretize the domain
dx = (xmax-xmin)/N;
x = xmin - dx : dx : xmax + dx; % we need ghost nodes

% set intitial condition
u0 = sin(2*pi*x); %u0 = exp(-200*(x - 0.25).^2);
plot(x,u0); % check the init cond.
u= u0;
unp1 = u0;

% loop through time
nsteps = tmax/dt;
for n=1 : nsteps
    %calculate boundary conditions
    u(1) = u(3);
    u(N+3) = u(N+1);
    for i = 2 : N+2
        unp1(i) = u(i) - v*dt/dx*(u(i) - u(i-1));
    end
    % update t and u
    t = t + dt;
    u = unp1;
    % calculate exact solution
    exact = sin(2*pi*(x-v*t));%exact = exp(-200*(x - 0.25 - v*t).^2);
    plot(x,exact,'r-');
    hold on
    plot(x,u, 'bo-','markerfacecolor','b');
    hold off
    axis([xmin xmax -0.5 1.5])
    ylim auto
    xlabel('x','fontsize', 16)
    ylabel('U(t,x)','fontsize',16)
    title(sprintf('time = %1.3f', t),'fontsize', 16)
    shg
    pause(dt);
end
