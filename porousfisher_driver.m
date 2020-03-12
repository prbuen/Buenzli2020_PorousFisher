% Nondimensionalised Porous-Fisher model applied to tissue engineering pore infilling.
%
% Authors: Pascal R Buenzli, Matthew J Simpson, 2019-2020

clear all

alpha2 = 0.02; % alpha^2: 0.01, 0.02, 0.05, 0.1, 0.2
u_initial = 0.5; % 0.2, 0.5, 0.8
u_c = 0.5; % 0.01 (leading edge); 0.2, 0.5, 0.8. NB: influences calculation of pore_area

show_plots = 1; % run twice to make last plot same size as the others (bug?)
show_plots = 0; % to calculate Tbridge

L = 1;
dx = L/50; % for speed
%dx = L/150; % for precision
N = round(L/dx+1);
u0 = zeros(N,N);

for i=1:N
    u0(i,1)=u_initial;
    u0(i,N)=u_initial;
    u0(1,i)=u_initial;
    u0(N,i)=u_initial;
end


u0=reshape(u0,[],1);

if show_plots
    tspan = 0:1:6;
else    
    tspan = [0 6];
end
[t,u] = ode45(@(t,u) porousfisher_odefn(t,u,N,alpha2,dx), tspan, u0);
% [t,u] = ode15s(@(t,u) porousfisher_odefn(t,u,N,alpha,dx), tspan, u0);

x=0:dx:L;
y=0:dx:L;

u = reshape(u, [], N,N); % to have u(t, x, y)
size_u = size(u)

if show_plots
    % plot still frames (x,y) -> u(t,x,y)
    x=0:dx:L;
    y=0:dx:L;
    for n = 1:length(tspan)
        % axis equal
        axis image
        subplot(1,length(tspan),n)
        hold on
        colormap gray
        contourf(x,y,reshape(u(n,:,:),N,N),0:0.01:1, 'LineColor', 'none');
        contour(x,y,reshape(u(n,:,:),N,N), 0:0.1:1, 'LineColor', 'black');
        contour(x,y,reshape(u(n,:,:),N,N),[u_c, u_c], 'r', 'LineWidth', 2);
        caxis([0 1]);
        colorbar
        hold off
        title(sprintf('t=%d', t(n)))
    end
else
    % determine Tbridge
    Tbridge = -1;
    for n=1:size_u(1)
        pore_area = 0;
        for i=1:N
            for j=1:N
                if u(n,i,j) < u_c
                    pore_area = pore_area + 1;
                end
            end
        end
        pore_area = pore_area * dx^2;
        if pore_area < 1e-10
            Tbridge = t(n)
            break
        end
    end
end