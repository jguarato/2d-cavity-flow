% 2D lid-driven cavity flow: incompressible Navier-Stokes equations
% Solved with finite volume method and the fractional step method 
% Developed by: Jessica Guarato
% Last modified on: November, 2017

clc; clear all; close all;

%% Setup

% Cavity                                                % Units
Lx              = 1;                                    % [m]
Ly              = 1;                                    % [m]
U               = 1;                                    % [m/s]

% Fluid properties
rho             = 1;                                    % [kg/m^3]
mu              = 1;                                    % [Pa*s]

% Number of volumes
Nx              = 100;
Ny              = 100;

% Time step
dt              = 1e-5;
tfinal          = 0.5;
ctfinal         = 999999;

% Gravity
gx              = 0;
gy              = -9.81;

%% Solving problem

nx = Nx; ny = Ny;
dx = Lx/nx;
dy = Ly/ny;

ngx = nx + 2;
ngy = ny + 2;

nu = mu/rho;

% Initial conditions
t = 0;
ct = 1;

u = zeros(ngy+1,ngx+1);
v = zeros(ngy+1,ngx+1);
p = zeros(ngy,ngx);

for j = 1:ngy+1
    for i = 1:ngx+1

        if (j == ngy+1 || j == ngy) % Sliding lid 
            u(j,i) = U;
            v(j,i) = 0;
        end

    end
end

u1 = u;
v1 = v;
p1 = p;

residue = 1;
while (ct <= ctfinal && t <= tfinal && residue > 1e-6)
    
    tstart = tic;
    
    u_pred = zeros(ngy,ngx);
    v_pred = zeros(ngy,ngx);
    
for j = 2:ngy
    for i = 2:ngx
        
        % Updating values
        p(j,i) = p1(j,i);
        u(j,i) = u1(j,i);
        v(j,i) = v1(j,i);

        % Velocities at faces
        up = u(j,i);
        ue = 0.5*(u(j,i+1) + u(j,i));
        uw = 0.5*(u(j,i-1) + u(j,i));
        un = 0.5*(u(j+1,i) + u(j,i));
        us = 0.5*(u(j-1,i) + u(j,i));

        vp = v(j,i);
        ve = 0.5*(v(j,i+1) + v(j,i));
        vw = 0.5*(v(j,i-1) + v(j,i));
        vn = 0.5*(v(j+1,i) + v(j,i));
        vs = 0.5*(v(j-1,i) + v(j,i));
        
        % Derivatives at faces
        dudx_e = (u(j,i+1) - u(j,i))/dx;
        dudx_w = (u(j,i) - u(j,i-1))/dx;
        dudy_n = (u(j+1,i) - u(j,i))/dy;
        dudy_s = (u(j,i) - u(j-1,i))/dy;
        
        dvdx_e = (v(j,i+1) - v(j,i))/dx;
        dvdx_w = (v(j,i) - v(j,i-1))/dx;
        dvdy_n = (v(j+1,i) - v(j,i))/dy;
        dvdy_s = (v(j,i) - v(j-1,i))/dy;
        
        % Pressure at faces
        pe = p(j,i);
        pw = p(j,i-1);
        pn = p(j,i);
        ps = p(j-1,i);

        % Predictor step for u
        u_pred(j,i) = up + ( -1/rho*(pe - pw)/dx + gx ...
                      + nu*((dudx_e - dudx_w)/dx + (dudy_n - dudy_s)/dy) ...
                      - (ue*ue - uw*uw)/dx - (un*vn - us*vs)/dy )*dt;
        
        % Predictor step for v
        v_pred(j,i) = vp + ( -1/rho*(pn - ps)/dy + gy ...
                      + nu*((dvdx_e - dvdx_w)/dx + (dvdy_n - dvdy_s)/dy) ...
                      - (ve*ue - vw*uw)/dx -(vn*vn - vs*vs)/dy )*dt;
               
    end
end

for j = 1:ngy+1
    for i = 1:ngx+1
        
        % Boundary conditions for velocities
        if (j == ngy+1 || j == ngy) % Sliding lid 
            u_pred(j,i) = U;
            v_pred(j,i) = 0;
            
        elseif (i == 1 || i == 2 || i == ngx+1 || i == ngx || j == 1 || j == 2) % Wall condition
            u_pred(j,i) = 0;
            v_pred(j,i) = 0;
            
        end
        
    end
end

Ap = sparse(zeros(ngx*ngy,ngx*ngy));
Bp = sparse(zeros(ngx*ngy,1));

for j = 2:ngy-1
    for i = 2:ngx-1
    
        m = (j-1)*ngx + i;

        % Boundary conditions for pressure
        % Southwest edge:
        if (m == 1)

            Ap(m,m) = -(dy/dx + dx/dy);
            Ap(m,m+1) = dy/dx;
            Ap(m,m+ngx) = dx/dy;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;
                        
        % Southeast edge:
        elseif (m == ngx)

            Ap(m,m) = -(dy/dx + dx/dy);
            Ap(m,m-1) = dy/dx;
            Ap(m,m+ngx) = dx/dy;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;               
        
        % Northwest edge:
        elseif m == (ngx*ngy - ngx + 1)

            Ap(m,m) = -(dy/dx + dx/dy);
            Ap(m,m-ngx) = dx/dy;
            Ap(m,m+1) = dy/dx;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;                 
        
        % Northeast edge:
        elseif m == ngx*ngy
            
            Ap(m,m) = -(dy/dx + dx/dy);
            Ap(m,m-ngx) = dx/dy;
            Ap(m,m-1) = dy/dx;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;  
                        
        % West boundary:
        elseif (m ~= 1 && m ~= (ngx*ngy - ngx + 1) && rem((m - 1),ngx) == 0)
            
            Ap(m,m) = -(dy/dx + 2*dx/dy);
            Ap(m,m-ngx) = dx/dy;
            Ap(m,m+1) = dy/dx;
            Ap(m,m+ngx) = dx/dy;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;

        % East boundary:
        elseif (m ~= ngx && m ~= ngx*ngy && rem(m,ngx) == 0)

            Ap(m,m) = -(dy/dx + 2*dx/dy);
            Ap(m,m-ngx) = dx/dy;
            Ap(m,m-1) = dy/dx;
            Ap(m,m+ngx) = dx/dy;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;

        % South boundary:
        elseif (m > 1 && m < ngx)

            Ap(m,m) = -(2*dy/dx + dx/dy);
            Ap(m,m-1) = dy/dx;
            Ap(m,m+1) = dy/dx;
            Ap(m,m+ngx) = dx/dy;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;

        % North boundary:
        elseif (m > (ngx*ngy - ngx + 1) && m < ngx*ngy)

            Ap(m,m) = -(2*dy/dx + dx/dy);
            Ap(m,m-ngx) = dx/dy;
            Ap(m,m-1) = dy/dx;
            Ap(m,m+1) = dy/dx;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;
        
        % Inside
        else

            Ap(m,m) = -(2*dy/dx + 2*dx/dy);
            Ap(m,m-ngx) = dx/dy;
            Ap(m,m-1) = dy/dx;
            Ap(m,m+1) = dy/dx;
            Ap(m,m+ngx) = dx/dy;
            Bp(m,1) = rho*((u_pred(j,i+1) - u_pred(j,i))*dy + ...
                            (v_pred(j+1,i) - v_pred(j,i))*dx)/dt;

        end
    end 
end

for m = 1:(ngx*ngy)
    
    if Ap(m,m) == 0
        Ap(m,m) = 1;
    end
    
end

% Solving the linear system
p_pred_m = Ap\Bp;

for j = 2:ngy
    for i = 2:ngx
        
        % Rearranging in matrix form
        p_pred(j,i) = p_pred_m((j-1)*ngx + i);

        pe = p_pred(j,i);
        pw = p_pred(j,i-1);
        pn = p_pred(j,i);
        ps = p_pred(j-1,i);
        
        % Obtaining new values
        p1(j,i) = p(j,i) + p_pred(j,i);
        u1(j,i) = u_pred(j,i) - 1/rho*(pe - pw)*dt/dx;
        v1(j,i) = v_pred(j,i) - 1/rho*(pn - ps)*dt/dy;
        
    end
end

for j = 1:ngy+1
    for i = 1:ngx+1
        
        % Boundary conditions: velocities
        if (j == ngy+1 || j == ngy) % Sliding lid 
            u1(j,i) = U;
            v1(j,i) = 0;
            
        elseif (i == 1 || i == 2 || i == ngx+1 || i == ngx || j == 1 || j == 2) % Wall condition
            u1(j,i) = 0;
            v1(j,i) = 0;

        end
        
    end
end

% Calculating the vorticity
w = zeros(ngy,ngx);
for j = 2:ngy-1
    for i = 2:ngx-1
        
        un = (u(j+1,i+1) + u(j+1,i) + u(j,i+1) + u(j,i))*0.25;
        us = (u(j,i+1) + u(j,i) + u(j-1,i+1) + u(j-1,i))*0.25;
        ve = (v(j+1,i+1) + v(j,i+1) + v(j+1,i) + v(j,i))*0.25;
        vw = (v(j+1,i) + v(j,i) + v(j+1,i-1) + v(j,i-1))*0.25;
       
        dudy = (un - us)/dy;
        dvdx = (ve - vw)/dx;
           
        w(j,i) = dvdx - dudy;
      
    end
end

    residue = max(max(abs(u1 - u)));
    telapsed = toc(tstart);
    clc;
    fprintf('\n ============================');
    fprintf('\n Iter = %06d',ct);
    fprintf('\n t = %05e',t);
    fprintf('\n residue max = %05e',residue);
    fprintf('\n elapsed time = %05e',telapsed);
    fprintf('\n ');
    
    t = t + dt;
    ct = ct + 1;
    %keyboard;
end

%% Results

x = 0:dx:Lx;
y = 0:dy:Ly;

% Mesh
figure
[X,Y] = meshgrid(x,y); surf(X,Y,ones(ny+1,nx+1));
view(2); axis equal; axis ([0 Lx 0 Ly]);
xlabel('x [m]'); ylabel('y [m]');

figure
%contourf(x,y,u);
imagesc(x,y,u1(2:end-1,2:end-1)); set(gca,'YDir','normal');
colormap(jet); colorbar;
axis equal; axis ([0 Lx 0 Ly]);
xlabel('x [m]'); ylabel('y [m]'); title('Velocity u');

figure
%contourf(x,y,u);
imagesc(x,y,v1(2:end-1,2:end-1)); set(gca,'YDir','normal');
colormap(jet); colorbar;
axis equal; axis ([0 Lx 0 Ly]);
xlabel('x [m]'); ylabel('y [m]'); title('Velocity v');

figure
%contourf(x,y,u);
imagesc(x,y,w(2:end-1,2:end-1)); set(gca,'YDir','normal');
colormap(jet); colorbar;
axis equal; axis ([0 Lx 0 Ly]);
xlabel('x [m]'); ylabel('y [m]'); title('Vorticity');

figure
%contourf(x,y,u);
imagesc(x,y,p1(2:end-1,2:end-1)); set(gca,'YDir','normal');
colormap(jet); colorbar;
axis equal; axis ([0 Lx 0 Ly]);
xlabel('x [m]'); ylabel('y [m]'); title('Pressure');

