clc;
%clear all;
%% Inputs

Ra = 8000; Pr = 0.71;
Ly = 20;
Lx = 1;
Nx = (Lx*40)+1; 
Ny = (Ly*40)+1;		% No. of grids along X and Y axes
bita = sqrt(0.6);

% Initialization of N.D. velocities, Pressure and Temperature matrices
U = zeros(Nx,Ny);
V = zeros(Nx,Ny);
P = zeros(Nx,Ny);
T = zeros(Nx,Ny); 		% Dimensionless temperarure
T_temp = zeros(Nx,Ny);
Nu = zeros(Ny);

%Setting initial values / conditions
P(:,:) = 0;
T_temp(:,:)= 0;
T(2:Nx-1,2:Ny-1) = 0;

%Boundary Conditions
U(1,:) = 0; V(1,:) = 0; 	% BC at Left verticle side
U(Nx,:) = 0; V(Nx,:) = 0;	% BC at Right verticle side
U(:,Ny) = 0; V(:,Ny) = 0; 	% BC at Top horizontal side
U(:,1) = 0; V(:,1) = 0; 	% BC at Bottom horizontal side
T(1,:) = -0.5;
T(Nx,:) = 0.5;
%% Eval of Inputs
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dt = 0.005;

% grid spacing along X and Y axes and time step
x = 0:dx:Lx;
y = 0:dy:Ly;

% setting initial values
iter = 1; 			% initial iteration value
error = 1; 			% inital error value for each iteration
error1 = zeros(1,2); 		% initializing 1D array for storing all error values
tic; 				% timer starts

while(error>1e-8 && iter<40000)
	% progm should exit loop if any of two conditions violated
	for j = 2:Ny-1; 	% Loop for y
		for i = 2:Nx-1;		 % Loop for x
P(i,j) = P(i,j)-((0.5*dt)/(dx*bita^2))*(U(i+1,j)-U(i-1,j))-(0.5*dt/(dy*bita^2))*(V(i,j+1 )-V(i,j-1));

U(i,j ) = U(i,j)-(0.5*dt/dx)*(U(i+1,j)^2-U(i-1,j)^2)-(0.5*dt/dy)*(U(i,j+1)*V(i,j+1)-U(i,j-1)*V(i,j-1))-(0.5*dt/dx)*(P(i+1,j)-P(i-1,j))+(sqrt(Pr/Ra)*(dt/dx^2))*(U(i+1,j)-2*U(i,j)+U(i-1,j))+(sqrt(Pr/Ra)*(dt/dy^2))*(U(i,j+1)-2*U(i,j)+U(i,j-1));

V(i,j)= V(i,j)-(0.5*dt/dy)*(V(i,j+1)^2-V(i,j-1)^2)-(0.5*dt/dx)*(U(i+1,j)*V(i+1,j)-U(i-1,j)*V(i-1,j))-(0.5*dt/dy)*(P(i,j+1)-P(i,j-1))+(sqrt(Pr/Ra)*(dt/dy^2))*(V(i,j+1)-2*V(i,j)+V(i,j-1))+(sqrt(Pr/Ra)*(dt/dx^2))*(V(i+1,j)-2*V(i,j)+V(i-1,j))+ dt*T(i,j);

T(i,j)= T(i,j)-(0.5*dt/dx)*(T(i+1,j)*U(i+1,j)-T(i-1,j)*U(i-1,j))-(0.5*dt/dy)*(T(i,j+1)*V(i,j+1)-T(i,j-1)*V(i,j-1))+(dt/(sqrt(Ra*Pr)*dx^2))*(T(i+1,j)-2*T(i,j)+T(i-1,j))+(dt/(sqrt(Ra*Pr)*dy^2))*(T(i,j+1)-2*T(i,j)+T(i,j-1));
		end
	end
	error = max(max(abs(T-T_temp)));
	error1(1,iter) = error;
	%% Neuman Boundary Conditions
	T(:,1) = T(:,2);
	T(:,Ny) = T(:,Ny-1);
	% pressure at boundary is equal to pressure at adjacent node
	P(:,1) = P(:,2);
	P(:,Ny) = P(:,Ny-1);
	P(1,:) = P(2,:);
	P(Nx,:) = P(Nx-1,:);
	T_temp = T ;
	iter = iter + 1;
end
toc 				% CPU time

Wz = zeros(Nx,Ny); 		% initializing vorticity matrix
for i = 2:Nx-1; 		% for inner grids along X-axis
	for j = 2:Ny-1; 	% for inner grids along Y-axis
		% vorticity
		Wz(i,j) = (0.5*(V(i+1,j)-V(i-1,j))/dx)-(0.5*(U(i,j+1)-U(i,j-1))/dy);
	end
end

for j=1:Ny
	Nu(j) = Nx*(T(Nx,2)-T(Nx-1,j))/(T(Nx,2)-T(1,2));
end
%% Analysis of Outputs or Post Processing
quiver(U',V') 			% velocity field
title({'AR = 20 & Ra = 8000','Velocity Vector Plot'})
pbaspect([6 12 1])
% pbaspect([6 8 1])
ax=gca;
ax.FontSize=20;
set(gca,'xtick',[])
set(gca,'ytick',[])

% figure(2);
% contourf(U',40) 		% filled contours of U velocity
% title({'AR = 20 & Ra = 8000','U velocity Contours'})
% pbaspect([6 8 1])
% ax=gca;
% ax.FontSize=20;

figure(3);
contourf(V',40) 		% filled contours of V velocity
title({'AR = 20 & Ra = 8000','V velocity Contours'})
pbaspect([6 12 1])
% pbaspect([6 8 1])
ax=gca;
ax.FontSize=20;
set(gca,'xtick',[])
set(gca,'ytick',[])
% 
% figure(4);
% contourf(T',40) 		% filled contours of Temperature profile
% title({'AR = 20 & Ra = 6500','Temperature Contours'})
% % pbaspect([6 12 1])
% pbaspect([6 8 1])
% ax=gca;
% ax.FontSize=20;

% figure(5);
% contourf(Wz',40) 		% filled contours of Vorticity
% title('Vorticity Contour Graph')
% % pbaspect([6 12 1])
% 
Uy = U((Nx-1)/2,:); 		% U velocity on Y-axis at horizontal centre
figure(6);
plot(Uy', y/Ly)
title({'AR = 20 & Ra = 8000','U velocity profile at Horizontal Centre'});
xlabel('U velocity');
ylabel('Length along Y-axis');
pbaspect([6 12 1])
ax=gca;
ax.FontSize=20;
set(gca,'xlabel',[])
set(gca,'ylabel',[])
% 
% Vx=V(:,(Ny-1)/2); 		% V velocity on X-axis at verticle centre
% figure(7);
% plot(x,Vx')
% title('V velocity profile at Verticle Centre');
% xlabel('Length along X-axis');
% ylabel('V velocity');
% 
% figure(8);
% semilogy(error1)
% title('Error Graph');
% 
% figure(9);
% plot(y/Ly,T(21,:))
% title({'AR = 20 & Ra = 6500','Temperature Profile Along Vertical Axis at Horizontal Mid'});
% xlabel('Dimensionless Height');
% ylabel('Dimensionless Temperature');
% % pbaspect([12 6 1])
% pbaspect([6 8 1])
% ax=gca;
% ax.FontSize=20;
% 
figure(10);
plot(y/Ly,Nu(:,1))
title({'AR = 20 & Ra = 8000','Local Nusselt Number'});
xlabel('Non-dimensional Height');
ylabel('Nu');
pbaspect([12 6 1])
ax=gca;
ax.FontSize=25;
set(gca,'xlabel',[])
% set(gca,'ylabel',[])
% 
% disp('Min value of U on verticle centreline : ')
% min(Uy)
% disp('Corresponding coordinates on Y-axis : ')
% y(Uy==min(Uy))
% 
% disp('Max value of U on verticle centreline : ')
% max(Uy)
% disp('Corresponding coordinates on Y-axis : ')
% y(Uy==max(Uy))
% 
% disp('Min value of V on horizontal centreline : ')
% min(Vx)
% disp('Corresponding coordinates on X-axis : ')
% x(Vx==min(Vx))
% 
% disp('Max value of V on horizontal centreline : ')
% max(Vx)
% disp('Corresponding coordinates on X-axis : ')
% x(Vx==max(Vx))
% 
% disp('Max value of Error during convergence : ')
% max(error1)
% disp('Min value of Error during convergence : ')
% min(error1)
% 
% disp('Average value of Nusselt Number : ') 
% sum(Nu(:,1))/max(size(Nu)) 