% Joshua Ferrigno
% Department of Mechanical and Aerospace Engineering - UT Arlington

clear all; close all; clc; format shortEng;

dt = .01;
tf = 3;
tspan = 0 : dt : tf;
frame = 1;
body.n = 5;
body.units.l={'meter'};
body.units.m={'kilogram'};
body.units.t={'seconds'};

for i = 1 : body.n
    body.x(:,i) = rand(3,1);
%     body.xdot(:,i) = 2*rand(3,1)-1;
body.xdot(:,i) = [0;0;0];
    body.theta(:,i) = 2*pi*rand(3,1);
    body.omega(:,i) = 2*rand(3,1)-1;
    body.m(i) = rand;
end

xo = [body.x body.xdot] ;

[t,x] = ode45(@(t,y) Dyn(t,y,body), tspan, xo);

movi(x,frame,body.n)

function L = movi(x,frame,n)
    for i = 1 : size(x,1)
        for j = 1 : n
            plot3(x(i,3*(j-1)+1),x(i,3*(j-1)+2),x(i,3*(j-1)+3),'o','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
%             plot3(
            axis([0 frame 0 frame 0 frame])
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            hold on
            grid on
        end  
        hold off
    pause(.01)
    end
end

function F = force(x,body)
    eps = 1e-9;
    G = 1;
    n = body.n;
    % assign the positions to n-body position vectors, r
    for i = 1 : n
        r(:,i) = x(1 + 3*(i-1) : 3 + 3*(i-1));
    end
    for i = 1 : n
        m(i) = body.m(i);
    end  
    % calculate the force from body i, to bodies j
    for i = 1 : n
        for j = 1 : n
            if i ~= j 
            F(:,i,j) = -(G * m(i) * m(j) * (r(:,i) - r(:,j))) / (norm(r(:,i) - r(:,j))^3 + eps);
            end
        end
    end
    % sum the forces acting on each body by bodies j 
    F=sum(F,3);
end

function dxdt = Dyn(~,x,body)
    % n bodies * second order DE * 3d space = state variable count
    body.F = force(x,body);
    statvarc=(body.n*6);
    dxdt = zeros(statvarc,1);
    for i = 1 : statvarc/2
        dxdt(i) = x(statvarc/2 + i);
    end
    for i = 1 : statvarc/6
        dxdt((3*i-2) + statvarc/2: (3*i) + statvarc/2) = body.F(1:3,i)/body.m(i);
    end
end
    
