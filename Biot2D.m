% function Biot2D
clc
close all
clear all
clear all -globals
Globals2D

N = 5;
K1D = 8;
c_flag = 0;
FinalTime = 1;
CFL = .75;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

BuildPeriodicMaps2D(2,2);

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

global Vq Pq
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;


%% params setup

U = cell(8,1);
for i = 1:8
    U{i} = zeros(Np,K);
end
U{1} = exp(-10^2*(x.^2+y.^2));
% U{2} = exp(-10^2*(x.^2+y.^2));
% U{5} = exp(-10^2*(x.^2+y.^2));

% randomly chosen Es, Ev
Es = [ones(3)+2*eye(3) zeros(3,1)
    zeros(1,3) 2];
Ev = [2 0 -1 0 ;0 2 0 -1; -1 0 2 0; 0 -1 0 2];

A0 =blkdiag(Es,Ev);
global invA0
invA = inv(A0);
for i = 1:size(invA,1)
    for j = 1:size(invA,2)
        invA0{i,j} = invA(i,j).*(1+.25*sin(pi*xq).*sin(pi*yq));
    end
end


%%

time = 0;

% Runge-Kutta residual storage
res = cell(8,1);
for fld = 1:8
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = CFL*2/CNh;

% outer time step loop
tstep = 0;
figure
% colormap(gray)
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        rhs = RHS2D(U,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        for fld = 1:length(U)
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs{fld};
            U{fld} = U{fld} + rk4b(INTRK) * res{fld};
        end                
        
    end;
    
    if 1 
        clf
        pp = U{1};
        vv = Vp*pp;
        %         vv = abs(vv);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
%         caxis([-.1 .2])
%         axis([0 1 0 1 -10 10])
        %         PlotField2D(N+1, x, y, pp); view(2)
        title(sprintf('time = %f',time))
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

% axis off
% view(0,0)


function rhs = RHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

% Define field differences at faces
Nfld = length(U);
dU = cell(Nfld,1);
Ux = cell(Nfld,1);
Uy = cell(Nfld,1);
for fld = 1:Nfld
    dU{fld} = zeros(Nfp*Nfaces,K); dU{fld}(:) = U{fld}(vmapP)-U{fld}(vmapM);
    Ur = Dr*U{fld}; Us = Ds*U{fld};    
    Ux{fld} = rx.*Ur + sx.*Us + .5*LIFT*(Fscale.*dU{fld}.*nx);
    Uy{fld} = ry.*Ur + sy.*Us + .5*LIFT*(Fscale.*dU{fld}.*ny);    
end

idx = [7 5 1 6 2 4 1 1];
signx = [1 -1 0 -1 -1 -1 1 0];
idy = [8 1 6 5 4 3 1 1];
signy = [1 0 -1 -1 -1 -1 0 1];

% penalty fluxes: An^2 * [[Q]] = An* (An*[[Q]])
f{1} = nx.*dU{7} + ny.*dU{8};
f{2} = -nx.*dU{5};
f{3} = -ny.*dU{6};
f{4} = -(ny.*dU{5} + nx.*dU{6});
f{5} = -(nx.*dU{2} + ny.*dU{4});
f{6} = -(ny.*dU{3} + nx.*dU{4});
f{7} = nx.*dU{1};
f{8} = ny.*dU{1};

% apply An again to An.*[[Q]]
f{1} = nx.*f{7} + ny.*f{8};
f{2} = -nx.*f{5};
f{3} = -ny.*f{6};
f{4} = -(ny.*f{5} + nx.*f{6});
f{5} = -(nx.*f{2} + ny.*f{4});
f{6} = -(ny.*f{3} + nx.*f{4});
f{7} = nx.*f{1};
f{8} = ny.*f{1};

tau = .25; % choose more carefully!

global Vq Pq
rr = cell(Nfld,1);
for fld = 1:Nfld    
    rhsx = signx(fld)*Ux{idx(fld)};
    rhsy = signy(fld)*Uy{idy(fld)};
    rr{fld} = rhsx + rhsy - tau*.5*LIFT*(Fscale.*f{fld});
    rr{fld} = Vq*rr{fld};
end

global invA0
rhs = cell(Nfld,1);
for ii = 1:Nfld
    rhs{ii} = zeros(Np,K);   
    for jj = 1:Nfld
        rhs{ii} = rhs{ii} - Pq*(invA0{ii,jj}.*rr{jj}); % negative sign
    end
end

end

