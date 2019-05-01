function [Q,v] = General_Q_Matrix_Newest(m_0,mu_func,sig_func,lx,ux,gridMethod,center, GridMultParam)
%GENERATES VARIANCE GRID
%
% variance process: dv_t = mu(v_t)dt + sig(v_t)dW_t
% mu_func: function handle, mu(v)
% sig_func: function handle, sig(v)
% lx: lower variance bound
% ux: upper variance bound
% m_0: number of states (grid points)
% gridMethod: 1 for uniform spacing
%             2 for Mijatovic and Pistorious
% center: applies to gridMethod 2,3,4,5: determines where the grid points should
%           cluster.. e.g. set to v_0 or estimated mean of variance process
% gridMult: alphabar = gridmult*(v_{m_0) - v_1)

nx = m_0; dx = (ux-lx)/nx;

if gridMethod == 1  % uniform grid , doesn't include lower bound
    v = lx + (1:nx)'*dx; 
elseif gridMethod == 2  %%Mijatovic and Pistorious
    v = zeros(m_0,1);
    v(1) = lx; v(m_0) = ux;
    mid = floor(m_0/2); 
    for k = 2: mid
        v(k) = center + sinh( (1-(k-1)/(mid -1))*asinh(v(1)-center) );
    end
    for k = mid+1:m_0-1
        v(k) = center + sinh( ((k - mid)/(mid))*asinh(v(m_0) - center) );
    end
elseif gridMethod == 3 %%Carr and Mayo
    x=0:1/(m_0-1):1; %map [0,1] to new grid
    alpha = GridMultParam*(ux-lx);
    c1 = asinh((lx-center)/alpha);
    c2 = asinh((ux-center)/alpha);
    v  = center + alpha*(c2*sinh(c2.*x+c1.*(1-x)));
elseif gridMethod == 4 %% Tavella and Randall
    v = zeros(m_0,1);
    v(1)   = lx;
    v(m_0) = ux;
    alpha = GridMultParam*(v(m_0) - v(1));
    c1 = asinh((v(1) - center)/alpha);
    c2 = asinh((v(m_0) - center)/alpha);
    v(2:m_0-1) = center + alpha*sinh(c2/m_0*(2:m_0-1) + c1*(1 - (2:m_0-1)/m_0));
elseif gridMethod == 5 %% Tavella and Randall PLUS we put center (e.g. v0) on grid
    tol = 1e-4; 
    v = zeros(m_0,1);
    v(1)   = lx;
    v(m_0) = ux;
    alpha = GridMultParam*(v(m_0) - v(1));
    c1 = asinh((v(1) - center)/alpha);
    c2 = asinh((v(m_0) - center)/alpha);
    vtil = zeros(m_0-1,1);
    vtil(2:m_0-2) = center + alpha*sinh(c2/(m_0 - 1)*(2:m_0-2) + c1*(1 - (2:m_0-2)/(m_0-1)));
    nnot_til = 1;
    while vtil(nnot_til) < center %Find the index left of the center point
        nnot_til = nnot_til + 1;
    end
    nnot_til = nnot_til - 1;
    v(2:nnot_til) = vtil(2:nnot_til);
    v(nnot_til + 2: m_0 - 1) = vtil(nnot_til + 1: m_0 - 2);
    if center - vtil(nnot_til)<tol %if you already have a point very close to center
        v(nnot_til) = center;
        v(nnot_til+1) = (center + vtil(nnot_til+1))/2;
    elseif vtil(nnot_til+1) - center <tol
        v(nnot_til+2) = center;
        v(nnot_til + 1) = (v(nnot_til+2) + v(nnot_til))/2;
    else
        v(nnot_til + 1) = center;
    end
end
    
%%% Now Generate Q Matrix
Q = zeros(m_0,m_0);
mu_vec = mu_func(v);
mu_plus = max(0,mu_vec);
mu_minus = max(0,-mu_vec);
sig2 = sig_func(v).^2;

if gridMethod == 0 %uniform grid
    for i=2:m_0-1
        temp = max(sig2(i) - dx*(mu_minus(i) + mu_plus(i)),0)/(2*dx^2);
        Q(i,i-1) = mu_minus(i)/dx + temp;% j = i-1
        Q(i,i+1) = mu_plus(i)/dx + temp;% j = i+1
        Q(i,i) = -Q(i,i-1) - Q(i,i+1);% j = i
    end
    Q(1,2) = abs(mu_vec(1))/dx;
    Q(m_0,m_0-1) = abs(mu_vec(m_0))/dx;       
else %nonuniform grids
    H = diff(v);
    for i=2:m_0-1
        HD = H(i-1); % down step
        HU = H(i);   % up step
        AA = max(sig2(i) - (HU*mu_plus(i) + HD*mu_minus(i)),0)/(HU+HD);
        Q(i,i-1) = (mu_minus(i) + AA)/HD; 
        Q(i,i+1) = (mu_plus(i) + AA)/HU;
        Q(i,i) = -Q(i,i-1) - Q(i,i+1);
    end
%     HU = H(1);   %Q(1,2) = abs(mu_vec(1))/HU;
%     HD = H(m_0-1);   %Q(m_0,m_0-1) = abs(mu_vec(m_0))/HD; 
    
% % %     Q(1,2) = Q(2,3);
% % %     Q(m_0,m_0-1) = Q(m_0-1,m_0-2);

%New Boundary Behavior
    HD = H(1);  HU = H(1);  %set up and down equal 
    AA = max(sig2(1) - (HU*mu_plus(1) + HD*mu_minus(1)),0)/(HU+HD);
    Q(1,2) = (mu_plus(1) + AA)/HU; 
    
    HD = H(m_0-1);  HU = H(m_0-1);  %set up and down equal 
    AA = max(sig2(m_0) - (HU*mu_plus(m_0) + HD*mu_minus(m_0)),0)/(HU+HD);
    Q(m_0,m_0-1) = (mu_plus(m_0) + AA)/HU;  
    
end
Q(1,1) = -Q(1,2);
Q(m_0,m_0) = - Q(m_0,m_0-1);



end

