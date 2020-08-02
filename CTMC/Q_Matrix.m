function [Q,v, c_index] = Q_Matrix(m_0, mu_func,sig_func,lx,ux,gridMethod,center, GridMultParam)
%GENERATES VARIANCE GRID
%
% variance process: dv_t = mu(v_t)dt + sig(v_t)dW_t
% mu_func: function handle, mu(v)
% sig_func: function handle, sig(v)
% lx: lower grid bound
% ux: upper grid bound
% m_0: number of states (grid points)
% gridMethod: 1 for uniform spacing
%             2 for Mijatovic and Pistorious
% center: applies to gridMethod 2,3,4,5: determines where the grid points should
%           cluster.. e.g. set to v_0 or estimated mean of variance process
% gridMult: alphabar = gridmult*(v_{m_0) - v_1)

c_index = -1; % NOTE: need to return index of center point, it will return -1 if we havent coded for the case yet (e.g. uniform isnt done)
Q = zeros(m_0,m_0);

if gridMethod == 1  % uniform grid , doesn't include lower bound
    nx = m_0; dx = (ux-lx)/nx;
    v = lx + (1:nx)'*dx; 
    
elseif gridMethod == 5 || gridMethod == 8 %% Tavella and Randall PLUS we put center (e.g. v0) on grid
    % TODO: MAKE MORE EFFICIENT... we need to build grid without searching over it, and always include center point
    % Build left half, add center, build right half: grid = [left_half; center; right_half]
    tol = 1e-6; 
    v(1)   = lx;
    v(m_0) = ux;
    alpha = GridMultParam*(v(m_0) - v(1));
    c1 = asinh((v(1) - center)/alpha);
    c2 = asinh((v(m_0) - center)/alpha);
    vtil = zeros(m_0-1,1);
    vtil(2:m_0-2) = center + alpha*sinh(c2/(m_0 - 1)*(2:m_0-2) + c1*(1 - (2:m_0-2)/(m_0-1)));
    %vtil(1:m_0-1) = center + alpha*sinh(c2/(m_0 - 1)*(1:m_0-1) + c1*(1 - (1:m_0-1)/(m_0-1)));
    nnot_til = 2;
    while vtil(nnot_til) < center %Find the index left of the center point
        nnot_til = nnot_til + 1;
    end
    nnot_til = nnot_til - 1;
    v(2:nnot_til) = vtil(2:nnot_til);
    v(nnot_til + 2: m_0 - 1) = vtil(nnot_til + 1: m_0 - 2);

    if center - vtil(nnot_til)<tol %if you already have a point very close to center
        c_index = nnot_til; 
        v(c_index) = center;
        v(nnot_til+1) = (center + vtil(nnot_til+1))/2;
    elseif vtil(nnot_til+1) - center <tol
        c_index = nnot_til+2;
        v(c_index) = center;
        v(nnot_til + 2) = (v(nnot_til+2) + v(nnot_til))/2;
    else
        c_index = nnot_til + 1;
        v(c_index) = center;
    end
    
elseif gridMethod == 7
    
    v = zeros(1,m_0);
    v(m_0) = ux;
    v(1)   = lx;
    alpha = GridMultParam*(ux - lx);
    c1 = asinh((lx - center)/alpha);
    c2 = asinh((ux - center)/alpha);
    
    v(2:m_0-1) = center + alpha*sinh(c2/(m_0)*(2:m_0-1) + c1*(1 - (2:m_0-1)/(m_0)));
    %v(1:m_0) = center + alpha*sinh(c2/(m_0)*(1:m_0) + c1*(1 - (1:m_0)/(m_0)));
    c_index = 1;
    while v(c_index) < center 
        c_index = c_index + 1;
    end

    if center ~= 0
        ratio = center / v(c_index);
        v = v*ratio;  % multiply it so that center lines up with v(c_index)
    else
        v = v + center - v(c_index);
    end
    
end


%%% Now Generate Q Matrix
mu_vec = mu_func(v);
mu_plus = max(0,mu_vec);
mu_minus = max(0,-mu_vec);
sig2 = sig_func(v).^2;


if gridMethod == 1 %uniform grid
    for i=2:m_0-1
        temp = max(sig2(i) - dx*(mu_minus(i) + mu_plus(i)),0)/(2*dx^2);
        Q(i,i-1) = mu_minus(i)/dx + temp;% j = i-1
        Q(i,i+1) = mu_plus(i)/dx + temp;% j = i+1
        Q(i,i) = -Q(i,i-1) - Q(i,i+1);% j = i
    end
    Q(1,2) = abs(mu_vec(1))/dx;
    Q(m_0,m_0-1) = abs(mu_vec(m_0))/dx; 
elseif gridMethod == 8
    
    H = diff(v);
    HD = H(1);  HU = H(1);  %set up and down equal 
    AA = max(sig2(1) - (HU*mu_plus(1) + HD*mu_minus(1)),0)/(HU+HD);
    Q(1,2) = (mu_plus(1) + AA)/HU; 
    Q(1,1) = -Q(1,2);

    HD = H(m_0-1);  HU = H(m_0-1);  %set up and down equal 
    AA = max(sig2(m_0) - (HU*mu_plus(m_0) + HD*mu_minus(m_0)),0)/(HU+HD);
    Q(m_0,m_0-1) = (mu_plus(m_0) + AA)/HU;  
    Q(m_0,m_0) = - Q(m_0,m_0-1);
    
    for i=2:m_0-1
        dvU = v(i-1) - v(i); dvD = v(i+1) - v(i);
        C = [1 1 1; dvU 0 dvD; dvU^2 0 dvD^2];
        z = [0; mu_vec(i); sig2(i)];
        hrow = C\z;
        Q(i,i-1) = hrow(1);  Q(i,i) = hrow(2); Q(i,i + 1) = hrow(3);
    end
    
    
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

