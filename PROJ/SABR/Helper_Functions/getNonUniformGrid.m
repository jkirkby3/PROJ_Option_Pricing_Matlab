function v = getNonUniformGrid(m_0, lx, ux, gridMethod, center, manualPoint, gridMultParam)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
    alpha = 0.8*(ux-lx);
    c1 = asinh((lx-center)/alpha);
    c2 = asinh((ux-center)/alpha);
    v  = center + alpha*(c2*sinh(c2.*x+c1.*(1-x)));
elseif gridMethod == 4 %% Tavella and Randall
    v = zeros(m_0,1);
    v(1)   = lx;
    v(m_0) = ux;
    alpha = gridMultParam*(v(m_0) - v(1));
    c1 = asinh((v(1) - center)/alpha);
    c2 = asinh((v(m_0) - center)/alpha);
    v(2:m_0-1) = center + alpha*sinh(c2/m_0*(2:m_0-1) + c1*(1 - (2:m_0-1)/m_0));
elseif gridMethod == 5 %% Tavella and Randall PLUS we put manualPoint (e.g. v0) on grid
    tol = 1e-7; 
    v = zeros(m_0,1);
    v(1)   = lx;
    v(m_0) = ux;
    alpha = gridMultParam*(v(m_0) - v(1));
    c1 = asinh((v(1) - center)/alpha);
    c2 = asinh((v(m_0) - center)/alpha);
    vtil = zeros(m_0-1,1);
    vtil(2:m_0-2) = center + alpha*sinh(c2/(m_0 - 1)*(2:m_0-2) + c1*(1 - (2:m_0-2)/(m_0-1)));
    nnot_til = 1;
    while vtil(nnot_til) < manualPoint %Find the index left of the manualPoint point
        nnot_til = nnot_til + 1;
    end
    nnot_til = nnot_til - 1;
    v(2:nnot_til) = vtil(2:nnot_til);
    v(nnot_til + 2: m_0 - 1) = vtil(nnot_til + 1: m_0 - 2);
    if manualPoint - vtil(nnot_til)<tol %if you already have a point very close to manualPoint
        v(nnot_til) = manualPoint;
        v(nnot_til+1) = (manualPoint + vtil(nnot_til+1))/2;
    elseif vtil(nnot_til+1) - manualPoint <tol
        v(nnot_til+2) = manualPoint;
        v(nnot_til + 1) = (v(nnot_til+2) + v(nnot_til))/2;
    else
        v(nnot_til + 1) = manualPoint;
    end
end

end

