function Q = getGenerator_Q_MatrixOnly(v, mu_func, sig_func, gridMethod)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%% Now Generate Q Matrix
m_0 = length(v);
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
    HU = H(1);   Q(1,2) = abs(mu_vec(1))/HU;
    HD = H(m_0-1);   Q(m_0,m_0-1) = abs(mu_vec(m_0))/HD;  
end
Q(1,1) = -Q(1,2);
Q(m_0,m_0) = - Q(m_0,m_0-1);

end

