function H = channel(N,Kc)
    AzoD = pi*(-.5+rand(Kc,1));
    EloD = pi*(-.5+rand(Kc,1));
    gama= sqrt(0.5)*(randn(Kc,1)+1j*randn(Kc,1));
    pi_sin_EloD=pi*sin(EloD);
    we  = pi_sin_EloD.*sin(AzoD);
    wa  = pi_sin_EloD.*cos(AzoD);
    E_streering_vector = exp(1j*we*(0:1:N-1)).';
    A_streering_vector = exp(1j*wa*(0:1:N-1)).';
    H   = E_streering_vector*diag(gama)*A_streering_vector';
end