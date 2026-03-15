[U,S,V] = svd(H);
F_Opt = V(:,1:Ns); % 基于SVD的预编码
W_Opt = ((1/sqrt(SNR))*(F_Opt'*H'*H*F_Opt+Ns/SNR*eye(Ns))\(F_Opt'*H'))'; % MMSE接收机