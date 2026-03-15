[IndexAt,IndexAr] = findSteeringVector(H,At,Ar,Ns);
F_BS = []; W_BS = [];
for n = 1:Ns 
    F_BS = [F_BS At(:,IndexAt(n))];
    W_BS = [W_BS Ar(:,IndexAr(n))];
end