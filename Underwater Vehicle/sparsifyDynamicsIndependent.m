function Xi = sparsifyDynamicsIndependent(Theta,dXdt,lambda,n)

% Initial guess using least-squares
Xi = pinv(Theta)*dXdt;  

% state-dependent lambda is sparsifying knob
for k=1:10
    for ind = 1:n    % n is state dimension
        if ind == 7 
            keyboard
        end
        smallinds = (abs(Xi(:,ind))<lambda(ind));   % find small coefficients
        Xi(smallinds,ind)=0;                        % and threshold 
        biginds = ~smallinds;
        % Regress dynamics onto remaining terms to find sparse Xi 
        Xi(biginds,ind) = pinv(Theta(:,biginds))*dXdt(:,ind); 
    end
end