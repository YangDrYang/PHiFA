function [phi0,d_phi,b,c] = loadSAParameters_lorbeer2018experimental(filename,tau)
    [fileID, errormsg] = fopen(filename);
    if fileID<0
        disp(errormsg);
        return;
    end
    fgets(fileID);  % Ignore first line
    A = fscanf(fileID,'%f %f %f %f %f %f %f %f %f %f',[10 inf]);
    A = sortrows(A');
    if tau>A(end,1) || tau<A(1,1)
        disp('Pulse length outside of data');
        return;
    end
    phi0(1) = interp1q(A(:,1),A(:,2),tau);
    phi0(2) = interp1q(A(:,1),A(:,3),tau);
    
    d_phi(1) = interp1q(A(:,1),A(:,4),tau);
    d_phi(2) = interp1q(A(:,1),A(:,5),tau);
    
    b(1) = interp1q(A(:,1),A(:,6),tau);
    b(2) = interp1q(A(:,1),A(:,7),tau);
    b = b*10^-6; % given in N/MW
    
    c(1) = interp1q(A(:,1),A(:,8),tau);
    c(2) = interp1q(A(:,1),A(:,9),tau);
end