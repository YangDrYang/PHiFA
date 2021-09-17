function [C,S] = unnorm_grav_coef(C,S)
for i=0:size(C,1)-1
    for j=0:i

        temp = 1;

        if j==0
            k = 1;
        else
            k = 2;
            for m=i-j+1:i+j
                temp = temp*sqrt(m);
            end
        end

        N = sqrt(1/(k*(2.0*i+1)))*temp;

        C(i+1,j+1) = C(i+1,j+1)/N;
        S(i+1,j+1) = S(i+1,j+1)/N;

    end
end
end