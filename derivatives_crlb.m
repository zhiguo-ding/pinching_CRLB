clear
close all
 d=6;
N =16;
Nbar  =4;

delta_vec = [0.1: 0.1 :10];
for di = 1 : length(delta_vec)
    delta = delta_vec(di);

    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;sum5=0;
    for i =1 : Nbar
        for n = 1 : Nbar
            sum1 = sum1 + (n-1/2)^2/((n-1/2)^2*delta+(i-1/2)^2*delta+d^2/delta)^4*...
                ((n-1/2)^2+(i-1/2)^2-d^2/delta^2)^2;
            sum2 = sum2 + (n-1/2)^2/((n-1/2)^2*delta+(i-1/2)^2*delta+d^2/delta)^3;
            sum3 = sum3 + (n-1/2)^2/((n-1/2)^2+(i-1/2)^2+d^2/delta^2)^2;
            sum4 = sum4 + (n-1/2)^2/((n-1/2)^2*delta+(i-1/2)^2*delta+d^2/delta)^3 ...
                *((n-1/2)^2+(i-1/2)^2-d^2/delta^2);
            sum5 = sum5 +  (n-1/2)^2/((n-1/2)^2+(i-1/2)^2+delta)^2;

        end
    end
    crb2rd(di) = 6*sum1-4*d^2/delta^3*sum2;
    crb(di) = 1/sum3;
    crb1st(di) = -2*sum4;
    crbnew(di) = 1/sum5;
end
plot(delta_vec,crb1st)

 