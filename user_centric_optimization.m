clear


ct = 500000;
snr=10000;
Mvec=[1:5]; 
height = 6; %d 
N=4;
 f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength  
Ke=2;
delta=0;
Nwg = 2; %number of waveguides 
Ntilde = N/Nwg;
reso = 1000;

delta_vec = [0.1: 0.1 : 10];
for di = 1 : length(delta_vec)
    delta = delta_vec(di);
    D = (Ntilde)*delta;
    for m = 1 : Nwg %the y-axis locations of the waveguides (pinching antennas)
        betay(m,1) = -D/2+(m-1)*D/Nwg+D/2/Nwg;
    end

    area_length = D;
    for m = 1 : N/Nwg %the x-axis locations of the waveguides (pinching antennas)
        na = N/Nwg;%number of antennas on each waveguide
        betax(m,1) = -area_length/2+(m-1)*area_length/na+area_length/2/na;
    end
    area_center = 0;
    betax = betax + area_center;
    
    pin_antenna = [];
    for nwg = 1 : Nwg
        na = N/Nwg;%number of antennas on each waveguide
        temp = [betax betay(nwg)*ones(N/Nwg,1)];
        pin_antenna = [pin_antenna ; temp];
    end
 
 
        xm = 0;
        ym = 0;        
        loc = [xm ym]; 

        %pinching antennas
        m=1;
            xm=loc(m,1);ym=loc(m,2);
            dis_all = (xm-pin_antenna(:,1)).^2+ym^2+height^2;
            [xa,xb] = sort(dis_all,"ascend");
            sum1 = 0; sum2 = 0;
            for nind = 1 : N
                n = xb(nind);
                xnpin = pin_antenna(n,1);
                ynpin = pin_antenna(n,2);
                sum1 = sum1 + (xm-xnpin)^2/((xm-xnpin)^2+(ym-ynpin)^2+height^2)^2;
                sum2 = sum2 + (ym-ynpin)^2/((xm-xnpin)^2+(ym-ynpin)^2+height^2)^2;
            end 
            crlb_pin(di) = Ke/(2*Ke+1)/sum1+Ke/(2*Ke+1)/sum2;
        
end

%analysis
delta = sqrt(2)*height;
    D = (Ntilde)*delta;
    for m = 1 : Nwg %the y-axis locations of the waveguides (pinching antennas)
        betay(m,1) = -D/2+(m-1)*D/Nwg+D/2/Nwg;
    end

    area_length = D;
    for m = 1 : N/Nwg %the x-axis locations of the waveguides (pinching antennas)
        na = N/Nwg;%number of antennas on each waveguide
        betax(m,1) = -area_length/2+(m-1)*area_length/na+area_length/2/na;
    end
    area_center = 0;
    betax = betax + area_center;
    
    pin_antenna = [];
    for nwg = 1 : Nwg
        na = N/Nwg;%number of antennas on each waveguide
        temp = [betax betay(nwg)*ones(N/Nwg,1)];
        pin_antenna = [pin_antenna ; temp];
    end
 
 
        xm = 0;
        ym = 0;        
        loc = [xm ym]; 

        %pinching antennas
        m=1;
            xm=loc(m,1);ym=loc(m,2);
            dis_all = (xm-pin_antenna(:,1)).^2+ym^2+height^2;
            [xa,xb] = sort(dis_all,"ascend");
            sum1 = 0; sum2 = 0;
            for nind = 1 : N
                n = xb(nind);
                xnpin = pin_antenna(n,1);
                ynpin = pin_antenna(n,2);
                sum1 = sum1 + (xm-xnpin)^2/((xm-xnpin)^2+(ym-ynpin)^2+height^2)^2;
                sum2 = sum2 + (ym-ynpin)^2/((xm-xnpin)^2+(ym-ynpin)^2+height^2)^2;
            end 
            crlb_ana = Ke/(2*Ke+1)/sum1+Ke/(2*Ke+1)/sum2;

%end
plot(delta_vec,crlb_pin, delta_vec, crlb_ana*ones(length(delta_vec),1) ) 