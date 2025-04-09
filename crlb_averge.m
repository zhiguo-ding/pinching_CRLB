clear

ct = 500000;
snr=10000;
Mvec=[1:5];
D = 10;
Dl = D*4;
height = 3; %d 
 
M=1;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength  
Ke=0.01;

delta=0;
Nwg = 2; %number of waveguides 
reso = 1000;
for m = 1 : Nwg %the y-axis locations of the waveguides (pinching antennas)
    betay(m,1) = -D/2+(m-1)*D/Nwg+D/2/Nwg;
end

Nvec = [4:4:20];
for ni = 1 : length(Nvec)
    N = Nvec(ni);
    Nused = N;
    area_length = Dl;
    for m = 1 : N/Nwg %the x-axis locations of the waveguides (pinching antennas)
        na = N/Nwg;%number of antennas on each waveguide
        betax(m,1) = -area_length/2+(m-1)*area_length/na+area_length/2/na;
    end
    area_center = 0;%Dl/4;
    betax = betax + area_center;
    
    pin_antenna = [];
    for nwg = 1 : Nwg
        na = N/Nwg;%number of antennas on each waveguide
        temp = [betax betay(nwg)*ones(N/Nwg,1)];
        pin_antenna = [pin_antenna ; temp];
    end
    %pin_antenna(1:N/2,1) = pin_antenna(1:N/2,1) + 1;
    
    
    conv_antenna = [[0:lambda:(N-1)*lambda]' D/2*ones(N,1)];
    theta_vec = [0:2*pi/N:2*pi*(N-1)/N];
    rmin = lambda/4/sin(pi/N);
    conv_antenna = [[rmin*sin(theta_vec)]' [rmin*cos(theta_vec)]'];
    
     
    for i = 1 : ct 
        loc = zeros(1,2);
        loc(:,1) = Dl*rand(1,1)-Dl/2; %length        
        loc(:,2) = D*rand(1,1)-D/2; %width, 
        loc = max(0.5*ones(1,2),loc);


        %pinching antennas
        for m = 1 : M
            xm=loc(m,1);ym=loc(m,2);
            dis_all = (xm-pin_antenna(:,1)).^2+ym^2+height^2;
            [xa,xb] = sort(dis_all,"ascend");
            sum1 = 0; sum2 = 0;
            for nind = 1 : Nused
                n = xb(nind);
                xnpin = pin_antenna(n,1);
                ynpin = pin_antenna(n,2);
                sum1 = sum1 + (xm-xnpin)^2/((xm-xnpin)^2+(ym-ynpin)^2+height^2)^2;
                sum2 = sum2 + (ym-ynpin)^2/((xm-xnpin)^2+(ym-ynpin)^2+height^2)^2;
            end 
            crlb_pinx(i) = Ke/(2*Ke+1)/sum1+Ke/(2*Ke+1)/sum2; 

            sumz1 = 0; sumz2 = 0;
            for n = 1 : N
                xncov = conv_antenna(n,1);
                yncov = conv_antenna(n,2);
                sumz1 = sumz1 + (xm-xncov)^2/((xm-xncov)^2+(ym-yncov)^2+height^2)^2;
                sumz2 = sumz2 + (ym-yncov)^2/((xm-xncov)^2+(ym-yncov)^2+height^2)^2;
            end
            crlb_convx(i) = Ke/(2*Ke+1)/sumz1+Ke/(2*Ke+1)/sumz2; 
        end
         
    end
    crlb_pin(ni) = mean(crlb_pinx);
    crlb_conv(ni) = mean(crlb_convx);

end
 
plot(Nvec,crlb_conv,Nvec,crlb_pin)