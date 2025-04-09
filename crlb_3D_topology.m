clear
close all

ct = 500000;
snr=10000;
Mvec=[1:5];
D = 10;
Dl = D*4;
height = 3; %d 
N=20;
M=1;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength  
Ke=2;
Nused = N;
delta=0;
Nwg = 2; %number of waveguides 
reso = 1000;
for m = 1 : Nwg %the y-axis locations of the waveguides (pinching antennas)
    betay(m,1) = -D/2+(m-1)*D/Nwg+D/2/Nwg;
end

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

stepx = Dl/reso;
xvec = [-Dl/2: stepx: Dl/2];
stepy = D/reso;
yvec = [-D/2: stepy: D/2];
for i = 1 : length(xvec) 
    for j = 1 : length(yvec)
        xm = xvec(i);
        ym = yvec(j);        
        loc = [xm ym]; 

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
            crlb_pin(i,j) = Ke/(2*Ke+1)/sum1+Ke/(2*Ke+1)/sum2;
            if crlb_pin(i,j) > 101
                crlb_pin(i,j) = 101;
            end

            sumz1 = 0; sumz2 = 0;
            for n = 1 : N
                xncov = conv_antenna(n,1);
                yncov = conv_antenna(n,2);
                sumz1 = sumz1 + (xm-xncov)^2/((xm-xncov)^2+(ym-yncov)^2+height^2)^2;
                sumz2 = sumz2 + (ym-yncov)^2/((xm-xncov)^2+(ym-yncov)^2+height^2)^2;
            end
            crlb_conv(i,j) = Ke/(2*Ke+1)/sumz1+Ke/(2*Ke+1)/sumz2;
            if crlb_conv(i,j) > 2001
                crlb_conv(i,j) = 2001;
            end
        end
    end
end
 
%end
mesh(yvec,xvec, crlb_pin)
figure
mesh(yvec,xvec, crlb_conv)