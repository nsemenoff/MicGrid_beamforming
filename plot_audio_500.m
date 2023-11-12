clear 

fd = 48000;
fs = 2000;
load data_2000

N = 1024;

x = [-22 -12 -5 0 5 12 22; ...
     -22 -12 -5 0 5 12 22; ...
     -22 -12 -5 0 5 12 22; ...
     -22 -12 -5 0 5 12 22; ...
     -22 -12 -5 0 5 12 22; ...
     -22 -12 -5 0 5 12 22; ...
     -22 -12 -5 0 5 12 22]/100 + 0.22; % meters
 
y = [-22 -22 -22 -22 -22 -22 -22; ...
     -12 -12 -12 -12 -12 -12 -12; ...
      -5  -5  -5  -5  -5  -5  -5; ...
       0   0   0   0   0   0   0; ...
       5   5   5   5   5   5   5; ...
      12  12  12  12  12  12  12; ...
      22  22  22  22  22  22  22]/100 + 0.22; % meters

  plot(x(:),y(:),'*')
  
  %%
  writerObj = VideoWriter('test2000.avi'); %Attempt to create an avi
  open(writerObj);
  
  c=340;
  Fx = floor(fs/fd*N)+1;
  phi = -50:1:50;
  res = zeros(length(phi));
  old_res = res;
  for i=1:N:size(s,1)-N
    signal = s(i:i+N-1,:);
    for j=1:size(s,2)
        f = fft(signal(:,j));
        a = abs(f);
        a(1) = 0;
        [m,Fx] = max(a);
        disp(Fx/length(a)*fd)
        z(j) = f(Fx);
    end
    for k1=1:length(phi)
        for k2=1:length(phi)
            angle1 = phi(k1);
            angle2 = phi(k2);
            for k3=1:size(s,2)
                dx = x(k3)/c*sind(angle1);
                dy = y(k3)/c*sind(angle2);
                delay = sqrt(dx^2 + dy^2);
                delay = dx + dy;
                V(k3) = exp(-1i*2*pi*fs*delay);
            end
            tmp = z*V'/10e7;
            res(k1,k2) = abs(tmp);
        end
    end
    old_res = 0.9*old_res + 0.1*res;
    res = old_res;
    res = res/max(max(res));
    res(1,1) = 1;
    surf(res)
    view(2)
    drawnow;
    
    img = getframe;
    writeVideo(writerObj,img)
  end

  close(writerObj);
  