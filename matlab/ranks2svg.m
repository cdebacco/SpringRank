% Code by Daniel Larremore
% Santa Fe Institute
% larremore@santafe.edu
% http://danlarremore.com
% v3
function [energy] = ranks2svg(A,s,filename)
[r,c,v] = find(A);
energy = (s(r)-s(c)-1).^2;
energy = lin(energy/max(energy),0.05,0.3);

% wid = 800; % Must be at least 400
hei = 800; % Must be at least 400
aspectRatio = 0.4;
wid = hei*aspectRatio;
s = max(s)-s; %flip orientation so top ranked shows up highest
y = (s-min(s))/(max(s)-min(s))*(7/8)*hei+hei/16;

% colors
c_fill    = [10,10,10]; %color - circle fill
% c_fill = 'none';
c_stroke  = [10,10,10]; %color - circle stroke
c_down    = [41,143,170]; %color - edges down
c_up      = [170,48,41]; %color - edges up
s_min = 1; %edge stroke min
s_max = 4; %edge stroke max
r_min = 1; %circle min radius
r_max = 3; %circle max radius
swirl = 0.75;

% bind sizes of circles to k.
Q = A;  
Q(1:size(Q,1)+1:end)=0; % kill diagonal
k = sum(Q,2); % bind to out degree. 
k = lin(k,r_min,r_max);

% transform edge weights in v to range [s_min,s_max]
w = lin(v,s_min,s_max);

fid = fopen(filename,'w');
fprintf(fid,'<svg width="%i" height="%i" id="vis">',wid,hei);
fprintf(fid,'\n');

% Paths
for i=1:length(v)
    fr = y(r(i));
    to = y(c(i));
    anchor_y = swirl*fr+(1-swirl)*to;
    anchor_x = (wid/2)-sign(fr-to)*abs(fr-to)*aspectRatio;
    if sign(fr-to) < 0
        rgb = getColor_identity(c_down,energy);
        opacity = 0.15;
        fprintf(fid,'<path d="M%f %f Q %f %f %f %f" stroke-width="%f" stroke="rgb(%i,%i,%i)" stroke-opacity="%f" fill="none"/>',...
            wid/2,fr,...
            anchor_x,anchor_y,wid/2,to,...
            w(i),...
            rgb(1),...
            rgb(2),...
            rgb(3),...
            opacity);
    else
        rgb = getColor_identity(c_up,energy);
        opacity = 0.15;
        fprintf(fid,'<path d="M%f %f Q %f %f %f %f" stroke-width="%f" stroke="rgb(%i,%i,%i)" stroke-opacity="%f" fill="none"/>',...
            wid/2,fr,...
            anchor_x,anchor_y,wid/2,to,...
            w(i),...
            rgb(1),...
            rgb(2),...
            rgb(3),...
            opacity);
    end
    fprintf(fid,'\n');
end

% Circles
n = size(A,1);
for i=1:n
    if strcmp(c_fill,'none')==1
        fprintf(fid,'<circle class="node" r="%f" style="fill: none; stroke: rgb(%i,%i,%i); stroke-width:1" cx="%f" cy="%f"/>',...
        k(i),c_stroke(1),c_stroke(2),c_stroke(3),wid/2,y(i) );
    else
    fprintf(fid,'<circle class="node" r="%f" style="fill: rgb(%i,%i,%i); stroke: rgb(%i,%i,%i); stroke-width:1" cx="%f" cy="%f"/>',...
        k(i),c_fill(1),c_fill(2),c_fill(3),c_stroke(1),c_stroke(2),c_stroke(3),wid/2,y(i) );
    end
    fprintf(fid,'\n');
end

fprintf(fid,'</svg>');

fclose(fid);

end

function [y] = lin(x,m,M)
y = (x - min(x))/(max(x)-min(x))*(M-m)+m;
end
function [rgb] = getColor_linearToWhite(rgb_base,scalar)
M = max(rgb_base);
rgb = round(rgb_base+(M-rgb_base)*scalar);
end
function [rgb] = getColor_identity(rgb_base,scalar)
rgb = rgb_base;
end
