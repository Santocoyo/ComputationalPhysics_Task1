clc;
clear all;
format long;

function local = X0(m,bL,bU)
  n = size(bL,2);
  local = zeros(m,n);
  for i = 1:m
    for j = 1:n
      local(i,j) = bL(j) + (bU(j) - bL(j))*rand();
    endfor
  endfor
end

x0      = 0;
y0      = 0;
r1      = 1.08913;
r2      = 0.42255;
r3      = 0.96444;
r4      = 0.58781;
rcx     = 0.39137;
rcy     = 0.42950;
theta0  = 0.00000;
theta2  = X0(10,0,2*pi);
%theta2 = 4;
theta4  = X0(10,0,2*pi);

delta = 0.00001;

r1_rand = X0(100, r1-delta, r1+delta);
r2_rand = X0(100, r2-delta, r2+delta);
r3_rand = X0(100, r3-delta, r3+delta);
r4_rand = X0(100, r4-delta, r4+delta);
rcx_rand = X0(100, rcx-delta, rcx+delta);
rcy_rand = X0(100, rcy-delta, rcy+delta);

l1 = r1_rand./r2_rand;
l2 = r1_rand./r3_rand;
l3 = (r4_rand.^2-r1_rand.^2-r2_rand.^2-r3_rand.^2)/(2*r2.*r3);

for ii = 1:size(theta2)(1)
  ka(:,ii) = cos(theta2)(ii) - l1  + l2 .* cos(theta2)(ii) + l3;
  kb(:,ii) = -2*sin(theta2)(ii);
  kc(:,ii) = l1 + (l2-1) .* cos(theta2)(ii) + l3;

  theta3 = 2*atan2(-kb-sqrt(kb.^2-4*ka.*kc),2*ka);
  
  rgenx(:,ii) = x0 + r2_rand*cos(theta2(ii)+theta0) + rcx_rand.*cos(theta3(:,ii)-theta0)-rcy_rand.*sin(theta3(:,ii)-theta0);
  rgeny(:,ii) = y0 + r2_rand*sin(theta2(ii)+theta0) + rcx_rand.*sin(theta3(:,ii)-theta0)+rcy_rand.*cos(theta3(:,ii)-theta0);
  
  rgenx_avg(ii) = mean(rgenx(:,ii));
  rgeny_avg(ii) = mean(rgeny(:,ii));
  
  rgenx_std(ii) = std(rgenx(:,ii));
  rgeny_std(ii) = std(rgeny(:,ii));
endfor





%rgenx_rand = X0(100, rgenx-delta, rgenx+delta);
%rgeny_rand = X0(100, rgeny-delta, rgeny+delta);

clf;
hold on;
scatter(rgenx_avg, rgeny_avg, 15, 'r');
scatter(rgenx_avg+rgenx_std, rgeny_avg, 5, 'b');
scatter(rgenx_avg-rgenx_std, rgeny_avg, 5, 'b');
scatter(rgenx_avg, rgeny_avg+rgeny_std, 5, 'b');
scatter(rgenx_avg, rgeny_avg-rgeny_std, 5, 'b');
hold off;