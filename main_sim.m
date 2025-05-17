clear; clc; close all;

% PARAMETERS
Lbase      = 3; H = 0;
base1      = [-Lbase/2; H];
base2      = [ Lbase/2; H];
ell        = [1;1;1];
m          = [4;2;1];

ke    = 30; be        = 14;    % PD gains
gamma = 0.2; Gamma    = 0.8;   % barrier thresholds
ka    = 100;           % avoidance gain
eps_p = 5e-2;          % perturbation magnitude
vel_thresh = 5e-2;     % deadlock velocity threshold

theta1_0   = [11*pi/18; -pi/3;  -pi/3];
theta2_0   = [11*pi/18;  2*pi/3; -pi/6];
theta1_des = [10*pi/12;  -pi/2;  -pi/6];
theta2_des = [13*pi/18;  -pi/3;  -pi/8];

dt    = 0.005;
T     = 3;
time  = 0:dt:T;
n     = 3;
NA    = 2;
Ntot  = n*NA;


V_hist   = zeros(2, numel(time));          
err_hist = zeros(Ntot, numel(time), 2);    

% RUN SIMULATIONS
theta_hist = zeros(Ntot, numel(time), 2);
tol_ang = 1e-2;
tol_vel = 1e-2;

for simIdx = 1:2
  useAvoid = (simIdx == 2);
  theta    = zeros(Ntot, numel(time));
  thdot    = zeros(Ntot, numel(time));
  theta(:,1) = [theta1_0; theta2_0];
  endStep = numel(time);
  
  for k = 1:numel(time)-1
    th  = theta(:,k);
    thd = thdot(:,k);
    th1 = th(1:n);     d1 = thd(1:n);
    th2 = th(n+1:end); d2 = thd(n+1:end);

    % dynamics 
    [M1,C1] = planarDynamics(th1,d1);
    [M2,C2] = planarDynamics(th2,d2);

    %nominal PD torque
    u_nom1 = -ke*(th1 - theta1_des) - be*d1;
    u_nom2 = -ke*(th2 - theta2_des) - be*d2;
    u_nom  = [u_nom1; u_nom2];

    % avoidance torque
    u_avoid = zeros(Ntot,1);
    if useAvoid
      ends  = computePlanarEndpoints(th, ell, base1, base2);
      pairs = nchoosek(1:2*n,2);
      for j = 1:size(pairs,1)
        i1 = pairs(j,1); i2 = pairs(j,2);
        P1 = ends(i1,1:2)'; Q1 = ends(i1,3:4)';
        P2 = ends(i2,1:2)'; Q2 = ends(i2,3:4)';
        Phi = segmentDistance(P1,Q1,P2,Q2);
        if Phi < Gamma
          % central difference grad Phi
          gradPhi = zeros(Ntot,1);
          delta   = 1e-6;
          for i = 1:Ntot
            thp = th; thm = th;
            thp(i) = thp(i) + delta;
            thm(i) = thm(i) - delta;
            ep = computePlanarEndpoints(thp, ell, base1, base2);
            em = computePlanarEndpoints(thm, ell, base1, base2);
            Phi_p = segmentDistance(ep(i1,1:2)', ep(i1,3:4)', ep(i2,1:2)', ep(i2,3:4)');
            Phi_m = segmentDistance(em(i1,1:2)', em(i1,3:4)', em(i2,1:2)', em(i2,3:4)');
            gradPhi(i) = (Phi_p - Phi_m)/(2*delta);
          end
          factor   = 2*(Phi - Gamma)/((Gamma - gamma)*(Phi - gamma)^3);
          u_avoid = u_avoid + ka * factor * gradPhi;
        end
      end
    end

    % perturbation torque
    u_pert = zeros(Ntot,1);
    target = [theta1_des; theta2_des];
    for i = 1:Ntot
      if abs(thd(i)) < vel_thresh
        u_pert(i) = eps_p * sign(target(i) - th(i));
      end
    end

    % combine into total torque
    u_total = u_nom - u_avoid + u_pert;

    % integrate
    u1    = u_total(1:n);
    u2    = u_total(n+1:end);
    thdd1 = M1 \ (u1 - C1*d1);
    thdd2 = M2 \ (u2 - C2*d2);

    thdot(:,k+1) = thd + dt*[thdd1; thdd2];
    theta(:,k+1) = th  + dt*thdot(:,k+1);

    % Lyapunov candidate
    V1 = 0.5*d1' * M1 * d1 + 0.5*ke * norm(th1 - theta1_des)^2;
    V2 = 0.5*d2' * M2 * d2 + 0.5*ke * norm(th2 - theta2_des)^2;
    V_hist(simIdx, k+1) = V1 + V2;
    err_hist(:,k+1,simIdx) = theta(:,k+1) - target;

    err1 = max(abs(theta(1:3,k+1)   - theta1_des));
    err2 = max(abs(theta(4:6,k+1)   - theta2_des));
    vel1 = max(abs(thdot(1:3,k+1)));
    vel2 = max(abs(thdot(4:6,k+1)));
    if err1 < tol_ang && err2 < tol_ang && vel1 < tol_vel && vel2 < tol_vel
      fprintf('Sim %d converged at t = %.3f s\n', simIdx, time(k+1));
      endStep = k+1;
      break;
    end
  end

  
  if endStep < numel(time)
    theta_hist(:,:,simIdx) = [theta(:,1:endStep), repmat(theta(:,endStep),1,numel(time)-endStep)];
    V_hist(simIdx,endStep+1:end)    = V_hist(simIdx,endStep);
    err_hist(:,endStep+1:end,simIdx) = repmat(err_hist(:,endStep,simIdx),1,numel(time)-endStep);
  else
    theta_hist(:,:,simIdx) = theta;
  end
end

% ANIMATE & SNAPSHOTS
snapDir = fullfile(pwd,'se525');
if ~exist(snapDir,'dir')
  mkdir(snapDir);
end
videoFilename = fullfile(snapDir, 'manipulators_animation.mp4');
v = VideoWriter(videoFilename, 'MPEG-4');
v.FrameRate = 20;    
open(v);

nextSnapshot = 0;
figure('Color','w','Position',[100 200 900 400]);
for k = 1:5:length(time)
  clf;
  for simIdx = 1:2
    subplot(1,2,simIdx); hold on;
    if simIdx==1, title('Without Avoidance');
    else           title('With Avoidance'); end
    thk = theta_hist(:,k,simIdx);
    plotArm(thk(1:3), ell, base1, 'r');
    plotArm(thk(4:6), ell, base2, 'b');
    xlabel('X'); ylabel('Y');
    axis equal; xlim([-3 3]); ylim([-1 4]);
    grid on;
    text(-2.8,3.8,sprintf('t = %.2f', time(k)));
  end
  drawnow;
  frame = getframe(gcf);
  writeVideo(v, frame);
  if time(k)>=nextSnapshot
    fname = fullfile(snapDir,sprintf('snapshot_%03.0fms.svg',nextSnapshot*1000));
    print(gcf,'-dsvg',fname,'-r300');
    nextSnapshot = nextSnapshot + 0.05;
  end
end
close(v)
% LYAPUNOV CANDIDATE
figure('Name','Lyapunov Candidate','Color','w');
plot(time,V_hist(1,:),'-r','LineWidth',1.5); hold on;
plot(time,V_hist(2,:),'-b','LineWidth',1.5);
xlabel('Time (s)'); ylabel('V(t)');
legend('Without Avoidance','With Avoidance','Location','northeast');
title('Lyapunov Candidate Over Time');

% JOINT ANGLE CONVERGENCE
figure('Name','Joint Angle Convergence','Color','w');
theta_des_all = repmat([theta1_des; theta2_des],1,numel(time));
for j=1:Ntot
  subplot(Ntot,1,j);
  plot(time,theta_hist(j,:,2),'-b','LineWidth',1.2); hold on;
  plot(time,theta_des_all(j,:), '--k','LineWidth',1);
  ylabel(sprintf('\\theta_{%d}',j));
  if j==1, legend('With Avoidance','Desired','Location','best'); end
  if j==Ntot, xlabel('Time (s)'); end
end
sgtitle('Joint Angles Converging to Desired (With Avoidance)');

function [M,C] = planarDynamics(th, thd)
  t2 = th(2); t3 = th(3);
  d1 = thd(1); d2 = thd(2); d3 = thd(3);
  M11 = cos(t2+t3)+4*cos(t2)+cos(t3)+23/4;
  M12 = 0.5*cos(t2+t3)+2*cos(t2)+cos(t3)+7/4;
  M13 = 0.5*cos(t2+t3)+0.5*cos(t3)+1/4;
  M22 =            cos(t3)           +7/4;
  M23 =            0.5*cos(t3)      +1/4;
  M33 =                         1/4;
  M   = [M11 M12 M13; M12 M22 M23; M13 M23 M33];

  C11 = -d2*(0.5*sin(t2+t3)+2*sin(t2)) -0.5*d3*(sin(t2+t3)+sin(t3));
  C12 = -d1*(0.5*sin(t2+t3)+2*sin(t2)) -d2*(0.5*sin(t2+t3)+2*sin(t2)) -0.5*d3*(sin(t2+t3)+sin(t3));
  C13 = -0.5*(d1+d2+d3)*sin(t3);
  C21 =  0.5*d1*(sin(t2+t3)+sin(t3)) -0.5*d3*sin(t3);
  C22 = -0.5*d3*sin(t3);
  C23 =  0;
  C31 =  0.5*(d1+d2)*sin(t3);
  C32 =  0; C33=0;
  C   = [C11 C12 C13; C21 C22 C23; C31 C32 C33];
end

function ends = computePlanarEndpoints(th, ell, b1, b2)
  n = numel(ell);
  ends = zeros(2*n,4);
  for a=1:2
    base = (a==1)*b1 + (a==2)*b2;
    the  = th((a-1)*n+(1:n));
    pos  = base;
    for j=1:n
      nxt = pos + ell(j)*[cos(sum(the(1:j))); sin(sum(the(1:j)))];
      idx = (a-1)*n + j;
      ends(idx,:) = [pos' nxt'];
      pos = nxt;
    end
  end
end

function d = segmentDistance(P1,Q1,P2,Q2)
  u = Q1-P1; v = Q2-P2; w = P1-P2;
  a = dot(u,u); b = dot(u,v); c=dot(v,v);
  d0=dot(u,w); e=dot(v,w); D=a*c-b^2;
  if abs(D)>1e-8
    s=(b*e-c*d0)/D; t=(a*e-b*d0)/D;
  else
    s=0; t=0;
  end
  s = max(0,min(1,s)); t = max(0,min(1,t));
  C1 = P1 + s*u; C2 = P2 + t*v;
  d = norm(C1 - C2);
end

function plotArm(th, ell, base, col)
  pts = base';
  for i=1:numel(ell)
    ang = sum(th(1:i));
    pts(end+1,:) = pts(end,:) + ell(i)*[cos(ang), sin(ang)];
  end
  plot(pts(:,1),pts(:,2),'-o','Color',col,'LineWidth',2,'MarkerFaceColor',col);
end
