clear all
clc
close all

format long
% ===============================
% Quadrotor Parameters
% ===============================
Jxx = 6.86e-5;  
Jyy = 9.2e-5;  
Jzz = 1.366e-4;
m   = 0.068; 
kt  = 0.01;  
kq  = 7.8263e-4; 
b   = 0.062/sqrt(2);
g   = 9.81;

t1 = (Jyy - Jzz)/Jxx;
t2 = (Jzz - Jxx)/Jyy;
t3 = (Jxx - Jyy)/Jzz;

dt   = 0.01;   % sampling step (10ms)
tend = 10;     % total simulation time
th   = 1e-5;

% ===============================
% Control & physical limits
% ===============================
Tmax = 2.0*m*g;
nmax = sqrt(Tmax/(4*kt));
txymax = (Tmax/4)*2*b;
tzmax = 2*kq*nmax^2;

% ===============================
% Reference sets for variety
% ===============================
phir_list   = deg2rad([5, -5]);
thetar_list = deg2rad([0, 5]);
psir_list   = deg2rad([0, 10]);
zr_list     = [-5, -7];

run_id = 1;

for phi_ref = phir_list
    for theta_ref = thetar_list
        for psi_ref = psir_list
            for zr = zr_list
                fprintf('Simulating run %d: phir=%.1f°, thetar=%.1f°, psir=%.1f°, zr=%.1f m\n', ...
                    run_id, rad2deg(phi_ref), rad2deg(theta_ref), rad2deg(psi_ref), zr);

                % ===============================
                % Initial Conditions
                % ===============================
                x=0; y=0; z=0; 
                u=0; v=0; w=0; 
                p=0; q=0; r=0; 
                phi=0; theta=0; psi=0;
                sump=0; sumt=0; sumpsi=0; sumz=0;

                % ===============================
                % Data storage
                % ===============================
                steps = round(tend/dt);
                tim = zeros(steps,1);
                hstore = zeros(steps,1); xstore = zeros(steps,1); ystore = zeros(steps,1);
                phistore = zeros(steps,1); thetastore = zeros(steps,1); psistore = zeros(steps,1);
                pstore = zeros(steps,1); qstore = zeros(steps,1); rstore = zeros(steps,1);
                ustore = zeros(steps,1); vstore = zeros(steps,1); wstore = zeros(steps,1);
                txstore = zeros(steps,1); tystore = zeros(steps,1); tzstore = zeros(steps,1);
                Tstore = zeros(steps,1);

                % ===============================
                % Main Simulation Loop
                % ===============================
                for i = 1:steps
                    % Control gains
                    k2 = 0.1; k1 = 1.0; ki = 0.4*0.01;
                    k21 = 0.1; k11 = 1.0; ki1 = 0.4*0.01;
                    k22 = 0.1; k12 = 1.0; ki2 = 0.4*0.01;

                    % Roll loop
                    sump = sump + (phi_ref - phi);
                    pr = k1*(phi_ref - phi) + ki*sump*dt;
                    tx = k2*(pr - p);
                    tx = max(min(tx, txymax), -txymax);
                    if abs(tx) < th, tx = 0; end

                    % Pitch loop
                    sumt = sumt + (theta_ref - theta);
                    qr = k11*(theta_ref - theta) + ki1*sumt*dt;
                    ty = k21*(qr - q);
                    ty = max(min(ty, txymax), -txymax);
                    if abs(ty) < th, ty = 0; end

                    % Yaw loop
                    sumpsi = sumpsi + (psi_ref - psi);
                    rref = k12*(psi_ref - psi) + ki2*sumpsi*dt;
                    tz = k22*(rref - r);
                    tz = max(min(tz, tzmax), -tzmax);
                    if abs(tz) < th, tz = 0; end

                    % Altitude hold
                    kv = -1.0; kz1 = 2.0; kz2 = 0.15;
                    sumz = sumz + (zr - z);
                    vzr = kz1*(zr - z) + kz2*sumz*dt;
                    T = kv*(vzr - w);
                    T = min(max(T, 0.1*m*g), Tmax);

                    % Thrust-Torque conversion
                    B = [tx; ty; tz; T];
                    A = [b*kt,-b*kt,-b*kt,b*kt;
                         b*kt,b*kt,-b*kt,-b*kt;
                         kq,-kq,kq,-kq;
                         kt,kt,kt,kt];
                    C = A \ B;
                    n1 = sqrt(max(C(1),0)); n2 = sqrt(max(C(2),0));
                    n3 = sqrt(max(C(3),0)); n4 = sqrt(max(C(4),0));

                    % Dynamics integration
                    pdot = t1*q*r + tx/Jxx - 2*p;
                    qdot = t2*p*r + ty/Jyy - 2*q;
                    rdot = t3*p*q + tz/Jzz - 2*r;
                    p = p + pdot*dt; q = q + qdot*dt; r = r + rdot*dt;

                    theta = max(min(theta, deg2rad(89.9)), deg2rad(-89.9));

                    cphi = cos(phi); sphi = sin(phi);
                    ctheta = cos(theta); stheta = sin(theta);
                    cpsi = cos(psi); spsi = sin(psi);
                    ttheta = stheta / max(ctheta, 1e-6);

                    phidot = p + sphi*ttheta*q + cphi*ttheta*r;
                    thetadot = cphi*q - sphi*r;
                    psidot = (sphi*q + cphi*r) / max(ctheta, 1e-6);

                    phi = phi + phidot*dt;
                    theta = theta + thetadot*dt;
                    psi = psi + psidot*dt;

                    fz = -T;
                    udot = r*v - q*w - g*stheta - 0.1*u;
                    vdot = p*w - r*u + g*ctheta*sphi - 0.1*v;
                    wdot = q*u - p*v + g*ctheta*cphi - 0.1*w + fz/m;

                    u = u + udot*dt;
                    v = v + vdot*dt;
                    w = w + wdot*dt;

                    xdot = (cpsi*ctheta)*u + (cpsi*stheta*sphi - spsi*cphi)*v + (spsi*sphi + cpsi*stheta*cphi)*w;
                    ydot = (spsi*ctheta)*u + (spsi*stheta*sphi + cpsi*cphi)*v + (spsi*stheta*cphi - cpsi*sphi)*w;
                    zdot = -stheta*u + ctheta*sphi*v + ctheta*cphi*w;

                    x = x + xdot*dt;
                    y = y + ydot*dt;
                    z = z + zdot*dt;

                    vars = [x y z phi theta psi p q r u v w];
                    % Check for NaN or Inf values
                    if any(~isfinite(vars))
                        warning('NaN/Inf detected at step %d. Resetting invalid values to 0.', i);
                        % Replace NaN/Inf entries with 0
                        vars(~isfinite(vars)) = 0;
                        % Optionally, clamp values to avoid large blow-ups
                        vars = max(min(vars, 1e3), -1e3);
                    end
                    
                    % Unpack back into individual variables
                    x      = vars(1);
                    y      = vars(2);
                    z      = vars(3);
                    phi    = vars(4);
                    theta  = vars(5);
                    psi    = vars(6);
                    p      = vars(7);
                    q      = vars(8);
                    r      = vars(9);
                    u      = vars(10);
                    v      = vars(11);
                    w      = vars(12);

                    % Optionally wrap angles to [-pi, pi] to keep them bounded
                    phi   = wrapToPi(phi);
                    theta = wrapToPi(theta);
                    psi   = wrapToPi(psi);

                    % Store
                    tim(i) = i*dt;
                    hstore(i) = -z; 
                    xstore(i) = x; 
                    ystore(i) = y;
                    phistore(i) = rad2deg(phi); 
                    thetastore(i) = rad2deg(theta); 
                    psistore(i) = rad2deg(psi);
                    pstore(i) = p; qstore(i) = q; rstore(i) = r;
                    ustore(i) = u; vstore(i) = v; wstore(i) = w;
                    txstore(i) = tx; tystore(i) = ty; tzstore(i) = tz; 
                    Tstore(i) = T;
                end

                % ===============================
                % Save Data
                % ===============================
                T_data = table(tim, hstore, xstore, ystore, phistore, thetastore, psistore, ...
                               pstore, qstore, rstore, ustore, vstore, wstore, ...
                               txstore, tystore, tzstore, Tstore, ...
                               'VariableNames', {'Time_s','Height_m','X_m','Y_m','Roll_deg','Pitch_deg','Yaw_deg', ...
                                                 'p','q','r','u','v','w','tx','ty','tz','Thrust'});

                filename = sprintf('quad_data_run_%02d.csv', run_id);
                writetable(T_data, filename);
                fprintf('Saved %s (%d samples)\n\n', filename, length(tim));
                run_id = run_id + 1;
            end
        end
    end
end
