function [u_mpc] = main_mpc_undrivable_obstacle(t, input, params, controller_on, points_x_r, points_x_l, points_y_r, points_y_l, obstacle_detected, ds_x, ds_y)

%           1   2  3  4    5       6       7      8
% input = beta, r, v, yaw, df_ref, u_prev, x_cor, y_cor
% params = [lf, lr, m, I, ha]
persistent opti r x0 u_prev v x u as k0L k1L k2L k0R k1R k2R sfr sfl srr srl radfl1 radfr1 radfl2 radfr2 dsx_1 dsx_2 dsy_1 dsy_2 radius_1 radius_2 sD obs_on_1 obs_on_2 caf car bigslack

% vehicle parameters
lf   = params(1); % GC-front distance
lr   = params(2); % GC-rear distance
m    = params(3); % mass
I    = params(4); % inertia
ha   = params(5); % half of the axle
wheel_dist_left = ha;   % distance to a wheel on y-axis
wheel_dist_right = -ha;
rad_ius = 1.6;

if t == 0
    opti = casadi.Opti(); % definition of the controller object
    roadPadding = 0.2;
    % mpc stuff
    nx = 5; % number of states
    nu = 1; % number of inputs
    
    % pred. horizon
    N = 15;
    % basic oprimization vars 
    r = opti.parameter(nu, 1); % references for inputs
    x0 = opti.parameter(nx, 1); % initial state
    x = opti.variable(nx, N+1); % state variables
    u = opti.variable(nu, N); % input variables
    as = opti.variable(nu, N+1); % parameter for abs substitution
    u_prev = opti.parameter(nu, 1); % previous input
    % bordes curve param
    k2L = opti.parameter(); 
    k1L = opti.parameter();
    k0L = opti.parameter();
    k2R = opti.parameter();
    k1R = opti.parameter();
    k0R = opti.parameter();
    % radius of obstacle for each wheel
    radius_1 = opti.parameter();
    radius_2 = opti.parameter();
    obs_on_1 = opti.parameter();
    obs_on_2 = opti.parameter();
    %distance to obstacle for each wheel
    dsx_1 = opti.parameter();
    dsx_2 = opti.parameter();
    dsy_1 = opti.parameter();
    dsy_2 = opti.parameter();
    % slack vars
    sfl = opti.variable(N+1, 1);
    sfr = opti.variable(N+1, 1);
    srl = opti.variable(N+1, 1);
    srr = opti.variable(N+1, 1);
    radfl1 = opti.variable(N+1, 1);
    radfr1 = opti.variable(N+1, 1);
    radfl2 = opti.variable(N+1, 1);
    radfr2 = opti.variable(N+1, 1);
    sD = opti.variable(N, 1);

    % weights
    steering_wheel_track = 100;
    R = steering_wheel_track;
    Qs   = 1e10; % slack weight 
    Qo   = 1e3;
    Qenv_ar = 1e4;
    change_pen = 1e2;
    v = opti.parameter(); % velocity at the moment
    caf = opti.parameter();
    car = opti.parameter();
    Ts = 0.04;
    deltaFmax = 0.65;
    deltaFmin = -0.65;
    % slew restrictions
    delta_slew = deg2rad(240)*Ts; % rad/s
    obj = 0;
    
    A = [-(caf+car)/m/v, (car*lr-caf*lf)/m/v/v-1;
         (car*lr-caf*lf)/I, -(lf^2*caf+car*lr^2)/I/v];
    B = [caf/m/v;
         caf*lf/I];

    Ad = eye(2) + Ts*A;
    Bd = Ts*B;
    
    disf = sqrt(lf^2+wheel_dist_left^2);
    cos_angl_left_front = cos(atan2(wheel_dist_left,lf));
    sin_angl_left_front = sin(atan2(wheel_dist_left,lf));
    cos_angl_right_front = cos(atan2(wheel_dist_right,lf));
    sin_angl_right_front = sin(atan2(wheel_dist_right,lf));
    
    points_count = 1;
    bigslack = opti.variable(N+1, 2*(points_count+1));
    
    for k = 1:N
        % basic dynamics and inputs
        if (k == 1)
            % dynamic model beta,r,x,y,psi
            opti.subject_to(x(1:2, k) == Ad*x0(1:2) + Bd*u(k));
            opti.subject_to(x(3:5, k) == Ts*[v*1; v*x0(1); x0(2)]);
        else
            opti.subject_to(x(1:2, k) == Ad*x(1:2, k-1) + Bd*u(k));
            opti.subject_to(x(3:5, k) == x(3:5, k-1) + Ts*[v*1; v*(x(5, k-1) + x(1, k-1)); x(2)]);
        end
        if (k == 1)
            xfl = x(3,k)+disf*cos_angl_left_front;
            yfl = x(4,k)+disf*sin_angl_left_front;
            xfr = x(3,k)+disf*cos_angl_right_front;
            yfr = x(4,k)+disf*sin_angl_right_front;
        else
            xfl = x(3,k)+disf*(cos_angl_left_front-sin_angl_left_front*x(5,k-1));% cos(p+psi)
            yfl = x(4,k)+disf*(sin_angl_left_front+cos_angl_left_front*x(5,k-1));% sin(p+psi)
            xfr = x(3,k)+disf*(cos_angl_right_front-sin_angl_right_front*x(5,k-1));% cos(p+psi)
            yfr = x(4,k)+disf*(sin_angl_right_front+cos_angl_right_front*x(5,k-1));% sin(p+psi)   
        end
        
        % tracking objective + slow change
        opti.subject_to(-r + u(k) <= as(k));
        opti.subject_to(r - u(k) <= as(k));
        obj = obj + 10*R*as(k)+ 10*R*(r - u(k))^2;
   
        opti.subject_to(deltaFmin <= u(k));
        opti.subject_to( u(k) <= deltaFmax);
        if k == 1
            du = u(k) - u_prev;
        else
            du = u(k) - u(k-1);
        end
        % slew protection
        opti.subject_to(-delta_slew - sD(k) <= du <= delta_slew + sD(k));
        obj = obj + Qs*sD(k)^2;
        obj = obj + (change_pen*du)^2;

        opti.subject_to(((yfl-dsy_1)^2+(xfl-dsx_1)^2) >= obs_on_1*(radius_1 - radfl1(k))^2);
        opti.subject_to(((yfr-dsy_1)^2+(xfr-dsx_1)^2) >= obs_on_1*(radius_1 - radfr1(k))^2);
        opti.subject_to(((yfr-dsy_2)^2+(xfr-dsx_2)^2) >= obs_on_2*(radius_2 - radfr2(k))^2);
        opti.subject_to(((yfl-dsy_2)^2+(xfl-dsx_2)^2) >= obs_on_2*(radius_2 - radfl2(k))^2);
        opti.subject_to(radfl1(k) >= 0);
        opti.subject_to(radfr1(k) >= 0);
        opti.subject_to(radfl2(k) >= 0);
        opti.subject_to(radfr2(k) >= 0);
        obj = obj + obs_on_1*Qo*radfl1(k)^2 + obs_on_1*Qo*radfr1(k)^2 + obs_on_2*10^2*Qo*radfl2(k)^2 + obs_on_2*10^2*Qo*radfr2(k)^2;

        % line segment between left front and right front wheels
        weight = linspace(0,1,points_count+2);
        for i = 1:points_count + 1
            new_x = xfl*weight(i+1) + xfr*(1-weight(i+1));
            new_y = yfl*weight(i+1) + yfr*(1-weight(i+1));
            opti.subject_to(((new_y-dsy_1)^2+(new_x-dsx_1)^2) >= obs_on_1*(radius_1 - bigslack(k, i))^2);
            opti.subject_to(((new_y-dsy_2)^2+(new_x-dsx_2)^2) >= obs_on_2*(radius_2 - bigslack(k, points_count + 1 + i))^2);
            
            opti.subject_to(bigslack(k, i) >= 0);
            obj = obj + obs_on_1*Qo*bigslack(k, i)^2 + obs_on_2*Qo*bigslack(k, points_count + 1 + i)^2;
        end
        
        % lane keeping
        opti.subject_to( k0L - roadPadding >= yfl - k2L*xfl^2 - k1L*xfl - sfl(k));
        opti.subject_to(k0R + roadPadding <= yfr - k2R*xfr^2 - k1R*xfr + sfr(k));
        opti.subject_to(sfl(k) >= 0);
        opti.subject_to(sfr(k) >= 0);
        opti.subject_to(srl(k) >= 0);
        opti.subject_to(srr(k) >= 0);
        % slack objective
        obj = obj + Qenv_ar*sfl(k)^2 + Qenv_ar*sfr(k)^2+ Qenv_ar*srl(k)^2 + Qenv_ar*srr(k)^2;
        obj = obj + 100*u(k)^2;
        obj = obj + 50*x(1,k)^2;
    end
%   add constrains for the final state
    opti.subject_to(sfl(N+1) >= 0);
    opti.subject_to(sfr(N+1) >= 0);
    opti.subject_to(srl(N+1) >= 0);
    opti.subject_to(srr(N+1) >= 0);
    opti.subject_to(radfl1(N+1) >= 0);
    opti.subject_to(radfr1(N+1) >= 0);
    opti.subject_to(radfl2(N+1) >= 0);
    opti.subject_to(radfr2(N+1) >= 0);
%   x(k+1)
    opti.subject_to(x(3:5, N+1) == x(3:5, N) + Ts*[v*1; v*(x(5, N) + x(1, N)); x(2)]);

    xfl = x(3,N+1)+disf*(cos_angl_left_front-sin_angl_left_front*x(5,N));% cos(p+psi)
    yfl = x(4,N+1)+disf*(sin_angl_left_front+cos_angl_left_front*x(5,N));% sin(p+psi)
    xfr = x(3,N+1)+disf*(cos_angl_right_front-sin_angl_right_front*x(5,N));% cos(p+psi)
    yfr = x(4,N+1)+disf*(sin_angl_right_front+cos_angl_right_front*x(5,N));% sin(p+psi)
            
    %obstacle avoidance and lane keeping k+1
    opti.subject_to(((yfl-dsy_1)^2+(xfl-dsx_1)^2) >= obs_on_1*(radius_1 - radfl1(N+1))^2);
    opti.subject_to(((yfr-dsy_1)^2+(xfr-dsx_1)^2) >= obs_on_1*(radius_1 - radfr1(N+1))^2);
    opti.subject_to(((yfr-dsy_2)^2+(xfr-dsx_2)^2) >= obs_on_2*(radius_2 - radfr2(N+1))^2);
    opti.subject_to(((yfl-dsy_2)^2+(xfl-dsx_2)^2) >= obs_on_2*(radius_2 - radfl2(N+1))^2);
    opti.subject_to( k0L - roadPadding >= yfl - k2L*xfl^2 - k1L*xfl - sfl(N+1));
    opti.subject_to(k0R + roadPadding <= yfr - k2R*xfr^2 - k1R*xfr + sfr(N+1));
    obj = obj + obs_on_1*Qo*radfl1(N+1)^2 + obs_on_1*Qo*radfr1(N+1)^2 + obs_on_2*10^2*Qo*radfl2(N+1)^2 + obs_on_2*10^2*Qo*radfr2(N+1)^2;
    obj = obj + Qenv_ar*sfl(N+1)^2 + Qenv_ar*sfr(N+1)^2+ Qenv_ar*srl(N+1)^2 + Qenv_ar*srr(N+1)^2;
    obj = obj + 50*x(1,N+1)^2;
    
    weight = linspace(0,1,points_count+2);
    for i = 1:points_count + 1
        new_x = xfl*weight(i+1) + xfr*(1-weight(i+1));
        new_y = yfl*weight(i+1) + yfr*(1-weight(i+1));
        opti.subject_to(((new_y-dsy_1)^2+(new_x-dsx_1)^2) >= obs_on_1*(radius_1 - bigslack(N+1, i))^2);
        opti.subject_to(((new_y-dsy_2)^2+(new_x-dsx_2)^2) >= obs_on_2*(radius_2 - bigslack(N+1, points_count + 1 + i))^2);
        obj = obj + obs_on_1*Qo*bigslack(N+1, i)^2 + obs_on_2*Qo*bigslack(N+1, points_count + 1 + i)^2;
    end
    
    opti.minimize(obj);
    more_opt = struct('print_level',0,'print_timing_statistics','no',...
    'acceptable_iter',1,...
    'acceptable_tol', 1e0,...
    'constr_viol_tol', 1e-1,...
    'derivative_test_tol', 1,...
    'fast_step_computation', 'yes',...
    'magic_steps', 'yes',...
    'tiny_step_tol', 1e0,...
    'warm_start_init_point', 'yes',...
    'bound_relax_factor', 1e-2,...
    'acceptable_constr_viol_tol', 1e6,...
    'acceptable_compl_inf_tol', 1e6,...
    'fixed_variable_treatment', 'make_constraint',...
    'mumps_pivtol', 0);
    
    s_opts = struct('verbose',false,'print_time',false,'ipopt',more_opt);
        opti.solver('ipopt',s_opts); % dont forget about solver params
end

if(controller_on ~= 1 || 0 == input(3)) % check v_min and start_time boundary
    u_mpc = input(5); % sends the reference instead
else
    
    % line fitting
    Xl = [points_x_l(1); points_x_l(2);points_x_l(3);points_x_l(4);points_x_l(5);points_x_l(6);points_x_l(7);points_x_l(8);points_x_l(9);points_x_l(10)];
    Yl = [points_y_l(1); points_y_l(2); points_y_l(3); points_y_l(4); points_y_l(5); points_y_l(6); points_y_l(7); points_y_l(8); points_y_l(9); points_y_l(10)];
    Xr = [points_x_r(1); points_x_r(2);points_x_r(3);points_x_r(4);points_x_r(5);points_x_r(6);points_x_r(7);points_x_r(8);points_x_r(9);points_x_r(10)];
    Yr = [points_y_r(1); points_y_r(2); points_y_r(3); points_y_r(4); points_y_r(5); points_y_r(6); points_y_r(7); points_y_r(8); points_y_r(9); points_y_r(10)];
    fl = polyfit(Xl,Yl,2);
    fr = polyfit(Xr,Yr,2);
    % input parsing
    opti.set_value(v, input(3));
    opti.set_value(u_prev, input(6));
    opti.set_value(r, input(5));
    
    opti.set_value(k0L, fl(3));
    opti.set_value(k1L, fl(2));
    opti.set_value(k2L, fl(1));
    
    opti.set_value(k0R, fr(3));
    opti.set_value(k1R, fr(2));
    opti.set_value(k2R, fr(1));
    
    opti.set_value(caf, params(7));
    opti.set_value(car, params(8));
 
    if obstacle_detected(1) >= 1 
        opti.set_value(radius_1,rad_ius);
        opti.set_value(dsx_1,ds_x(1));
        opti.set_value(dsy_1,ds_y(1));
        opti.set_value(obs_on_1, 1);
    else
        opti.set_value(radius_1,0);
        opti.set_value(dsx_1,0);
        opti.set_value(dsy_1,0);
        opti.set_value(obs_on_1, 0);
    end
     
    if obstacle_detected(2) >= 1
        opti.set_value(radius_2,rad_ius);
        opti.set_value(dsx_2,ds_x(2));
        opti.set_value(dsy_2,ds_y(2));
        opti.set_value(obs_on_2, 1);
    else
        opti.set_value(radius_2,0);
        opti.set_value(dsx_2,0);
        opti.set_value(dsy_2,0);
        opti.set_value(obs_on_2, 0);
    end
    
    if obstacle_detected(2) < 1 && obstacle_detected(1) >= 1
        opti.set_value(radius_1,rad_ius);
        opti.set_value(dsx_1,ds_x(1));
        opti.set_value(dsy_1,ds_y(1));
        opti.set_value(obs_on_1, 1);
        opti.set_value(radius_2,rad_ius);
        opti.set_value(dsx_2,ds_x(1));
        opti.set_value(dsy_2,ds_y(1));
        opti.set_value(obs_on_2, 1);
    end
    
    if obstacle_detected(1) < 1 && obstacle_detected(2) >= 1
        opti.set_value(radius_1,rad_ius);
        opti.set_value(dsx_1,ds_x(2));
        opti.set_value(dsy_1,ds_y(2));
        opti.set_value(obs_on_1, 1);
        opti.set_value(radius_2,rad_ius);
        opti.set_value(dsx_2,ds_x(2));
        opti.set_value(dsy_2,ds_y(2));
        opti.set_value(obs_on_2, 1);
    end
    
    if obstacle_detected(2) < 1 && obstacle_detected(1) < 1
        opti.set_value(radius_1,0);
        opti.set_value(dsx_1,0);
        opti.set_value(dsy_1,0);
        opti.set_value(obs_on_1, 0);
        opti.set_value(radius_2,0);
        opti.set_value(dsx_2,0);
        opti.set_value(dsy_2,0);
        opti.set_value(obs_on_2, 0);
    end

    opti.set_initial(u(:), input(6));
    opti.set_value(x0, [input(1), input(2), 0, 0, 0]);
    
    %% obtain a solution
    try
       sol = opti.solve();
    catch
       deb = opti.debug;
    end
    u_mpc = sol.value(u(1));
end
