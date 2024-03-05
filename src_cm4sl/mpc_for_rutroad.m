function [u_mpc] = mpc_for_rutroad(t, input, params, controller_on, line_detection)

%           1   2  3  4    5       6       7      8
% input = beta, r, v, yaw, df_ref, u_prev, x_cor, y_cor
persistent opti r x0 u_prev v x u as k0LL k1LL k2LL k3LL k0LR k1LR k2LR k3LR k0RL k1RL k2RL k3RL k0RR k1RR k2RR k3RR sfrl sfll srrl srll sfrr sflr srrr srlr sol sD caf car

% vehicle parameters
lf   = params(1); % GC-front distance
lr   = params(2); % GC-rear distance
m    = params(3); % mass
I    = params(4); % inertia
ha   = params(5); % half of the axle
wheel_dist_left = ha;   % distance to a wheel on y-axis
wheel_dist_right = -ha;

if t == 0
    opti = casadi.Opti(); % definition of the controller object
    roadPadding = 0.2;
    % mpc stuff
    nx = 5; % number of states
    nu = 1; % number of inputs
    
     % pred. horizon
    N = 10;
    % basic oprimization vars 
    r = opti.parameter(nu, 1); % references for inputs
    x0 = opti.parameter(nx, 1); % initial state
    x = opti.variable(nx, N+1); % state variables
    u = opti.variable(nu, N); % input variables
    as = opti.variable(nu, N+1); % parameter for abs substitution
    u_prev = opti.parameter(nu, 1); % previous input
    % bordes curve param
    k3LL = opti.parameter(); 
    k2LL = opti.parameter(); 
    k1LL = opti.parameter();
    k0LL = opti.parameter();
    k3LR = opti.parameter(); 
    k2LR = opti.parameter(); 
    k1LR = opti.parameter();
    k0LR = opti.parameter();
    k3RL = opti.parameter();
    k2RL = opti.parameter();
    k1RL = opti.parameter();
    k0RL = opti.parameter();
    k3RR = opti.parameter();
    k2RR = opti.parameter();
    k1RR = opti.parameter();
    k0RR = opti.parameter();
    % slack vars
    sfll = opti.variable(N+1, 1);
    sflr = opti.variable(N+1, 1);
    sfrl = opti.variable(N+1, 1);
    sfrr = opti.variable(N+1, 1);
    srll = opti.variable(N+1, 1);
    srlr = opti.variable(N+1, 1);
    srrl = opti.variable(N+1, 1);
    srrr = opti.variable(N+1, 1);
    sD = opti.variable(N, 1);

    % weights
    steering_wheel_track = 1e3;
    R = steering_wheel_track;
    Qs   = 1e10; % slack weight 
    Qenv_ar =1e4;
    change_pen=1e2;
    v = opti.parameter(); % velocity at the moment
    caf = opti.parameter();
    car = opti.parameter();
    Ts = 0.05;%preview_dist/v/N;
    deltaFmax = 0.65;
    deltaFmin = -0.65;
    % slew restrictions
    delta_slew = deg2rad(120)*Ts; % rad/s
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
        obj = obj + R*as(k)+ R*(r - u(k))^2;
        
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
        % lane keeping
            
        opti.subject_to( k0LL - roadPadding >= yfl - k3LL*xfl^3 - k2LL*xfl^2 - k1LL*xfl - sfll(k));
        opti.subject_to( k0LR + roadPadding <= yfl - k3LR*xfl^3 - k2LR*xfl^2 - k1LR*xfl + sflr(k));
        opti.subject_to(k0RL - roadPadding >= yfr - k3RL*xfr^3 - k2RL*xfr^2 - k1RL*xfr - sfrl(k));
        opti.subject_to(k0RR + roadPadding <= yfr - k3RR*xfr^3 - k2RR*xfr^2 - k1RR*xfr + sfrr(k));
        % slack stuff
        opti.subject_to(sfll(k) >= 0);
        opti.subject_to(sfrl(k) >= 0);
        opti.subject_to(srll(k) >= 0);
        opti.subject_to(srrl(k) >= 0);
        opti.subject_to(sflr(k) >= 0);
        opti.subject_to(sfrr(k) >= 0);
        opti.subject_to(srlr(k) >= 0);
        opti.subject_to(srrr(k) >= 0);
        % slack objective
        obj = obj + Qenv_ar*sfll(k)^2 + Qenv_ar*sfrl(k)^2;
        obj = obj + Qenv_ar*sflr(k)^2 + Qenv_ar*sfrr(k)^2;
        obj = obj + 100*u(k)^2;
        obj = obj + 50*x(1,k)^2;
    end
%  % add constrains for the final state
    opti.subject_to(sfll(k) >= 0);
    opti.subject_to(sfrl(k) >= 0);
    opti.subject_to(sflr(k) >= 0);
    opti.subject_to(sfrr(k) >= 0);
%   x(k+1)
    opti.subject_to(x(3:5, N+1) == x(3:5, N) + Ts*[v*1; v*(x(5, N) + x(1, N)); x(2)]);

    xfl = x(3,N+1)+disf*(cos_angl_left_front-sin_angl_left_front*x(5,N));% cos(p+psi)
    yfl = x(4,N+1)+disf*(sin_angl_left_front+cos_angl_left_front*x(5,N));% sin(p+psi)
    xfr = x(3,N+1)+disf*(cos_angl_right_front-sin_angl_right_front*x(5,N));% cos(p+psi)
    yfr = x(4,N+1)+disf*(sin_angl_right_front+cos_angl_right_front*x(5,N));% sin(p+psi)
            
    %obstacle avoidance and lane keeping k+1
    opti.subject_to( k0LL - roadPadding >= yfl - k3LL*xfl^3 - k2LL*xfl^2 - k1LL*xfl - sfll(N+1));
    opti.subject_to( k0LR + roadPadding <= yfl - k3LR*xfl^3 - k2LR*xfl^2 - k1LR*xfl + sflr(N+1));
    opti.subject_to(k0RL - roadPadding >= yfr - k3RL*xfr^3 - k2RL*xfr^2 - k1RL*xfr - sfrl(N+1));
    opti.subject_to(k0RR + roadPadding <= yfr - k3RR*xfr^3 - k2RR*xfr^2 - k1RR*xfr + sfrr(N+1));
    obj = obj + Qenv_ar*sfll(N+1)^2 + Qenv_ar*sfrl(N+1)^2;
    obj = obj + Qenv_ar*sflr(N+1)^2 + Qenv_ar*sfrr(N+1)^2;
    obj = obj + 50*x(1,N+1)^2;

    
    opti.minimize(obj);

    more_opt = struct('print_level', 0, ...
        'print_timing_statistics', 'no', ...
        'acceptable_iter', 1, ...
        'acceptable_tol', 2e0, ...
        'constr_viol_tol', 1e-1, ...
        'derivative_test_tol', 1, ...
        'fast_step_computation', 'yes', ...
        'magic_steps', 'yes', ...
        'tiny_step_tol', 1e0, ...
        'warm_start_init_point', 'yes', ...
        'bound_relax_factor', 1e-2, ...
        'acceptable_constr_viol_tol', 1e6, ...
        'acceptable_compl_inf_tol', 1e6, ...
        'fixed_variable_treatment', 'make_parameter_nodual', ...
        'expect_infeasible_problem', 'no', ...
        'linear_solver', 'ma27', ...
        'nlp_scaling_method', 'gradient-based', ...
        'linear_scaling_on_demand', 'yes', ...
        'least_square_init_primal', 'yes', ...
        'nlp_lower_bound_inf', -1e12, ...
        'nlp_upper_bound_inf', 1e12,...
        'mu_strategy', 'adaptive',...
        'adaptive_mu_globalization', 'kkt-error',...
        'mu_oracle', 'loqo',...
        'tol', 2e0, ...
        'dual_inf_tol', 1e1, ...
        'compl_inf_tol', 1e6, ...
        'mu_target', 1e-1,...
        'accept_after_max_steps', 10,...
        'adaptive_mu_kkterror_red_iters', 2);
    
    s_opts = struct('verbose',false,'print_time',false,'expand',true,'ipopt',more_opt);
    opti.solver('ipopt',s_opts); % dont forget about solver params
end

if(controller_on ~= 1 || 0 == input(3)) % check v_min and start_time boundary
    u_mpc = input(5);
else
    
    %line fitting
    main_left = polyfit(line_detection(:,1),line_detection(:,2),3);
    wheel_left_left = polyfit(line_detection(:,3),line_detection(:,4),3);
    wheel_left_right = polyfit(line_detection(:,5),line_detection(:,6),3);
    wheel_right_left = polyfit(line_detection(:,7),line_detection(:,8),3);
    wheel_right_right = polyfit(line_detection(:,9),line_detection(:,10),3);
    main_right = polyfit(line_detection(:,11),line_detection(:,12),3);
    % input parsing
    opti.set_value(v, input(3));
    opti.set_value(u_prev, input(6));
    opti.set_value(r, input(5));
    
    % here it is changeble to go IN and OUT of ruts
    
    % In ruts
    opti.set_value(k0LL, wheel_left_left(4));
    opti.set_value(k1LL, wheel_left_left(3));
    opti.set_value(k2LL, wheel_left_left(2));
    opti.set_value(k3LL, wheel_left_left(1));
    opti.set_value(k0LR, wheel_left_right(4));
    opti.set_value(k1LR, wheel_left_right(3));
    opti.set_value(k2LR, wheel_left_right(2));
    opti.set_value(k3LR, wheel_left_right(1));
    
    opti.set_value(k0RL, wheel_right_left(4));
    opti.set_value(k1RL, wheel_right_left(3));
    opti.set_value(k2RL, wheel_right_left(2));
    opti.set_value(k3RL, wheel_right_left(1));
    opti.set_value(k0RR, wheel_right_right(4));
    opti.set_value(k1RR, wheel_right_right(3));
    opti.set_value(k2RR, wheel_right_right(2));
    opti.set_value(k3RR, wheel_right_right(1));
    
    % out ruts
%     opti.set_value(k0LL, main_left(4));
%     opti.set_value(k1LL, main_left(3));
%     opti.set_value(k2LL, main_left(2));
%     opti.set_value(k3LL, main_left(1));
%     opti.set_value(k0LR, wheel_left_left(4));
%     opti.set_value(k1LR, wheel_left_left(3));
%     opti.set_value(k2LR, wheel_left_left(2));
%     opti.set_value(k3LR, wheel_left_left(1));
%     
%     opti.set_value(k0RL, wheel_left_right(4));
%     opti.set_value(k1RL, wheel_left_right(3));
%     opti.set_value(k2RL, wheel_left_right(2));
%     opti.set_value(k3RL, wheel_left_right(1));
%     opti.set_value(k0RR, wheel_right_left(4));
%     opti.set_value(k1RR, wheel_right_left(3));
%     opti.set_value(k2RR, wheel_right_left(2));
%     opti.set_value(k3RR, wheel_right_left(1));
%     
    opti.set_value(caf, params(7));
    opti.set_value(car, params(8));

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

