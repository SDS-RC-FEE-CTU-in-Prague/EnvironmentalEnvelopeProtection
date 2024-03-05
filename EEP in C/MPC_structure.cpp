//
// Created by efrem on 11/1/2023.
//
#include <casadi/casadi.hpp>
#include <list>
#include "Obstacle_structure.cpp"

class MPC_structure {
    // ok, looking at the documentation https://web.casadi.org/docs/, probably all this code can be presented a bit easier
public:
    // mpc stuff
    unsigned int nx = 5; // number of states
    unsigned int nu = 1; // number of inputs
    unsigned int max_obstacle_count;
    //pred. horizon
    unsigned int N;
    casadi::Opti opti = casadi::Opti();
    // basic optimization vars
    casadi::MX r = opti.parameter(nu, 1); // references for inputs
    casadi::MX x0 = opti.parameter(nx, 1); // initial state
    casadi::MX x; // state variables
    casadi::MX u; // input variables
    casadi::MX as; // parameter for abs substitution
    casadi::MX u_prev = opti.parameter(nu, 1); // previous input
    // borders curve param
    casadi::MX k3L = opti.parameter(2);
    casadi::MX k2L = opti.parameter(2); // second parameter is used for rut road
    casadi::MX k1L = opti.parameter(2);
    casadi::MX k0L = opti.parameter(2);
    casadi::MX k3R = opti.parameter(2);
    casadi::MX k2R = opti.parameter(2);
    casadi::MX k1R = opti.parameter(2);
    casadi::MX k0R = opti.parameter(2);
    // slack vars for lines
    casadi::MX sfl; //the second row is used for rut roads
    casadi::MX sfr;
    // radius of obstacle for each front wheel
    casadi::MX radius; //will be defined in the constructor
    casadi::MX obs_weights; // weighting factors (take as priorities) for each obstacle
    // the following parameters can be used for additional weight changing
    casadi::MX obs_on; // each obstacle has it own instance: 0 - no obstacle, 1 - obstacle
    casadi::MX undrivable_on; // each obstacle has it own instance: 0 - drivable, 1 - undrivable (we add additional constraint)
    // distance to obstacle for each wheel
    casadi::MX dsx;
    casadi::MX dsy;
    // slack vars for obstacles
    casadi::MX sofl; // left front wheel
    casadi::MX sofr; // right front wheel
    casadi::MX somid; // middle part of the front axle (could be extended for more points in the future)
    // slack for fast input change
    casadi::MX sdu;
    // additional variables
    casadi::MX v = opti.parameter(); // velocity at the moment
    casadi::MX caf = opti.parameter(); // tire stiffness
    casadi::MX car = opti.parameter();

    MPC_structure(unsigned int prediction_horizon, unsigned int max_obstacle_count) {
        this->N = prediction_horizon;
        this->x = opti.variable(nx, N + 1);
        this->u = opti.variable(nu, N);
        this->as = opti.variable(nu, N + 1);
        this->sfl = opti.variable(N + 1, 2); //the second row is used for rut roads
        this->sfr = opti.variable(N + 1, 2);
        this->sdu = opti.variable(N, 1);
        this->max_obstacle_count = max_obstacle_count;
        this->radius = opti.parameter(max_obstacle_count);
        this->obs_weights = opti.parameter(max_obstacle_count);
        this->obs_on = opti.parameter(max_obstacle_count);
        this->undrivable_on = opti.parameter(max_obstacle_count);
        this->dsx = opti.parameter(max_obstacle_count);
        this->dsy = opti.parameter(max_obstacle_count);
        this->sofl = opti.variable(N + 1, max_obstacle_count);
        this->sofr = opti.variable(N + 1, max_obstacle_count);
        this->somid = opti.variable(N + 1, max_obstacle_count);
    }

    // vehParams: lf, lr, m, I, h
    // weights: R, Qs, Qenv_ar, Qdu, Qu, Qx (only for beta)
    void initMpc(const double vehParams[5], const double *weights, double roadPadding) {
        // vehicle parameters parsing
        double lf = vehParams[0]; // GC-front distance
        double lr = vehParams[1]; // GC-rear distance
        double m = vehParams[2]; // mass
        double I = vehParams[3]; // inertia
        double ha = vehParams[4]; // half of the axle
        double wheel_dist_left = ha;   // distance to a wheel on y-axis
        double wheel_dist_right = -ha;
        double Ts = 0.05; // sampling time
        double deltaFmax = 0.65;
        double deltaFmin = -0.65;
        // slew restrictions
        double delta_slew = 240.0 / 360.0 * 2.0 * 3.1416 * Ts;

        // weights parsing
        double R = weights[0]; //steering wheel tracking
        double Qs = weights[1]; // slack weight
        double Qenv_ar = weights[2];
        double Qdu = weights[3];
        double Qu = weights[4];
        double Qx = weights[5]; // used only for beta

        casadi::MX obj = 0; //main objective

        //general vehicle dynamics
        casadi::MX a11 = 1 + Ts * -(caf + car) / m / v;
        casadi::MX a12 = Ts * ((car * lr - caf * lf) / m / v / v - 1);
        casadi::MX a21 = Ts * (car * lr - caf * lf) / I;
        casadi::MX a22 = 1 + Ts * -(lf * lf * caf + car * lr * lr) / I / v;
        casadi::MX b1 = Ts * caf / m / v;
        casadi::MX b2 = Ts * caf * lf / I;

        //additional values for wheel positions calculation
        double disf = sqrt(pow(lf, 2) + pow(wheel_dist_left, 2));
        double cos_angl_left_front = cos(atan2(wheel_dist_left, lf));
        double sin_angl_left_front = sin(atan2(wheel_dist_left, lf));
        double cos_angl_right_front = cos(atan2(wheel_dist_right, lf));
        double sin_angl_right_front = sin(atan2(wheel_dist_right, lf));

        //local variables
        casadi::MX du = opti.variable();
        casadi::MX xfl = opti.variable();
        casadi::MX yfl = opti.variable();
        casadi::MX xfr = opti.variable();
        casadi::MX yfr = opti.variable();

        for (unsigned int k = 0;
             k < N + 1; ++k) { //it goes also for N step due to overlapping for persistent feasibility
            // basic dynamics and inputs
            if (k == 0) {
                // dynamic model beta, r, x, y, psi
                opti.subject_to(x(0, k) == a11 * x0(0) + a12 * x0(1) + b1 * u(k));
                opti.subject_to(x(1, k) == a21 * x0(0) + a22 * x0(1) + b2 * u(k));
                opti.subject_to(x(2, k) == Ts * v);
                opti.subject_to(x(3, k) == Ts * v * x0(0));
                opti.subject_to(x(4, k) == Ts * x0(1));
            } else if (k == N) {
                //   x(k+1) assume input equal to the previous step (for persistent feasibility)
                opti.subject_to(x(0, k) == a11 * x(0, k - 1) + a12 * x(1, k - 1) + b1 * u(k - 1));
                opti.subject_to(x(1, k) == a21 * x(0, k - 1) + a22 * x(1, k - 1) + b2 * u(k - 1));
                opti.subject_to(x(2, k) == x(2, k - 1) + Ts * v);
                opti.subject_to(x(3, k) == x(3, k - 1) + Ts * v * (x(4, k - 1) + x(0, k - 1)));
                opti.subject_to(x(4, k) == x(4, k - 1) + Ts * x(1));
            } else {
                opti.subject_to(x(0, k) == a11 * x(0, k - 1) + a12 * x(1, k - 1) + b1 * u(k));
                opti.subject_to(x(1, k) == a21 * x(0, k - 1) + a22 * x(1, k - 1) + b2 * u(k));
                opti.subject_to(x(2, k) == x(2, k - 1) + Ts * v);
                opti.subject_to(x(3, k) == x(3, k - 1) + Ts * v * (x(4, k - 1) + x(0, k - 1)));
                opti.subject_to(x(4, k) == x(4, k - 1) + Ts * x(1));
            }
            // wheel positions
            if (k == 0) {
                xfl = x(2, k) + disf * cos_angl_left_front;
                yfl = x(3, k) + disf * sin_angl_left_front;
                xfr = x(2, k) + disf * cos_angl_right_front;
                yfr = x(3, k) + disf * sin_angl_right_front;
            } else {
                xfl = x(2, k) + disf * (cos_angl_left_front - sin_angl_left_front * x(4, k - 1));
                yfl = x(3, k) + disf * (sin_angl_left_front + cos_angl_left_front * x(4, k - 1));
                xfr = x(2, k) + disf * (cos_angl_right_front - sin_angl_right_front * x(4, k - 1));
                yfr = x(3, k) + disf * (sin_angl_right_front + cos_angl_right_front * x(4, k - 1));
            }
            // tracking objective
            if (k != N) { // there is no change for the last iteration
                opti.subject_to(-r + u(k) <= as(k));
                opti.subject_to(r - u(k) <= as(k));
                opti.subject_to(as(k) >= 0);
                obj = obj + R * as(k) + R * (r - u(k)) * (r - u(k));
            }
            // input constraints and penalisation
            if (k != N) {
                opti.subject_to(deltaFmin <= u(k));
                opti.subject_to(u(k) <= deltaFmax);
                if (k == 0) {
                    du = u(k) - u_prev;
                } else {
                    du = u(k) - u(k - 1);
                }
                // slew protection
                opti.subject_to(-delta_slew - sdu(k) <= du <= delta_slew + sdu(k));
                obj = obj + Qs * sdu(k) * sdu(k);
                obj = obj + (Qdu * du) * (Qdu * du);
                obj = obj + Qu * u(k) * u(k);
                opti.subject_to(sdu(k) >= 0);
            }

            // obstacle avoidance
            for (unsigned int i = 0; i < max_obstacle_count; ++i) { //slacks go here
                opti.subject_to(sofl(k, i) >= 0);
                opti.subject_to(sofr(k, i) >= 0);
                opti.subject_to(somid(k, i) >= 0);
                obj = obj + obs_weights(i) * obs_on(i) * sofl(k, i) * sofl(k, i);
                obj = obj + obs_weights(i) * obs_on(i) * sofr(k, i) * sofr(k, i);
                obj = obj + obs_weights(i) * obs_on(i) * undrivable_on(i) * somid(k, i) * somid(k, i);

                opti.subject_to(
                        (pow(yfl - dsy(i), 2) + pow(xfl - dsx(i), 2)) >= obs_on(i) * pow(radius(i) - sofl(k, i), 2));
                opti.subject_to(
                        (pow(yfr - dsy(i), 2) + pow(xfr - dsx(i), 2)) >= obs_on(i) * pow(radius(i) - sofr(k, i), 2));
                opti.subject_to(
                        (pow((yfl + yfr) / 2 - dsy(i), 2) + pow((xfl + xfr) / 2 - dsx(i), 2)) >=
                        obs_on(i) * undrivable_on(i) * pow(radius(i) - somid(k, i), 2));
            }

            // lane keeping
            // 0 is always left side of a wheel (!), 1 - right
            opti.subject_to(sfl(k, 0) >= 0); // left wheel, left line
            opti.subject_to(sfr(k, 0) >= 0); // right wheel, right line
            opti.subject_to(sfl(k, 1) >= 0); // left wheel, right line
            opti.subject_to(sfr(k, 1) >= 0); // right wheel, left line
            obj = obj + Qenv_ar * sfl(k, 0) * sfl(k, 0) + Qenv_ar * sfr(k, 0) * sfr(k, 0);
            obj = obj + Qenv_ar * sfl(k, 1) * sfl(k, 1) + Qenv_ar * sfr(k, 1) * sfr(k, 1);

            opti.subject_to(k0L(0) - roadPadding >=
                            yfl - k3L(0) * xfl * xfl * xfl - k2L(0) * xfl * xfl - k1L(0) * xfl - sfl(k, 0));
//            //the following row is used only for road rut, and can be switched off when k2L, k0L -> -Inf, k1L -> 0
            opti.subject_to(k0L(1) + roadPadding <=
                            yfl - k3L(1) * xfl * xfl * xfl - k2L(1) * xfl * xfl - k1L(1) * xfl + sfl(k, 1));
//
            opti.subject_to(k0R(0) + roadPadding <=
                            yfr - k3R(0) * xfr * xfr * xfr - k2R(0) * xfr * xfr - k1R(0) * xfr + sfr(k, 0));
//            //the following row is used only for road rut, and can be switched off when k2R, k0R -> +Inf, k1R -> 0
            opti.subject_to(k0R(1) - roadPadding >=
                            yfr - k3R(1) * xfr * xfr * xfr - k2R(1) * xfr * xfr - k1R(1) * xfr - sfr(k, 1));

            // this is just a penalization of fast changes of the sideslip
            obj = obj + Qx * x(0, k) * x(0, k);
        }

        // final function
        opti.minimize(obj);

        // solver params
        casadi::Dict opts_dict = casadi::Dict();
        opts_dict["ipopt.print_level"] = 0;
        opts_dict["ipopt.print_timing_statistics"] = "no";
        opts_dict["ipopt.acceptable_iter"] = 1;
        opts_dict["ipopt.acceptable_tol"] = 2e0;
        opts_dict["ipopt.constr_viol_tol"] = 1e-1;
        opts_dict["ipopt.derivative_test_tol"] = 1;
        opts_dict["ipopt.fast_step_computation"] = "yes";
        opts_dict["ipopt.magic_steps"] = "yes";
        opts_dict["ipopt.tiny_step_tol"] = 1e0;
        opts_dict["ipopt.warm_start_init_point"] = "yes";
        opts_dict["ipopt.bound_relax_factor"] = 1e-2;
        opts_dict["ipopt.acceptable_constr_viol_tol"] = 1e6; //1e6
        opts_dict["ipopt.acceptable_compl_inf_tol"] = 1e6; // 1e6
        opts_dict["ipopt.fixed_variable_treatment"] = "make_parameter_nodual";
        opts_dict["ipopt.expect_infeasible_problem"] = "no";
        opts_dict["ipopt.linear_solver"] = "ma27";
        opts_dict["ipopt.nlp_scaling_method"] = "gradient-based";
        opts_dict["ipopt.linear_scaling_on_demand"] = "yes";
        opts_dict["ipopt.least_square_init_primal"] = "yes";
        opts_dict["ipopt.nlp_lower_bound_inf"] = -1e12;
        opts_dict["ipopt.nlp_upper_bound_inf"] = 1e12;
        opts_dict["ipopt.mu_strategy"] = "adaptive";
        opts_dict["ipopt.adaptive_mu_globalization"] = "kkt-error";
        opts_dict["ipopt.mu_oracle"] = "loqo";
        opts_dict["ipopt.tol"] = 2e0;
        opts_dict["ipopt.dual_inf_tol"] = 1e1;
        opts_dict["ipopt.compl_inf_tol"] = 1e6;
        opts_dict["ipopt.mu_target"] = 1e-1;
        opts_dict["ipopt.accept_after_max_steps"] = 10;
        opts_dict["ipopt.adaptive_mu_kkterror_red_iters"] = 2;
        opts_dict["verbose"] = 0;
        opts_dict["print_time"] = 0;
        opts_dict["expand"] = 1;


        opti.solver("ipopt", opts_dict);
    }

    double solveMPC(double beta_measured, double r_measured, double v_measured, double df_ref, double df_prev,
                    const double params[2], const double roadParams[12], bool rutRoadOn,
                    std::list<Obstacle_structure> obstacles) {

        // input parsing
        opti.set_value(v, v_measured);
        opti.set_value(u_prev, df_prev);
        opti.set_value(r, df_ref);
        opti.set_value(caf, params[0]);
        opti.set_value(car, params[1]);
        opti.set_initial(u, df_ref);
        opti.set_value(x0, {beta_measured, r_measured, 0.0, 0.0, 0.0});

        // road params setting
        opti.set_value(k0L(0), roadParams[0]);
        opti.set_value(k1L(0), roadParams[1]);
        opti.set_value(k2L(0), roadParams[2]);
        opti.set_value(k3L(0), roadParams[3]);
        if (rutRoadOn) {
            opti.set_value(k0L(1), roadParams[4]);
            opti.set_value(k1L(1), roadParams[5]);
            opti.set_value(k2L(1), roadParams[6]);
            opti.set_value(k3L(1), roadParams[7]);
        } else { // set the right boundary for the right wheel (unreachable)
            opti.set_value(k0L(1), -100); // assume as -Inf (maybe we can use something from casadi)
            opti.set_value(k1L(1), 0);
            opti.set_value(k2L(1), 0); // assume as -Inf
            opti.set_value(k3L(1), 0);
        }

        opti.set_value(k0R(0), roadParams[8]);
        opti.set_value(k1R(0), roadParams[9]);
        opti.set_value(k2R(0), roadParams[10]);
        opti.set_value(k3R(0), roadParams[11]);
        if (rutRoadOn) {
            opti.set_value(k0R(1), roadParams[12]);
            opti.set_value(k1R(1), roadParams[13]);
            opti.set_value(k2R(1), roadParams[14]);
            opti.set_value(k3R(1), roadParams[15]);
        } else { // set the left boundary for the left wheel (unreachable)
            opti.set_value(k0R(1), 100); // assume as +Inf
            opti.set_value(k1R(1), 0);
            opti.set_value(k2R(1), 0); // assume as +Inf
            opti.set_value(k3R(1), 0);
        }

        // obstacles params settings
        std::vector<double> radius_vec;
        std::vector<double> obs_weights_vec;
        std::vector<double> obs_on_vec;
        std::vector<double> undrivable_on_vec;
        std::vector<double> dsx_vec;
        std::vector<double> dsy_vec;

        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            radius_vec.push_back(it->radius);
            obs_weights_vec.push_back(it->weight);
            obs_on_vec.push_back(it->obs_on);
            undrivable_on_vec.push_back(it->undrivable_on);
            dsx_vec.push_back(it->dsx);
            dsy_vec.push_back(it->dsy);
        }

        opti.set_value(radius, radius_vec);
        opti.set_value(obs_weights, obs_weights_vec);
        opti.set_value(obs_on, obs_on_vec);
        opti.set_value(undrivable_on, undrivable_on_vec);
        opti.set_value(dsx, dsx_vec);
        opti.set_value(dsy, dsy_vec);

        // obtain a solution
        double u_mpc = df_prev;
        try {
            auto sol = opti.solve();
            auto tmp = sol.value(u);
            u_mpc = tmp->front();
        }
        catch (int exNumber) {
            auto deb = opti.debug();
        }

        return u_mpc;
    }
};