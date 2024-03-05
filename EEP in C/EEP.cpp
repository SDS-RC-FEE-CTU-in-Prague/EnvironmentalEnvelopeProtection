//
// Created by efrem on 10/23/2023.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <casadi/casadi.hpp>
#include <chrono>
#include "environmental_envelope_math_operations.cpp"
#include "MPC_structure.cpp"

using namespace std::chrono;

//todo: rewrite comments
void parseCsvDataFile(const std::string &csvLine,
                      double *beta, double *r, double *v, double *df_ref, double *df_prev, double params[2],
                      double roadParams[16], double obstacle_detected[2], double ds_x[2], double ds_y[2]) {
    // K = number of iterations. The rest is the description from matlab
    // t - time [K x 1] (is used only for the first problem initialization),
    // input - physical inputs [K x 11]: beta, r, v, wf (not used), wr (not used), phi (not used), df_ref, wf_ref (not used), df_prev, x (not used), y (not used),
    // params - vehicle parameters [K x 8]: lf, lr, m, Iz, h, R, caf, car,
    // controller_on - logical controller on\off [K x 1],
    // points_x_r - ten x positions of right lane boundary wrt the vehicle initial position (0 in each new run of mpc) [K x 10],
    // points_x_l - ten x positions of left lane boundary wrt the vehicle initial position (0 in each new run of mpc) [K x 10],
    // points_y_r - ten y positions of right lane boundary wrt the vehicle initial position (0 in each new run of mpc) [K x 10],
    // points_y_l - ten y positions of left lane boundary wrt the vehicle initial position (0 in each new run of mpc) [K x 10],
    // obstacle_detected - flag that shows that wheel sensors detected an obstacle [K x 4], last two are unused,
    // ds_x - x positions of obstacle center wrt the vehicle initial position [K x 4], last two are unused,
    // ds_y - y positions of obstacle center wrt the vehicle initial position [K x 4], last two are unused

    std::string val;
    std::vector<double> parsed_line;
    std::stringstream s(csvLine);
    while (getline(s, val, ','))
        parsed_line.push_back(std::stod(val));

    // data parsing
    beta[0] = parsed_line[1];
    r[0] = parsed_line[2];
    v[0] = parsed_line[3];
    df_ref[0] = parsed_line[5];
    df_prev[0] = parsed_line[6];

    // caf, car
    params[0] = parsed_line[15];
    params[1] = parsed_line[16];

    // left line fitting and parameters setting
    std::vector<double> points_x;
    for (int i = 0; i < 10; i++)
        points_x.push_back(parsed_line[27 + i]);
    std::vector<double> points_y;
    for (int i = 0; i < 10; i++)
        points_y.push_back(parsed_line[47 + i]);
    calculateCubicParabolicCoefficients(points_x, points_y, roadParams[0], roadParams[1], roadParams[2], roadParams[3]);

    // right line fitting and parameters setting
    points_x.clear();
    points_y.clear();
    for (int i = 0; i < 10; i++)
        points_x.push_back(parsed_line[17 + i]);
    for (int i = 0; i < 10; i++)
        points_y.push_back(parsed_line[37 + i]);
    calculateCubicParabolicCoefficients(points_x, points_y, roadParams[8], roadParams[9], roadParams[10],
                                        roadParams[11]);

    // this is still ridiculous conversion due to stupid carmaker sensors
    if (parsed_line[57] >= 1) {
        obstacle_detected[0] = 1;
        ds_x[0] = parsed_line[61];
        ds_y[0] = parsed_line[65];
    } else {
        obstacle_detected[0] = 0;
        ds_x[0] = 0;
        ds_y[0] = 0;
    }

    if (parsed_line[58] >= 1) {
        obstacle_detected[1] = 1;
        ds_x[1] = parsed_line[62];
        ds_y[1] = parsed_line[66];
    } else {
        obstacle_detected[1] = 0;
        ds_x[1] = 0;
        ds_y[1] = 0;
    }

    if (parsed_line[58] < 1 && parsed_line[57] >= 1) {
        obstacle_detected[0] = 1;
        obstacle_detected[1] = 1;
        ds_x[0] = parsed_line[61];
        ds_x[1] = parsed_line[61];
        ds_y[0] = parsed_line[65];
        ds_y[1] = parsed_line[65];
    }

    if (parsed_line[57] < 1 && parsed_line[58] >= 1) {
        obstacle_detected[0] = 1;
        obstacle_detected[1] = 1;
        ds_x[0] = parsed_line[62];
        ds_x[1] = parsed_line[62];
        ds_y[0] = parsed_line[66];
        ds_y[1] = parsed_line[66];
    }

    if (parsed_line[58] < 1 && parsed_line[57] < 1) {
        obstacle_detected[0] = 0;
        obstacle_detected[1] = 0;
        ds_x[0] = 0;
        ds_x[1] = 0;
        ds_y[0] = 0;
        ds_y[1] = 0;
    }
}

void parseCsvCurveParams(const std::string &csvLine, double roadParams[16]) {
    std::string val;
    std::vector<double> parsed_line;
    std::stringstream s(csvLine);
    while (getline(s, val, ','))
        parsed_line.push_back(std::stod(val));

    // ll, lr, rl, rr
    roadParams[0] = parsed_line[0];
    roadParams[1] = parsed_line[1];
    roadParams[2] = parsed_line[2];
    roadParams[3] = parsed_line[3];

    roadParams[4] = parsed_line[4];
    roadParams[5] = parsed_line[5];
    roadParams[6] = parsed_line[6];
    roadParams[7] = parsed_line[7];

    roadParams[8] = parsed_line[12];
    roadParams[9] = parsed_line[13];
    roadParams[10] = parsed_line[14];
    roadParams[11] = parsed_line[15];

    roadParams[12] = parsed_line[8];
    roadParams[13] = parsed_line[9];
    roadParams[14] = parsed_line[10];
    roadParams[15] = parsed_line[11];
}

// dx and dy order: small caution; big caution
void parseCsvObstacles(const std::string &csvLine, double dx[2], double dy[2], double obstacle_detected[2]) {
    std::string val;
    std::vector<double> parsed_line;
    std::stringstream s(csvLine);
    while (getline(s, val, ','))
        parsed_line.push_back(std::stod(val));

    if (parsed_line[0] > 0) // small caution sign
    {
        obstacle_detected[0] = 1;
        dx[0] = parsed_line[0];
        dy[0] = parsed_line[1];
    } else {
        obstacle_detected[0] = 0;
        dx[0] = 0;
        dy[0] = 0;
    }

    if (parsed_line[2] > 0) // small caution sign
    {
        obstacle_detected[1] = 1;
        dx[1] = parsed_line[2];
        dy[1] = parsed_line[3];
    } else {
        obstacle_detected[1] = 0;
        dx[1] = 0;
        dy[1] = 0;
    }

}

void run_do30_p1_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_DO30_p1.csv");
    std::ofstream u_file(datafolder + "/u_DO30_p1.csv");
    std::ofstream t_file(datafolder + "/times_DO30_p1.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Brick 1", 0, 0, 0.7, 0, 0, 1e3);
    obstacles.emplace_back("Brick 2", 0, 0, 0.7, 0, 0, 1e3);

    std::string line;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

void run_do30_p2_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_DO30_p2.csv");
    std::ofstream u_file(datafolder + "/u_DO30_p2.csv");
    std::ofstream t_file(datafolder + "/times_DO30_p2.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Brick 1", 0, 0, 0.7, 0, 0, 1e3);
    obstacles.emplace_back("Brick 2", 0, 0, 0.7, 0, 0, 1e3);

    std::string line;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

void run_do70_p1_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_DO70_p1.csv");
    std::ofstream u_file(datafolder + "/u_DO70_p1.csv");
    std::ofstream t_file(datafolder + "/times_DO70_p1.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Brick 1", 0, 0, 0.7, 0, 0, 1e3);
    obstacles.emplace_back("Brick 2", 0, 0, 0.7, 0, 0, 1e3);

    std::string line;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

void run_do70_p2_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_DO70_p2.csv");
    std::ofstream u_file(datafolder + "/u_DO70_p2.csv");
    std::ofstream t_file(datafolder + "/times_DO70_p2.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Brick 1", 0, 0, 0.7, 0, 0, 1e3);
    obstacles.emplace_back("Brick 2", 0, 0, 0.7, 0, 0, 1e3);

    std::string line;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

void run_op2_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_OP2.csv");
    std::ifstream data2(datafolder + "/obstacles_OP2.csv");
    std::ofstream u_file(datafolder + "/u_OP2.csv");
    std::ofstream t_file(datafolder + "/times_OP2.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Road irregularity", 0, 0, 1.5, 0, 0, 1e3);
    obstacles.emplace_back("Carton box", 0, 0, 1.5, 0, 0, 1e6);

    std::string line, line2;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        std::getline(data2, line2);
        parseCsvObstacles(line2, ds_x, ds_y, obstacle_detected);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

void run_op1_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_OP1.csv");
    std::ifstream data2(datafolder + "/obstacles_OP1.csv");
    std::ofstream u_file(datafolder + "/u_OP1.csv");
    std::ofstream t_file(datafolder + "/times_OP1.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Road irregularity", 0, 0, 1.5, 0, 0, 1e3);
    obstacles.emplace_back("Carton box", 0, 0, 1.5, 0, 0, 1e6);

    std::string line, line2;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        std::getline(data2, line2);
        parseCsvObstacles(line2, ds_x, ds_y, obstacle_detected);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}


void run_uo_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_UO.csv");
    std::ofstream u_file(datafolder + "/u_UO.csv");
    std::ofstream t_file(datafolder + "/times_UO.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Car front", 0, 0, 2, 0, 1, 1e3);
    obstacles.emplace_back("Car rear", 0, 0, 2, 0, 1, 1e3);

    std::string line;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

void run_LK_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_LK.csv");
    std::ofstream u_file(datafolder + "/u_LK.csv");
    std::ofstream t_file(datafolder + "/times_LK.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Brick 1", 0, 0, 0.7, 0, 0, 1e3);
    obstacles.emplace_back("Brick 2", 0, 0, 0.7, 0, 0, 1e3);

    std::string line;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;

        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }
        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, false, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

void run_ruts_test(const std::string &datafolder) {
    // data streams
    std::ifstream data(datafolder + "/mpc_input_Ruts.csv");
    std::ifstream lineParams(datafolder + "/curves_ruts.csv");
    std::ofstream u_file(datafolder + "/u_ruts.csv");
    std::ofstream t_file(datafolder + "/time_ruts.csv");

    // parameter definition
    double vehParams[5] = {0.971, 1.566, 1463.0, 1967.8, 0.789};
    double weights[6] = {1e3, 1e10, 1e4, 1e2, 1e2, 5e1};
    double road_padding = 0.2;

    // initial obstacle definition
    std::list<Obstacle_structure> obstacles;
    obstacles.emplace_back("Brick 1", 0, 0, 0.7, 0, 0, 1e3);
    obstacles.emplace_back("Brick 2", 0, 0, 0.7, 0, 0, 1e3);

    std::string line, line2;
    double beta, r, v, df_ref, df_prev;
    double params[2] = {0};
    double roadParams[16] = {0};
    double obstacle_detected[2] = {0};
    double ds_x[2] = {0};
    double ds_y[2] = {0};

    MPC_structure mpc = MPC_structure(10, 2);
    mpc.initMpc(vehParams, weights, road_padding);

    double sol = 0;
    int counterForFirstRun = 0;
    while (std::getline(data, line)) {

        // first run to get everything into cash
        if (counterForFirstRun == 0) {
            parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
            int i = 0;
            std::list<Obstacle_structure>::iterator it;
            for (it = obstacles.begin(); it != obstacles.end(); ++it) {
                it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
                ++i;
            }
            sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, true, obstacles);
            ++counterForFirstRun;
        } else ++counterForFirstRun;
        std::getline(lineParams, line2);
        parseCsvDataFile(line, &beta, &r, &v, &df_ref, &df_prev, params, roadParams, obstacle_detected, ds_x, ds_y);
        parseCsvCurveParams(line2, roadParams);
        int i = 0;
        std::list<Obstacle_structure>::iterator it;
        for (it = obstacles.begin(); it != obstacles.end(); ++it) {
            it->updateObstacle(ds_x[i], ds_y[i], it->radius, obstacle_detected[i], it->undrivable_on, it->weight);
            ++i;
        }

        auto start = high_resolution_clock::now();
        sol = mpc.solveMPC(beta, r, v, df_ref, df_prev, params, roadParams, true, obstacles);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        u_file << sol << std::endl;
        t_file << std::to_string((double) duration.count() / 1000.0) << std::endl;
    }
    data.close();
    u_file.close();
    t_file.close();
}

int main() {
    auto start = high_resolution_clock::now();

    std::string datafolder = "../";

    run_LK_test(datafolder);
    run_ruts_test(datafolder); //in this test there are small deviations from the original windows test
    run_do30_p1_test(datafolder);
    run_do30_p2_test(datafolder);
    run_do70_p1_test(datafolder);
    run_do70_p2_test(datafolder);
    run_uo_test(datafolder);
    run_op2_test(datafolder);
    run_op1_test(datafolder);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    std::cout << std::to_string((double) duration.count() / 1000.0) << std::endl;
}
