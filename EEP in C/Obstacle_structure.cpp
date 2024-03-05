//
// Created by efrem on 11/1/2023.
//

#include <sstream>

class Obstacle_structure {
public:
    std::string obstacle_name; // tag name for the obstacle to keep and order
    double dsx; // x (along the vehicle) position wrt vehicle CG, positive for forward
    double dsy; // y (across the vehicle) position wrt vehicle CG, positive for left (?)
    double radius; // size of the obstacle
    double obs_on; // (0/1 - OFF/ON) is it observable? If no, it is not considered in the MPC
    double undrivable_on; // (0/1 - OFF/ON) is it undrivable? If yes, it creates additional constraints for the car avoidance
    double weight; // weighting factor, can be used for vehicle prioritization

    Obstacle_structure(std::string obstacle_name, double dsx, double dsy, double radius, double obs_on,
                       double undrivable_on, double weight) {
        this->obstacle_name = obstacle_name;
        this->dsx = dsx;
        this->dsy = dsy;
        this->radius = radius;
        this->obs_on = obs_on;
        this->undrivable_on = undrivable_on;
        this->weight = weight;
    }

    void updateObstacle(double new_dsx, double new_dsy, double new_radius, double new_obs_on, double new_undrivable_on,
                        double new_weight) {
        this->dsx = new_dsx;
        this->dsy = new_dsy;
        this->radius = new_radius;
        this->obs_on = new_obs_on;
        this->undrivable_on = new_undrivable_on;
        this->weight = new_weight;
    }
};
