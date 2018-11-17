/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;


const int PARTICLE_NUM = 100;
static default_random_engine generator;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
    //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    // GPS measurement uncertainty std[] [x [m], y [m], theta [rad]]

    // init particles
    // init particle's x,y position, theta and weight
    normal_distribution<double> distribution_x(x, std[0]); //set x generator at mean given by sensed x with stdev in x
    normal_distribution<double> distribution_y(y, std[1]); //set y generator at mean given by sensed y with stdev in y
    normal_distribution<double> distribution_theta(theta, std[2]); //set theta generator at mean given by sensed theta with stdev in theta

    num_particles = PARTICLE_NUM;

    for (int i = 0; i < num_particles; i++) {
        Particle newparticle; //make a new particle and add the Gaussian noise
        newparticle.id = i;
        newparticle.x = distribution_x(generator);
        newparticle.y = distribution_y(generator);
        newparticle.theta = distribution_theta(generator);
        newparticle.weight = 1.0;

        particles.push_back(newparticle); //put new particle in the array of particles
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // params 
    // delta_t: 时间间隔
    // std_pos: 标准差
    // velocity: 速度
    // yaw_rate: 角速度

    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    double mean_x;
    double mean_y;
    double mean_theta;

    for (int i = 0; i < num_particles; i++) {
        
        if (fabs(yaw_rate) < 0.00001) {
            // dx = v * t * cos(theta)
            mean_x = particles[i].x + delta_t * velocity * cos(particles[i].theta);
            // dy = v * t * sin(theta)
            mean_y = particles[i].y + delta_t * velocity * sin(particles[i].theta);
            mean_theta = particles[i].theta + yaw_rate;
        }
        else {
            //The equations for updating x, y and the yaw angle when the yaw rate is not equal to zero (Udacity SDCND Term2 Lesson 14. slide 8)
            mean_x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
            mean_y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
            //std::cout << "p_theta: " << particles[i].theta << std::endl;
            mean_theta = particles[i].theta + yaw_rate * delta_t;
            //std::cout << "mean_theta: " << mean_theta << std::endl;
        }

        // add random noise
        normal_distribution<double> distribution_x(mean_x, std_pos[0]); //set x generator at mean given by sensed x with stdev in x
        normal_distribution<double> distribution_y(mean_y, std_pos[1]); //set y generator at mean given by sensed y with stdev in y
        normal_distribution<double> distribution_theta(mean_theta, std_pos[2]); //set theta generator at mean given by sensed theta with stdev in theta

        //Update particle x,y,theta using Gaussian generator:
        particles[i].x = distribution_x(generator);
        particles[i].y = distribution_y(generator);
        particles[i].theta = distribution_theta(generator);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
    //   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation 
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html

    // purpose: update the weight for each particle

    std::vector<LandmarkObs> predicted;
    for (int i = 0; i < num_particles; i++) {
        predicted.clear();

        // get the particle current position and theta
        double x1 = particles[i].x;
        double y1 = particles[i].y;
        double theta = particles[i].theta;
        

        // observations is vector, each landmark has an observations.
        for (int j = 0; j < observations.size();j++) {
            // 坐标系变换
            double x2 = x1 + observations[j].x * cos(theta) + observations[j].y * sin(-theta); //convert observed (Local) coord to Global coord
            double y2 = y1 + observations[j].x * sin(theta) + observations[j].y * cos(theta); //convert observed (Local) coord to Global coord


            // 变换前后的距离超过sensor_range（传感器范围），说明粒子位置与车辆位置差距过大，直接跳过
            if (dist(x1, y1, x2, y2) <= sensor_range) {
                //create a prediction add to predicted list
                LandmarkObs pred;
                // Find nearest neighbour from map - using naive search, quite costly O(n)?, future considerations include quadtree search
                double lastmin = 1000.0;

                // 从地图中找出看距离变换后最近的坐标点，并与地图标定点计算距离差值
                //std::cout << "landmarklistsize: " << map_landmarks.landmark_list.size() << endl;
                for (int k = 0; k < map_landmarks.landmark_list.size();k++) {
                    double min = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, x2, y2); // The distance to nearest map landmark  
                    if (min < lastmin) {
                        lastmin = min;
                        pred.x = x2 - map_landmarks.landmark_list[k].x_f; //doing the subtraction here for efficiency later in the Gaussian
                        pred.y = y2 - map_landmarks.landmark_list[k].y_f; //doing the subtraction here for efficiency later in the Gaussian
                        pred.id = map_landmarks.landmark_list[k].id_i;
                    }
                }
                predicted.push_back(pred);
            }
        }
        //Multivariate Gaussian, simplify calc by collapsing variables:
        double sigxx = std::pow(std_landmark[0], 2);
        double sigyy = std::pow(std_landmark[1], 2);
        double coeff = 1.0 / (2.0 * M_PI*std_landmark[0] * std_landmark[1]); //1.7683

        double weight = 1.0; // initialize weight
        for (int l = 0; l < predicted.size(); l++) {
            double xdiffsq = std::pow(predicted[l].x, 2);
            double ydiffsq = std::pow(predicted[l].y, 2);
            weight *= coeff * exp(-0.5 * (xdiffsq / sigxx + ydiffsq / sigyy));
        }
        // Done calculating combined weight of all observations, assign to particle:
        particles[i].weight = weight;
    }


    // our resample method need not normalize weights
    // Normalize Weights:
    // double sum = 0;
    // for (int i = 0; i < num_particles;i++) {
    //     sum += particles[i].weight;
    // }
    // for (int i = 0; i < num_particles;i++) {
    //     // normalize
    //     particles[i].weight *= 1/sum;
    // }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight. 
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    // assemble a list of all weights

    // create temporary particle list
    static std::uniform_int_distribution<unsigned> uniform_int(0, num_particles);
    static std::uniform_real_distribution<> uniform_float(0.0, 1.0);
    std::vector<Particle> p_temp;
    
    // find the first index
    int index = uniform_int(generator);
    // std::cout << "first_index=" << index << std::endl;
    float beta = 0;
    // find max of w
    float mw = 0.0;
    for (int i = 0; i < num_particles; i++) {
        if (particles[i].weight > mw) {
            mw = particles[i].weight;
        }
    }
    // std::cout << "mw=" << mw << std::endl;

    for (int i = 0; i < num_particles; i++) {
        beta += 2 * mw * uniform_float(generator);
        // std::cout << "beta=" << beta << std::endl;
        while (beta > particles[index].weight) {

            // std::cout << "weight=" << particles[index].weight << std::endl;
            beta -= particles[index].weight;
            index = (index + 1) % num_particles;
        }
        p_temp.push_back(particles[index]);
        // std::cout << "p_temp index=" << index << std::endl;
    }

    for (int i = 0; i < num_particles; i++) {
        particles[i] = p_temp[i];
    }

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
