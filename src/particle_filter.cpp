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

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
   // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
   //   x, y, theta and their uncertainties from GPS) and all weights to 1.
   // Add random Gaussian noise to each particle.
   // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

   std::default_random_engine gen;
   // Create Gaussian distribution for x, y and theta.
   std::normal_distribution<double> dist_x(x, std[0]);
   std::normal_distribution<double> dist_y(y, std[1]);
   std::normal_distribution<double> dist_psi(theta, std[2]);

   num_particles = 100;
   // initialize particles
    for (size_t i = 0; i < num_particles; i++)
    {
      Particle p = {};
        p.id = i+1;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_psi(gen);
        p.weight = 1.0;
        particles.push_back(p);
    }

    is_initialized = true;
    std::cout << "Is initialized" << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
   // TODO: Add measurements to each particle and add random Gaussian noise.
   // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
   //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   //  http://www.cplusplus.com/reference/random/default_random_engine/


	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(0.0, std_pos[0]);
	std::normal_distribution<double> dist_y(0.0, std_pos[1]);
	std::normal_distribution<double> dist_theta(0.0, std_pos[2]);

	for (int i=0; i<num_particles; i++) {
		Particle* p = &particles[i];
		double xf, yf, thetaf;
		// Switch dynamics based on turn rate.
		if (fabs(yaw_rate) != 0) {
			// non-zero yaw-rate
			xf = p->x + velocity/yaw_rate * (sin(p->theta + yaw_rate*delta_t) - sin(p->theta));
			yf = p->y + velocity/yaw_rate * (cos(p->theta) - cos(p->theta + yaw_rate*delta_t));
			thetaf = p->theta + yaw_rate*delta_t;
		} else {
			// Zero Yaw-Rate Dynamics
			xf = p->x + velocity * delta_t * cos(p->theta);
			yf = p->y + velocity * delta_t * sin(p->theta);
			thetaf = p->theta;
		}
		// Overwrite particle state in particles vector with deterministic state + noise.
		p->x = xf + dist_x(gen);
		p->y = yf + dist_y(gen);
		p->theta = thetaf + dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
   // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
   //   observed measurement to this particular landmark.
   // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
   //   implement this method and use it as a helper during the updateWeights phase.



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
      std::vector<LandmarkObs> observations, Map map_landmarks) {
   // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
   //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
   //   according to the MAP'S coordinate system. You will need to transform between the two systems.
   //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
   //   The following is a good resource for the theory:
   //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   //   and the following is a good resource for the actual equation to implement (look at equation
   //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
   //   for the fact that the map's y-axis actually points downwards.)
   //   http://planning.cs.uiuc.edu/node99.html
	double weights_sum = 0;


    for (size_t i = 0; i < num_particles; i++) {
    	// find close landmarks (in sensor range)
    	Particle *p = &(particles[i]);
		double weight = 1.0;

    	std::vector<Map::single_landmark_s> close_landmarks;
		for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
			if (dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, p->x, p->y) <= sensor_range) {
				close_landmarks.push_back(map_landmarks.landmark_list[j]);
			}
		}



		//transform observations
		std::vector<LandmarkObs>::iterator it;

		for (it = observations.begin(); it < observations.end(); it++) {
			LandmarkObs obs = *it;
			double ox = p->x + cos(-p->theta) * obs.x + sin(-p->theta) * obs.y;
			double oy = p->y - sin(-p->theta) * obs.x + cos(-p->theta) * obs.y;

			double nearest_neighbor; //
			double closest_x; // x coordinate of closest associated landmark
			double closest_y; // y coordinate of closest associated landmark
			//perform nearest neighbor search
			for (int k=0; k<close_landmarks.size(); k++) {
				double temp_neighbor = dist(ox, oy, close_landmarks[k].x_f, close_landmarks[k].y_f);

				if (k == 0 || temp_neighbor < nearest_neighbor) {
					nearest_neighbor = temp_neighbor;
					closest_x = close_landmarks[k].x_f;
					closest_y = close_landmarks[k].y_f;
				}
			}

			double denom = 2 * 3.14 * std_landmark[0] * std_landmark[1];
			double expo_x = pow((ox - closest_x),2)/(2*pow(std_landmark[0],2));
			double expo_y = pow((oy - closest_y),2)/(2*pow(std_landmark[1],2));
			double guassian_prob = (1/denom) * exp(-(expo_x + expo_y));

			weight = weight * guassian_prob;
		}

		p->weight = weight;
		weights_sum = weights_sum + weight;
    }

	for (int i=0; i < num_particles; i++) {
		particles[i].weight = particles[i].weight / weights_sum;
		weights.push_back(particles[i].weight);
	}

}



void ParticleFilter::resample() {
   // TODO: Resample particles with replacement with probability proportional to their weight.
   // NOTE: You may find std::discrete_distribution helpful here.
   //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::default_random_engine gen;
	std::discrete_distribution<> distribution(weights.begin(), weights.end());
	std::vector<Particle> re_particles;

	for (int i = 0; i < num_particles; i++) {
		int weighted_index = distribution(gen);
		re_particles.push_back(particles[weighted_index]);

		//std::cout << weighted_index << std::endl;
		//std::cout << "Paritcle weight = " << particles[weighted_index].weight << std::endl;
		//particles[i] = particles[weighted_index];
	}
	particles = re_particles;
}

void ParticleFilter::write(std::string filename) {
   // You don't need to modify this file.
   std::ofstream dataFile;
   dataFile.open(filename, std::ios::app);
   for (int i = 0; i < num_particles; ++i) {
      dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
   }
   dataFile.close();
}
