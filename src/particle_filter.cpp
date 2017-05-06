#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

/**
 * void ParticleFilter::init(double x, double y, double theta, double std[])
 * Initializes particle filter by initializing particles to Gaussian
 * distribution around first position and all the weights to 1.
 * @param x Initial x position [m] (simulated estimate from GPS)
 * @param y Initial y position [m]
 * @param theta Initial orientation [rad]
 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {

   //Create Gaussian distribution for x, y and theta
   std::default_random_engine gen;
   std::normal_distribution<double> dist_x(x, std[0]);
   std::normal_distribution<double> dist_y(y, std[1]);
   std::normal_distribution<double> dist_psi(theta, std[2]);

   num_particles = 100;
   // Initialize particles
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
}

/*  ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
	Predicts the state for the next time step using the process model.
    * @param delta_t Time between time step t and t+1 in measurements [s]
    * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
*/
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(0.0, std_pos[0]);
	std::normal_distribution<double> dist_y(0.0, std_pos[1]);
	std::normal_distribution<double> dist_theta(0.0, std_pos[2]);

	// Add measurements and random Gaussian noise for each particle
	for (int i=0; i<num_particles; i++) {
		Particle* p = &particles[i];
		double xf, yf, thetaf;
		// Non-zero yaw-rate
		if (fabs(yaw_rate) != 0) {
			xf = p->x + velocity/yaw_rate * (sin(p->theta + yaw_rate*delta_t) - sin(p->theta));
			yf = p->y + velocity/yaw_rate * (cos(p->theta) - cos(p->theta + yaw_rate*delta_t));
			thetaf = p->theta + yaw_rate*delta_t;
		} else {
			// Zero yaw-rate
			xf = p->x + velocity * delta_t * cos(p->theta);
			yf = p->y + velocity * delta_t * sin(p->theta);
			thetaf = p->theta;
		}

		// Set new particle measurements with noise
		p->x = xf + dist_x(gen);
		p->y = yf + dist_y(gen);
		p->theta = thetaf + dist_theta(gen);
	}

}

/**
 * ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	std::vector<LandmarkObs> observations, Map map_landmarks)
 * Updates the weights for each particle based on the likelihood of the observed measurements.
 * @param sensor_range Range [m] of sensor
 * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
 *   standard deviation of bearing [rad]]
 * @param observations Vector of landmark observations
 * @param map Map class containing map landmarks
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	std::vector<LandmarkObs> observations, Map map_landmarks) {

	// Sum of weights (for normalization)
	double weights_sum = 0;

	// Iterate through particles, update weights
	for (int i=0; i<num_particles; i++) {
		double weight = 1.0;
    	Particle *p = &(particles[i]);

		// Fill vector with landmarks that are within sensor range
		std::vector<Map::single_landmark_s> sensed_landmarks;
		for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
			if (dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, p->x, p->y) <= sensor_range) {
				sensed_landmarks.push_back(map_landmarks.landmark_list[j]);
			}
		}

		// Convert observations from vehicle to map coordinates
		for (int j=0; j<observations.size(); j++) {

			// nearest neighbor x and y values
			double nearest_neighbor, closest_x, closest_y;
			double ox = p->x + cos(-p->theta) * observations[j].x + sin(-p->theta) * observations[j].y;
			double oy = p->y - sin(-p->theta) * observations[j].x + cos(-p->theta) * observations[j].y;

			for (int k=0; k<sensed_landmarks.size(); k++) {
				double temp_neighbor = dist(ox, oy, sensed_landmarks[k].x_f, sensed_landmarks[k].y_f);

				// Update nearest_neighbor if current value is smaller than saved
				if (k==0 || (temp_neighbor < nearest_neighbor)) {
					nearest_neighbor = temp_neighbor;
					closest_x = sensed_landmarks[k].x_f;
					closest_y = sensed_landmarks[k].y_f;
				}
			}

			// Assign weight of particle based on Multivariate Gaussian Probability
			double denom = 2 * 3.14 * std_landmark[0] * std_landmark[1];
			double expo_x = pow((ox - closest_x),2)/(2*pow(std_landmark[0],2));
			double expo_y = pow((oy - closest_y),2)/(2*pow(std_landmark[1],2));
			double guassian_prob = (1/denom) * exp(-(expo_x + expo_y));
			weight *= guassian_prob;
		}

		//Update weight of particle
		p->weight = weight;
		weights_sum += p->weight;
	}

	// Normalize weights and update weights vector with normalized weight
	for (int i=0; i < num_particles; i++) {
		particles[i].weight /= weights_sum;
		weights.push_back(particles[i].weight);
	}
}

/**
 * ParticleFilter::resample()
 * Resamples from the updated set of particles to form the new set of particles.
 */
void ParticleFilter::resample() {

	//Use discrete distribution on weights vector
	std::default_random_engine gen;
	std::discrete_distribution<> distribution(weights.begin(), weights.end());
	std::vector<Particle> re_particles; // vector of resampled particles

	// Fill re_particles with particles chosen by discrete distribution function
	for (int i = 0; i < num_particles; i++) {
		int weighted_index = distribution(gen);
		re_particles.push_back(particles[weighted_index]);
	}
	// reset particles to resampled particles and clear weights vector
	particles = re_particles;
	weights.clear();
}

/**
 * ParticleFilter::write(std::string filename)
 * Writes particle positions to a file.
 * @param filename File to write particle positions to.
 */
void ParticleFilter::write(std::string filename) {
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
