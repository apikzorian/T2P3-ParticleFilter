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
	num_particles = 100;

	std::default_random_engine gen;
	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);

	// Initialize particles.
	for (int i=0; i<num_particles; i++) {
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
	}
	is_initialized = true;
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

	// Add measurements and random Gaussian noise for each particle
	for (int i=0; i<num_particles; i++) {
		Particle* p = &particles[i];
		double xf, yf, thetaf;
		// Non-zero yaw
		if (fabs(yaw_rate) != 0) {
			// non-zero yaw-rate
			xf = p->x + velocity/yaw_rate * (sin(p->theta + yaw_rate*delta_t) - sin(p->theta));
			yf = p->y + velocity/yaw_rate * (cos(p->theta) - cos(p->theta + yaw_rate*delta_t));
			thetaf = p->theta + yaw_rate*delta_t;
		} else {
			// Zero yaw
			xf = p->x + velocity * delta_t * cos(p->theta);
			yf = p->y + velocity * delta_t * sin(p->theta);
			thetaf = p->theta;
		}

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

	// Sum of weights (for normalization)
	double weights_sum = 0;

	for (int i=0; i<num_particles; i++) {
		double weight = 1.0;

		// Fill vector with landmarks that are within sensor range
		std::vector<Map::single_landmark_s> sensed_landmarks;
		for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
			if (dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particles[i].x, particles[i].y) <= sensor_range) {
				sensed_landmarks.push_back(map_landmarks.landmark_list[j]);
			}
		}

		// Convert observations from vehicle to map coordinates
		for (int j=0; j<observations.size(); j++) {

			// nearest neighbor x and y values
			double nearest_neighbor, closest_x, closest_y;
			double ox = particles[i].x + cos(-particles[i].theta) * observations[j].x + sin(-particles[i].theta) * observations[j].y;
			double oy = particles[i].y - sin(-particles[i].theta) * observations[j].x + cos(-particles[i].theta) * observations[j].y;

			for (int k=0; k<sensed_landmarks.size(); k++) {
				double temp_neighbor = dist(ox, oy, sensed_landmarks[k].x_f, sensed_landmarks[k].y_f);

				// Update nearest_neighbor if current value is smaller than saved
				if (k==0 || (temp_neighbor < nearest_neighbor)) {
					nearest_neighbor = temp_neighbor;
					closest_x = sensed_landmarks[k].x_f;
					closest_y = sensed_landmarks[k].y_f;
				}
			} // end of nearest neighbor search loop

			// Assign weight of particle based on Multivariate Gaussian Probability
			double denom = 2 * 3.14 * std_landmark[0] * std_landmark[1];
			double expo_x = pow((ox - closest_x),2)/(2*pow(std_landmark[0],2));
			double expo_y = pow((oy - closest_y),2)/(2*pow(std_landmark[1],2));
			double guassian_prob = (1/denom) * exp(-(expo_x + expo_y));
			//Update weight
			weight *= guassian_prob;
		}

		particles[i].weight = weight;
		weights_sum += particles[i].weight;
	} // end of particle loop

	// Normalize weights
	for (int i=0; i < num_particles; i++) {
		particles[i].weight /= weights_sum;
		weights.push_back(particles[i].weight);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution


	//std::random_device rd;
	//std::mt19937 gen(rd());

	std::default_random_engine gen;
	std::discrete_distribution<> distribution(weights.begin(), weights.end());
	std::vector<Particle> re_particles;

	for (int i = 0; i < num_particles; i++) {
		int weighted_index = distribution(gen);
		re_particles.push_back(particles[weighted_index]);

		std::cout << "Particle weight_index = " << weighted_index << std::endl;
		std::cout << "Particle weight = " << particles[weighted_index].weight << std::endl;

	}
	std::cout << "Particles size = " << particles.size() << std::endl;
	std::cout << "Re_particles size = " << re_particles.size() << std::endl;

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
