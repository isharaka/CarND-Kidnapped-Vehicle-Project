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

void ParticleFilter::init(double x, double y, double theta, double std[]){
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	std::cout << "ParticleFilter::init" << std::endl;

	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i=0; i<num_particles; i++) {
		weights.push_back(1.0);

		Particle particle;

		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = weights[i];

		particles.push_back(particle);
	}

	for(int i = 0; i< weights.size(); i++) {
		std::cout << weights[i] << std::endl;
	}

	for(int i = 0; i< particles.size(); i++) {
		std::cout << particles[i].x << " "
				<< particles[i].y << " "
				<< particles[i].theta << " "<< std::endl;
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	for(int i = 0; i< particles.size(); i++) {
		if (std::abs(yaw_rate) > 0.001) {
			particles[i].x += velocity * (std::sin(particles[i].theta + yaw_rate*delta_t) - std::sin(particles[i].theta)) / yaw_rate;
			particles[i].y += velocity * (-std::cos(particles[i].theta - yaw_rate*delta_t) + std::cos(particles[i].theta)) / yaw_rate;
		} else {
			particles[i].x += velocity * std::sin(particles[i].theta) * delta_t;
			particles[i].y += velocity * std::cos(particles[i].theta) * delta_t;
		}
		particles[i].theta += yaw_rate*delta_t;
		//angle normalization
	    while (particles[i].theta> M_PI) particles[i].theta-=2.*M_PI;
	    while (particles[i].theta<-M_PI) particles[i].theta+=2.*M_PI;
	}

		std::cout << " delta_t " << delta_t << " velocity " << velocity << " yaw rate " << yaw_rate << std::endl;

	for(int i = 0; i< particles.size(); i++) {
		std::cout << "i " << i << " x " << particles[i].x << " y " << particles[i].y << " theta " << particles[i].theta << std::endl;
	}

}

void ParticleFilter::dataAssociation(/*std::vector<LandmarkObs> predicted*/const Map &map_landmarks, const std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int p = 0; p < num_particles; ++p) {
		for (int o = 0; o < observations.size(); ++o) {
			// transform each sensor measurement ,which is in local(car) coordinates to global(map) coordinates
			LandmarkObs transformed_observation;

			transformed_observation.x = particles[p].x 
									  	+ observations[o].x * std::cos(particles[p].theta)
										- observations[o].y * std::sin(particles[p].theta);
			transformed_observation.y = particles[p].y 
										+ observations[o].x * std::sin(particles[p].theta)
										+ observations[o].y * std::cos(particles[p].theta);

			//std::cout << " o " << o << " obs " << observations[o].x << ":" << observations[o].y <<
			//			" trans " << transformed_observation.x << ":" << transformed_observation.y << std::endl;

			double min_distance = std::numeric_limits<double>::max();
			int closest_landmark;

			for (int l = 0; l < map_landmarks.landmark_list.size(); ++l) {
				// Calculate euclidean distance between each landmark pos and sensor pos
				double distance = pow(map_landmarks.landmark_list[l].x_f - transformed_observation.x, 2)
								+ pow(map_landmarks.landmark_list[l].y_f - transformed_observation.y, 2);

				// if this distance is less than min distance,associate current particle with that landmark index
				if (distance < min_distance) 
					closest_landmark = map_landmarks.landmark_list[l].id_i;
			}

			particles[p].associations.push_back(closest_landmark);
			particles[p].sense_x.push_back(transformed_observation.x);
			particles[p].sense_y.push_back(transformed_observation.y);
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	dataAssociation(map_landmarks, observations);

	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];

	weights.clear();

	for (int p = 0; p < num_particles; ++p) {
		double weight = 1.0;
		for (int o = 0; o < observations.size(); ++o) {
			double mu_x = map_landmarks.landmark_list[particles[p].associations[o]].x_f;
			double mu_y = map_landmarks.landmark_list[particles[p].associations[o]].y_f;

			double x_obs = particles[p].sense_x[o];
			double y_obs = particles[p].sense_y[o];

			double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));
			double exponent = (pow(x_obs - mu_x,2))/(2 * pow(sig_x,2)) + (pow(y_obs - mu_y,2))/(2 * pow(sig_y,2));

			weight *= gauss_norm * exp(-exponent);
		}

		particles[p].weight = weight;
		weights.push_back(weight);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(weights.begin(), weights.end());

    std::vector<Particle> particles_new;

    for(int n=0; n<num_particles; ++n) {
        particles_new.push_back(particles[d(gen)]);
    }

    particles.clear();
    particles = particles_new;
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
