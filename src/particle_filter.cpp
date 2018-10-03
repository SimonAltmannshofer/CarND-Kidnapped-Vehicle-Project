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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	
	
	// create normal (Gaussian) distributions
	if (is_initialized) {
		return;
	}
	else {
		default_random_engine gen;
		normal_distribution<double> dist_x(0.0, std[0]);
		normal_distribution<double> dist_y(0.0, std[1]);
		normal_distribution<double> dist_theta(0.0, std[2]);
		double dist_on = 0.0;
		
		num_particles = 150;
		for(int i=0; i<num_particles; i++){
			weights.push_back(1.0);
			
			Particle new_particle;
			new_particle.id = i;
			new_particle.x = x + dist_on*dist_x(gen);
			new_particle.y = y + dist_on*dist_y(gen);
			new_particle.theta = theta + dist_on*dist_theta(gen);
			new_particle.weight = 1.0;
			
			particles.push_back(new_particle);
		}
		
		is_initialized = true;
	}
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	normal_distribution<double> dist_x(0.0, std_pos[0]);
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);
	double dist_on = 0.0;
	
	if(abs(yaw_rate) < 0.001){
		for(int i=0; i<num_particles; i++){
			particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta) + dist_on*dist_x(gen);
			particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta) + dist_on*dist_y(gen);
			particles[i].theta = particles[i].theta + dist_on*dist_theta(gen);
		}
	}
	else{
		for(int i=0; i<num_particles; i++){
			particles[i].x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + dist_x(gen);
			particles[i].y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) + dist_y(gen);
			particles[i].theta = particles[i].theta + yaw_rate*delta_t + dist_theta(gen);
		}
	}
	
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for (int i = 0; i < num_particles; i++) {
		double x_p = particles[i].x;
		double y_p = particles[i].y;
		double theta_p = particles[i].theta;

		
		vector<LandmarkObs> landmarks_in_range;
		vector<LandmarkObs> transformed_landmarks;
		vector<LandmarkObs> associated_landmarks;

		// transform observations
		for (unsigned int j = 0; j < observations.size(); j++) {
			int id_c = observations[j].id;
			double x_c = observations[j].x;
			double y_c = observations[j].y;

			double x_m = x_p + x_c * cos(theta_p) - y_c * sin(theta_p);
			double y_m = y_p + y_c * cos(theta_p) + x_c * sin(theta_p);

			LandmarkObs transformed_landmark;
			transformed_landmark.id = id_c;
			transformed_landmark.x = x_m;
			transformed_landmark.y = y_m;

			transformed_landmarks.push_back(transformed_landmark);
		}
		
		// get all landmarks in sensor range
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			int id_landmark = map_landmarks.landmark_list[j].id_i;
			double x_landmark = map_landmarks.landmark_list[j].x_f;
			double y_landmark = map_landmarks.landmark_list[j].y_f;

			double x_dist = x_landmark - x_p;
			double y_dist = y_landmark - y_p;
			double dist = sqrt(x_dist * x_dist + y_dist * y_dist);

			if (dist < sensor_range) {
				LandmarkObs landmark_in_range;
				landmark_in_range.id = id_landmark;
				landmark_in_range.x = x_landmark;
				landmark_in_range.y = y_landmark;

				landmarks_in_range.push_back(landmark_in_range);
			}
		}

		// get associations
		for (unsigned int j = 0; j < transformed_landmarks.size(); j++) {
			int id_associated = -1;
			double x_associated = 0.0;
			double y_associated = 0.0;
			double dist_min = 1000000;

			for (unsigned int k = 0; k < landmarks_in_range.size(); k++) {
				double x_delta = landmarks_in_range[k].x - transformed_landmarks[j].x;
				double y_delta = landmarks_in_range[k].y - transformed_landmarks[j].y;
				double dist = x_delta * x_delta + y_delta * y_delta;

				if (dist < dist_min) {
					id_associated = k;
					x_associated = transformed_landmarks[j].x;
					y_associated = transformed_landmarks[j].y;
					dist_min = dist;
				}
			}
			transformed_landmarks[j].id = id_associated;
			LandmarkObs associated_landmark;
			
			associated_landmark.id = id_associated;
			associated_landmark.x = x_associated;
			associated_landmark.y = y_associated;
			
			associated_landmarks.push_back(associated_landmark);
		}
		
		// calculate weights
		double std_x = std_landmark[0];
		double std_y = std_landmark[1];
		double denom = sqrt(2.0 * M_PI * std_x * std_y);
		double denom_1 = 2 * std_x * std_x;
		double denom_2 = 2 * std_y * std_y;

		double weight = 1.0;
		for (unsigned int j = 0; j < transformed_landmarks.size(); j++) {
			int id_m = transformed_landmarks[j].id;
			double x_m = transformed_landmarks[j].x;
			double y_m = transformed_landmarks[j].y;

			double x_predict = landmarks_in_range[id_m].x;
			double y_predict = landmarks_in_range[id_m].y;

			double x_delta = x_m - x_predict;
			double y_delta = y_m - y_predict;

			double part_1 = x_delta * x_delta / denom_1;
			double part_2 = y_delta * y_delta / denom_2;

			weight *= exp(-(part_1 + part_2)) / denom;
		}
		if (weight == 0) {
			particles[i].weight = 0.0000001;
			weights[i] = 0.0000001;
		} else {
			particles[i].weight = weight;
			weights[i] = weight;
		}
	}
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine generator;
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	vector<Particle> new_particles;
	
	
	for(int i=0; i<num_particles; i++){
		const int index = distribution(generator);
		//cout << "index: " << index << endl;
		
		Particle new_particle;
		new_particle.id = index;
		new_particle.x = particles[index].x;
		new_particle.y = particles[index].y;
		new_particle.theta = particles[index].theta;
		new_particle.weight = 1.0;
		new_particles.push_back(new_particle);
	}
	
	particles = new_particles;
	
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
	
	return particle;
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
