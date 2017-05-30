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
	num_particles = 1000;

	std::default_random_engine generator;
	std::normal_distribution<double> distribution_x(x,std[0]);
	std::normal_distribution<double> distribution_y(y,std[1]);
	std::normal_distribution<double> distribution_theta(theta,std[2]);

	for (size_t i=0;i<num_particles;i++)
	{
		Particle particle;
		particle.id = int(i);
		particle.x = distribution_x(generator);
		particle.y = distribution_y(generator);
		particle.theta = distribution_theta(generator);
		particle.weight = 1;
		particles.push_back(particle);
		weights.push_back(1);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::default_random_engine generator;
	std::normal_distribution<double> noise_x(0,std_pos[0]);
	std::normal_distribution<double> noise_y(0,std_pos[1]);
	std::normal_distribution<double> noise_theta(0,std_pos[2]);

	for (size_t i=0;i<num_particles;i++)
	{
		double x0 = particles[i].x;
		double y0 = particles[i].y;
		double theta0 = particles[i].theta;

		if (fabs(yaw_rate) < 0.001)
		{
			particles[i].x = x0+velocity*delta_t*cos(theta0)+noise_x(generator);
			particles[i].y = y0+velocity*delta_t*sin(theta0)+noise_y(generator);
			particles[i].theta = theta0+noise_theta(generator);
		}
		else
		{
			particles[i].x = x0+velocity/yaw_rate*(sin(theta0+yaw_rate*delta_t)-sin(theta0))+noise_x(generator);
			particles[i].y = y0+velocity/yaw_rate*(cos(theta0)-cos(theta0+yaw_rate*delta_t))+noise_y(generator);
			particles[i].theta = theta0+yaw_rate*delta_t+noise_theta(generator);
		}
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (size_t o=0;o<observations.size();o++)
	{
		observations[o].id = -1;
		LandmarkObs obs = observations[o];
		double min_dist = std::numeric_limits<double>::max();

		for (size_t p=0;p<predicted.size();p++)
		{
			LandmarkObs pred = predicted[p];
			double x_diff = obs.x-pred.x;
			double y_diff = obs.y-pred.y;
			double dist = x_diff*x_diff+y_diff*y_diff;

			if (dist<min_dist)
			{
				min_dist = dist;
				observations[o].id = p;
			}
		}
	}


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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	double normalizer = 1.0/(2.0*M_PI*std_landmark[0]*std_landmark[1]);

	double weights_normalizer = 0;
	for (size_t i=0;i<num_particles;i++)
	{
		Particle particle = particles[i];

		// Predict where landmarks are in VEHICLE coordinate system
		std::vector<LandmarkObs> landmarks_predicted;
		for (size_t p=0;p<map_landmarks.landmark_list.size();p++)
		{
			double x_diff = map_landmarks.landmark_list[p].x_f-particle.x;
			double y_diff = map_landmarks.landmark_list[p].y_f-particle.y;


			LandmarkObs landmark_predicted;
			landmark_predicted.id = map_landmarks.landmark_list[p].id_i;
			landmark_predicted.x = cos(particle.theta)*x_diff+sin(particle.theta)*y_diff;
			landmark_predicted.y = -sin(particle.theta)*x_diff+cos(particle.theta)*y_diff;

			double landmark_dist_squared = landmark_predicted.x*landmark_predicted.x+landmark_predicted.y*landmark_predicted.y;
			if (landmark_dist_squared <= sensor_range*sensor_range)
				landmarks_predicted.push_back(landmark_predicted);
		}

		// Associate nearest landmarks
		dataAssociation(landmarks_predicted, observations);

		// Transform observations into MAP coordinate system (used by simulator)
		std::vector<LandmarkObs> observations_transformed;
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;
		for (size_t o=0;o<observations.size();o++)
		{
			LandmarkObs obs = observations[o];

			if (obs.id >= 0) // If it is associated with a landmark
			{
				// Transform measurement from VEHICLE coordinate system to MAP coordinate system
				double transformed_x = particle.x+obs.x*cos(particle.theta)-obs.y*sin(particle.theta);
				double transformed_y = particle.y+obs.x*sin(particle.theta)+obs.y*cos(particle.theta);

				associations.push_back(landmarks_predicted[obs.id].id);
				sense_x.push_back(transformed_x);
				sense_y.push_back(transformed_y);
			}
		}
		particles[i] = SetAssociations(particle, associations, sense_x, sense_y);

		// Calculate particle weights
		particles[i].weight = 1;
		bool has_observations = false;
		for (size_t o=0;o<observations.size();o++)
		{
			if (observations[o].id >= 0)
			{
				double x = observations[o].x;
				double ux = landmarks_predicted[observations[o].id].x;
				double y = observations[o].y;
				double uy = landmarks_predicted[observations[o].id].y;
				double x_diff = x-ux;
				double y_diff = y-uy;
				particles[i].weight *= normalizer*exp(-(x_diff*x_diff/(2.0*std_landmark[0]*std_landmark[0])+y_diff*y_diff/(2.0*std_landmark[1]*std_landmark[1])));
				has_observations = true;
			}
		}
		if (!has_observations)
			particles[i].weight = 0;

		weights_normalizer += particles[i].weight; // Calculate weight normalizer
	}

	// Normalize particle weights
	for (size_t i=0;i<num_particles;i++)
	{
		weights[i] = particles[i].weight/weights_normalizer;
		particles[i].weight = weights[i];
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<int> d(weights.begin(),weights.end());

	std::vector<Particle> particles_resampled;

	for (size_t i=0;i<num_particles;i++)
	{
		particles_resampled.push_back(particles[d(gen)]);
	}

	particles = particles_resampled;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
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
