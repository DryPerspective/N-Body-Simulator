#include "Planet.h"

//Constructors
Planet::Planet(const std::string& inName, const double inMass, const PhysicsVector<3>& inPos, const PhysicsVector<3>& inVel) :
	m_name{ inName }, m_mass{ inMass }, m_position{ inPos }, m_velocity{ inVel }{}

Planet::Planet(const std::string& inName, const double inMass, const PhysicsVector<3>& inPos, const PhysicsVector<3>& inVel, const PhysicsVector<3>& inAcc) :
	m_name{ inName }, m_mass{ inMass }, m_position{ inPos }, m_velocity{ inVel }, m_acceleration{ inAcc }{}

Planet::Planet(const double inMass, const PhysicsVector<3>& inPos, const PhysicsVector<3>& inVel, const PhysicsVector<3>& inAcc) :
	m_mass{ inMass }, m_position{ inPos }, m_velocity{ inVel }, m_acceleration{ inAcc }{}

Planet::Planet(const std::string& inName) : m_name{ inName } {}


//Getters and setters
double Planet::getMass() const {
	return m_mass;
}
const PhysicsVector<3>& Planet::getPosition() const {
	return m_position;
}
const PhysicsVector<3>& Planet::getVelocity() const {
	return m_velocity;
}
const PhysicsVector<3>& Planet::getAcceleration() const {
	return m_acceleration;
}
const std::string& Planet::getName() const {
	return m_name;
}

void Planet::setMass(double inMass) {
	m_mass = inMass;
}
void Planet::setPosition(const PhysicsVector<3>& inPos) {
	m_position = inPos;
}
void Planet::setVelocity(const PhysicsVector<3>& inVel) {
	m_velocity = inVel;
}
void Planet::setAcceleration(const PhysicsVector<3>& inAcc) {
	m_acceleration = inAcc;
}
void Planet::setName(const std::string& inName) {
	m_name = inName;
}


//Planetary functions

//Calculate the acceleration on this planet caused by a single other planet. This is given by -(G M)/r^2 times the r unit vector, where r is the distance between the two objects.
PhysicsVector<3> Planet::calcAcceleration(const Planet& inPlanet) {
	double r{ (m_position - inPlanet.m_position).magnitude() };
	PhysicsVector<3> outVector { (m_position - inPlanet.m_position).getUnitVector() };
	outVector.scaleVector(-(G * inPlanet.m_mass) / pow(r, 2));
	return outVector;
}
//Calculate the total acceleration felt on a planet, as caused by every other planet in the system, then set acceleration accordingly.
void Planet::updateAccelerationEuler(const planetArray_t& inPlanets) {
	PhysicsVector<3> sumAcceleration{ 0,0,0 };
	for (const auto& planet : inPlanets) {
		if (this == &planet) continue;											//Skip over this planet if it appears in the array
		sumAcceleration += calcAcceleration(planet);
	}
	m_acceleration = sumAcceleration;
}
//Use the first order Euler approximation to update the position of a planet. Namely Position2 = Position1 + velocity * time step.
void Planet::updatePositionEuler(const double timeStep) {
	double newX{ m_position.getX() + (m_velocity.getX() * timeStep) };
	double newY{ m_position.getY() + (m_velocity.getY() * timeStep) };
	double newZ{ m_position.getZ() + (m_velocity.getZ() * timeStep) };
	m_position = { newX,newY,newZ };
}
//Use the first order Euler approximation to update the velocity of a planet. Namely Velocity2 = Velocity1 + acceleration * time step.
void Planet::updateVelocityEuler(const double timeStep) {
	double newX{ m_velocity.getX() + (m_acceleration.getX() * timeStep) };
	double newY{ m_velocity.getY() + (m_acceleration.getY() * timeStep) };
	double newZ{ m_velocity.getZ() + (m_acceleration.getZ() * timeStep) };
	m_velocity = { newX,newY,newZ };
}


/*  
* Algorithm implementation functions.
* These combine the individual functions above into a single, composite method.
*/
//First the Euler approximation. This one is fairly simple to implement.
void Planet::updateEuler(const planetArray_t& inPlanets, double timeStep) {
	updateAccelerationEuler(inPlanets);
	updatePositionEuler(timeStep);
	updateVelocityEuler(timeStep);
}
//Second, the Euler-Cromer method, also known as the semi-implicit Euler method. 
//Thanks to the way we implemented the position and velocity functions, this is very easy to implement.
void Planet::updateEulerCromer(const planetArray_t& inPlanets, double timeStep) {
	updateAccelerationEuler(inPlanets);
	updateVelocityEuler(timeStep);						//Note the swapped order of velocity and position functions vs Euler. This is because while the velocity function is the same,
	updatePositionEuler(timeStep);						//the position equation wants the (n+1)th velocity, rather than the nth.
}