#ifndef Planet_H
#define Planet_H

#include <string>
#include <cmath> //For squaring

#include "PhysicsVector.h"

/* 
* A planet object for the purposes of gravitational simulation. Potentially a misnomer as it should account for every space object, including the sun and any other non-planets.
* However calling it a "gravity field" or "space object" would be much more vague and wouldn't describe the purpose of this object as accurately.
* The decision was made to measure all quantities in base units (e.g. kg, m, m/s) rather than their scaled up counterparts, as the double type should comfortably contain the numbers,
* and it reduces the risk of errors from scale factoring and other issues.
* 
* Several functions are included to calculate and update a planet's status according to various approximations.
*/

class Planet
{
private:
	static constexpr double G{ 6.67408e-11 };				//The gravitational constant. Unless you're doing weird things with general relativity, this value is universally constant.
															//G is also very difficult to measure precisely, so 5 decimal places is all we get here.
	double			 m_mass{0};								//Planetary mass, measured in kg.
	PhysicsVector<3> m_position;							//Planetary position, measured in m.
	PhysicsVector<3> m_velocity;							//Planetary velocity, measured in m/s.
	PhysicsVector<3> m_acceleration;						//Planetary acceleration, measured in m/s^2.
	std::string		 m_name{"Unnamed Planet"};				//And of course the object's name.

public:
	//Constructors.
	Planet(const std::string& inName, const double inMass, const PhysicsVector<3>& inPos, const PhysicsVector<3>& inVel);
	Planet(const std::string& inName, const double inMass, const PhysicsVector<3>& inPos, const PhysicsVector<3>& inVel, const PhysicsVector<3>& inAcc);
	Planet(const double inMass,const PhysicsVector<3>& inPos,const PhysicsVector<3>& inVel,const PhysicsVector<3>& inAcc);
	Planet(const std::string& inName);

	//Virtual default destructor in case of inheritance.
	virtual ~Planet() = default;

	//Getters and setters
	double getMass() const;
	const PhysicsVector<3>& getPosition() const;
	const PhysicsVector<3>& getVelocity() const;
	const PhysicsVector<3>& getAcceleration() const;
	const std::string& getName() const;

	void setMass(double inMass);
	void setPosition(const PhysicsVector<3>& inPos);
	void setVelocity(const PhysicsVector<3>& inVel);
	void setAcceleration(const PhysicsVector<3>& inAcc);
	void setName(const std::string& inName);

	//To prevent confusion between a vector, the mathematical object of a number with direction, and std::vector, we use this alias.
	using planetArray_t = std::vector<Planet>;
	//Calculate the acceleration on this planet caused by a single other planet.
	PhysicsVector<3> calcAcceleration(const Planet& inPlanet);

	//Update the three quantities following the Euler method.
	void updateAccelerationEuler(const planetArray_t& inPlanets);
	void updatePositionEuler(const double timeStep);
	void updateVelocityEuler(const double timeStep);



	//The core algorithm functions.
	void updateEuler(const planetArray_t& inPlanets, double timeStep);
	void updateEulerCromer(const planetArray_t& inPlanets, double timeStep);



};




#endif