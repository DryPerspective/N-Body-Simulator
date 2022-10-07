// SolarSystem.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <string_view>	//For more efficient "views" of strings.
#include <charconv>		//To read string_views into numbers
#include <bitset>		//Used to track properly initialised components of a planet.
#include <array>		//Used to track how far along the simulation is


#include "PhysicsVector.h"
#include "Planet.h"

//To prevent confusion between a vector, the mathematical object of a number with direction, and std::vector, we use this alias.
using planetArray_t = std::vector<Planet>;
//To save having to dig through the utilities library namespaces we use an alias for our vector
using vector3D_t = Physics::PhysicsVector<3>;


//This function adds the default solar system to the simulation. Planetary data courtesy of NASA JPL.
//This function is only called if after reading the input file, the simulation contains no planets.
void importDefaultData(planetArray_t& inPlanets) {
	inPlanets.push_back(Planet("The Sun", 1.989e30, { 0,0,0 }, { 1.998619875971241, 1.177175852520643e1,-6.135600299763972e-2 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Mercury", 3.3011e23, { 1.275387239870491E+10,-6.680195324480709E+10,-6.616376210554786E+09 }, { 3.815800795678611E+04,1.123692837720359E+04,-2.583452372780768E+03 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Venus", 4.867e24, { -8.073224723501202E+10,7.027586666429530E+10,5.627818208653621E+09 }, { -2.299827401900994E+04,-2.669115882767952E+04,9.610940692989782E+02 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Earth", 5.972e24, { 4.788721549926552E+10,1.398390053760727E+11,-2.917617879798263E+07 }, { -2.869322295421606E+04,9.472398427890313E+03,-1.294094780725619E-00 }, { 0,0,0 }));
	inPlanets.push_back(Planet("The Moon", 734.9e20, { 4.749196053391321E+10,1.399182076993898E+11,-3.486943982706219E+07 }, { -2.890724003060377E+04,8.531016069261970E+03,8.300527233703736E+01 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Mars", 6.4171e23, { -2.360304784158461E+11,7.782743203688863E+10,7.409494561464485E+09 }, { -6.646816636079097E+03,-2.094094408471671E+04,-2.759397656641038E+02 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Jupiter", 1.89813e27, { -7.635337060440624E+11,2.666352191711917E+11, 1.596697237644111E+10 }, { -4.459151830811911E+03,-1.171879602036105E+04,1.485480013373461E+02 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Saturn", 5.68319e26, { -5.754602000703751E+11,-1.380800977297312E+12,4.691113811667019E+10 }, { 8.388118620089763E+03,-3.745812490969359E+03,-2.682504240279582E+02 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Uranus", 86.8103e24, { 2.828705362370189E+12, 9.657796340541244E+11,-3.305961929341555E+10 }, { -2.249907923122420E+03,6.127203368970902E+03,5.166083013695255E+01 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Neptune", 102.41e24, { 4.177286553745139E+12,-1.624410031732890E+12,-6.281810904534376E+10 }, { 1.934495516018552E+03,5.098519902111810E+03,-1.496666233625485E+02 }, { 0,0,0 }));
	inPlanets.push_back(Planet("Pluto", 1.308e22, { 1.263871593868758E+12,-4.769395770475431E+12,1.447666788459496E+11 }, { 5.347856858111191E+03,2.674281760600502E+02,-1.564505494419083E+03 }, { 0, 0, 0 }));
}



//This method calculates the position of the center of mass of the system. This is given by Sum(mass_n * position_n)/Sum(mass_n). 
//Since we're dealing with vectors, we need to do this on a per-component basis.
vector3D_t centreOfMass(const planetArray_t& Planets) {
	double comX{ 0 };
	double comY{ 0 };
	double comZ{ 0 };
	double totalMass{ 0 };
	for (auto planet : Planets) {
		comX += planet.getPosition().getX() * planet.getMass();		//Sum the individual mass*position terms for the numerator
		comY += planet.getPosition().getY() * planet.getMass();
		comZ += planet.getPosition().getZ() * planet.getMass();
		totalMass += planet.getMass();								//Sum the masses for the denominator
	}
	vector3D_t centreOfMass{ comX / totalMass, comY / totalMass, comZ / totalMass };
	return centreOfMass;	
}

//This function makes use of std::from_chars to read a double value from a string_view.
double readChars(const std::string_view& inString) {
	double outputNumber;

	//The format. from_chars needs to know ahead of time whether the value it's reading is in scientific format.
	//So we do a simple check - if the string contains an e then we assume it's in scientific. Invalid input is handled later.
	std::chars_format format;
	if (inString.find('e') == std::string_view::npos)format = std::chars_format::general;
	else format = std::chars_format::scientific;

	//And now to grab our result. Note that the actual value we're looking for is stored in outputNumber, whereas the result of from_chars is a struct containing data about how the process went.
	const auto result{ std::from_chars(inString.data(), inString.data() + inString.length(), outputNumber, format) };

	//Now onto error handling and ensuring that we got the right result.
	//As the success largely depends on the user entering the correct data in the input file, we throw exceptions when they do not.
	if (result.ptr == inString.data() + inString.length()) return outputNumber;
	else if (result.ec == std::errc::invalid_argument) {
		std::cerr << "Error in config file. Value: " << inString << " follows invalid format!";
		throw std::invalid_argument("Error: invalid line format in config.txt");
	}
	else if (result.ec == std::errc::result_out_of_range) {
		std::cerr << "Error in config file. Value " << inString << " goes out of range. Distance: " << result.ptr - inString.data() << '\n';
		throw std::range_error("Error: config.txt value outside of double range.");
	}


	return -1;
}

//This vector reads a string of the form "(e1, e2, e3)" and transforms it into a PhysicsVector object. To save having to pass that object by value, we pass it in and modify it in place.
//This seems like the best solution as the only application for this function is to modify already existent vector objects.
void readVector(std::string_view inString, vector3D_t& finalVector) {


	//First trim the brackets from the string, if they exist
	if (inString[0] == '(')inString.remove_prefix(1);
	if (inString[inString.length() - 1] == ')')inString.remove_suffix(1);

	//Then we delimit by the comma. If this fails we throw an exception.
	//First we ensure that there are the expected two commas.
	auto numberOfCommas{ std::count(inString.begin(),inString.end(),',') };
	if (numberOfCommas != 2) {
		std::cerr << "Error in config file. Line: " << inString << " does not contain the correct amount of commas to be read as a 3D vector";
		throw std::invalid_argument("Error: Expecting two commas to read a 3D vector.");
	}
	else {
		//Otherwise, we delimit at the first comma.
		auto firstComma{ inString.find_first_of(',') };
		std::string_view firstTerm{ inString.substr(0,firstComma) };
		double e1{ readChars(firstTerm) };
		
		//With term 1 extracted, we can remove it (and the first comma) from the string and focus on the second term
		inString.remove_prefix(firstComma+1);
		auto secondComma{ inString.find_first_of(',') };
		std::string_view secondTerm{ inString.substr(0,secondComma) };
		double e2{ readChars(secondTerm) };

		//Finally, we can remove the second term and all we will be left with in the string view is the third term.
		inString.remove_prefix(secondComma+1);
		double e3{ readChars(inString) };

		//Then we just set the vector passed in.
		finalVector = { e1,e2,e3 };

	}
}


int main()
{
	//The simulation configuration variables, measured in seconds.
	double timeStep{ 1 };		
	double totalLength{ 10 };	

	planetArray_t Planets{};

	//We avoid using the ConfigReader object from the Basic Utilities library, as it does not support multiple variables with the same name in the config file, and the config file
	//is more readable divided into easy blocks of name, mass, position, velocity, for each planet rather than having to contrive and hard-code a unique identifier for each planet's properties.

	//The data read into the simulation comes from a file called config.txt
	std::ifstream configFile("config.txt");

	//Next, we want to read through the config file to set our simulation variables and read in all our planets.
	//We need to declare a lot of variables ahead of time as we don't want them reset every loop.
	std::string inputLine;									//Each line of the file is stored in this variable
	std::bitset<4> initialisedComponents{ "0000" };			//This is used to track whether all four components needed to construct a planet have been initialised.
	std::string newName;									//The name of each new planet. Stored as a string - storing it as a string_view would change its content if the string it's viewing changes.
	double newMass;											//Its mass
	vector3D_t newPos;										//Position
	vector3D_t newVel;										//And velocity.

	//Then we loop through the file.
	while (getline(configFile, inputLine)) {
		//First we need to do some basic processing of our lines to get them into values we can read.
		inputLine.erase(std::remove_if(inputLine.begin(), inputLine.end(), isspace), inputLine.end());		//First trim off the whitespace
		if (inputLine[0] == '#' || inputLine.empty())continue;												//Then ignore commented and empty lines
		
		
		//Split the line using = as a delimiter.
		auto splitPos{ inputLine.find('=') };																

		//Then use string_view to get two views of the line - one for each side of the delimiter
		std::string_view lineBeforeEquals{ inputLine };														
		lineBeforeEquals.remove_suffix(inputLine.length() - splitPos);
		std::string_view lineAfterEquals{ inputLine };
		lineAfterEquals.remove_prefix(splitPos+1);
		
		//Once we have separated out our lines, we can start processing them. We start with our simulation constants.
		if (lineBeforeEquals == "timeStep") timeStep = readChars(lineAfterEquals);
		else if (lineBeforeEquals == "simulationLength")totalLength = readChars(lineAfterEquals);
		//If we get this far we are probably creating a new planet.			
		else if (lineBeforeEquals == "name") {
			newName = lineAfterEquals;
			initialisedComponents.set(0, true);
		}
		else if (lineBeforeEquals == "mass") {
			newMass = readChars(lineAfterEquals);
			initialisedComponents.set(1, true);
		}
		else if (lineBeforeEquals == "position") {
			readVector(lineAfterEquals, newPos);
			initialisedComponents.set(2, true);
		}
		else if (lineBeforeEquals == "velocity") {
			readVector(lineAfterEquals, newVel);
			initialisedComponents.set(3, true);
		}
		//If we skip past those and we can't identify what lineBeforeEquals says, we have a problem.
		else {
			std::cerr << "Error in config file: Line " << lineBeforeEquals << " does not match an expected value.\n";
			throw std::invalid_argument("Error in config file: Invalid expression");
		}

		

		//If all four planet variables have been properly set, our initialisedComponents bitset will be all true.
		if (initialisedComponents.all()) {			
			Planets.push_back(Planet(newName, newMass, newPos, newVel));	//So we make the new planet
			initialisedComponents.reset();									//And reset the bitset.
		}

	}

	std::cout << "Simulation time step : " << timeStep << '\t' << "Simulation total simulated length: "<<totalLength << '\n';

	if (Planets.size() == 0) {
		std::cout << "Planets could not be read from config.txt, or config.txt is empty. Adding default solar system...\n";
		importDefaultData(Planets);
		std::cout << "Default solar system loaded.\n";
	}
	else {
		std::cout << "Planets being simulated: " << Planets.size() << '\n';
	}

	//Create our outputfile
	std::string outputFileName{ "cppOutputFile.csv" };
	std::ofstream outputFile(outputFileName);

	//Write column headers to the output file
	for (auto& planet : Planets) {
		outputFile << planet.getName() << "X," << planet.getName() << "Y," << planet.getName() << "Z,";
	}
	outputFile << '\n';

	//As we are potentially simulating a lot of planets over a long period of time, it might be nice to know how far along the simulation is.
	std::array<double, 100> percentageMarkers;	//A measure of how far along the simulation is. entry [0] -> 1%, [1] -> 2% etc.
	std::array<bool, 100> hasbeenPrinted;		//A collection of bools to ensure that each percentage marker is only printed once.

	//Fill them with initial values
	for (int i = 0; i < 100; ++i) {
		percentageMarkers[i] = i * (totalLength / 100);
		hasbeenPrinted[i] = false;
	}

	std::cout << "Beginning simulation.\n";

	int currentPercent{ 0 };	//Used as a tracker to prevent needing to search the entire percentageMarkers for how far along we are every run.
	double currentLength{ 0 };
	while(currentLength<totalLength){
		//First, process how far along we are:
		if (currentLength > percentageMarkers[currentPercent] && hasbeenPrinted[currentPercent]==false && currentPercent<99) {	//currentPercent <99 to prevent access violation
			std::cout << currentPercent + 1 << "% complete.\n";
			hasbeenPrinted[currentPercent] = true;
			++currentPercent;
		}




		//In reality, the planets don't orbit the exact center of the sun. They orbit the system's joint center of mass.
		//By far the simplest way to implement this is set the center of mass at the origin of the system, and move everything else in the universe around to accommodate.
		vector3D_t CoM{ centreOfMass(Planets) };
		for (auto& planet : Planets) {
			planet.setPosition(planet.getPosition() - CoM);
		}

		//Update the planet following the Euler Cromer method.
		for (auto& planet : Planets) {
			planet.updateEulerCromer(Planets, timeStep);
		}

		//And write the updated data to the output file.
		for (const auto& planet : Planets) {
			outputFile << planet.getPosition().getX() << "," << planet.getPosition().getY() << ',' << planet.getPosition().getZ() << ',';
		}
		outputFile <<  '\n';
		currentLength += timeStep;
		
	}

	std::cout << "100% complete.\nData written to " << outputFileName << '\n';

	outputFile.close();
	

}

