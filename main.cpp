#include <iostream>

#include <json/json.h>
#include <fstream>
#include <string>

#include <Simulator.h>


int main(int argc, char* argv[]){


    std::cout << "Welcome to the Orbit Simulator" <<std::endl;

    // Read in inputs
    int tf = std::stoi(argv[1]);
    std::string spacecraft = argv[2];
    std::string central_body = argv[3];

    std::cout << "Loading Json inputs" << std::endl;

    Json::Value spacecraft_json, body_json;
    
    // Load data from json files
    try{
        std::ifstream input1(spacecraft);
        std::ifstream input2(central_body);

        input1 >> spacecraft_json;
        input2 >> body_json;
    }
    catch(const std::exception& e){

        std::cout << "Could not load json files" << std::endl;

        std::cout << e.what() << std::endl;
        return 1;
    }
    
    // Call simulator class
    std::cout << "Simulating spacecraft " << spacecraft_json["name"].asString() << " for " << std::to_string(tf) << "s." << std::endl;
    Simulator sim(tf, spacecraft_json, body_json);

    sim.simulate();

    // Call this function in python and generate some basic plots

    return 0;
}