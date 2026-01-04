#include "Simulator.h"

#include <iostream>
#include <json/json.h>
#include <fstream>
#include <string>
#include <cmath>

#include <filesystem>
#include <sys/stat.h>

#include <iomanip>

//#define BOOST_MATH_DISABLE_THREADS  // Disable threading in Boost.Math

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <Functions.h>

#include "Nrlmsise00.hpp"
#include <ctime>

using namespace boost::numeric::odeint;


Simulator::Simulator(double tf, Json::Value spacecraft, Json::Value central_body)
: atmos_flags(), nrlmsise00(atmos_flags)
{
    atmos_flags.fill(0);


    this->tf = tf;
    this->spacecraft = spacecraft;
    this->central_body = central_body;

    mu = G*central_body["mass"].asDouble();

    // Set up initial state
    state.resize(6);

    state << spacecraft["initial_condition"]["eci_position"][0].asDouble(),
             spacecraft["initial_condition"]["eci_position"][1].asDouble(),
             spacecraft["initial_condition"]["eci_position"][2].asDouble(),
             spacecraft["initial_condition"]["eci_velocity"][0].asDouble(),
             spacecraft["initial_condition"]["eci_velocity"][1].asDouble(),
             spacecraft["initial_condition"]["eci_velocity"][2].asDouble();

    
    abs_tol = spacecraft["abs_tol"].asDouble();
    rel_tol = spacecraft["rel_tol"].asDouble();

    burn_counter = 0;
    num_burns = spacecraft["burns"].size();

    // Use ctime to get the day of the year for the current unix timestep
    time_t total_seconds = static_cast<time_t>(spacecraft["date"].asDouble());
    struct tm *now = std::gmtime(&total_seconds);
    
    // Create a copy of the time struct and reset to Jan 1st 00:00:00
    struct tm start_of_year = *now;
    start_of_year.tm_mon = 0;      // January
    start_of_year.tm_mday = 1;     // 1st
    start_of_year.tm_hour = 0;
    start_of_year.tm_min = 0;
    start_of_year.tm_sec = 0;



    #ifdef _WIN32
        this->year_start_unix = static_cast<double>(_mkgmtime(&start_of_year));
    #else
        this->year_start_unix = static_cast<double>(timegm(&start_of_year));
    #endif

    
}

void Simulator::simulate(){

    // Define rk45 solver to work with Eigen library vectors
    typedef runge_kutta_dopri5<Eigen::VectorXd, double, Eigen::VectorXd, double, vector_space_algebra> error_stepper_type;

    // Wrap the ode function in a lambda
    auto ode_func = [this](const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t) {
        this->ode_function(x, dxdt, t);
    };


    // Integrate the ode function until tf or until collision event (detected in observe)
    // Also record states vs. time through observe function
    try{
        integrate_adaptive( make_controlled< error_stepper_type >( abs_tol , rel_tol ) , 
                    ode_func , state , 0.0 , tf , 0.01, [this](Eigen::VectorXd &state , double t ) {
                                this->observe( state , t );
                            });
                        }
    
    catch(const std::runtime_error& e){

        std::cout << "Spacecraft impacted the central body surface." << std::endl;

    }

    // Wtite outputs to csv
    write_output(time, states, derived_states);
}

void Simulator::observe(Eigen::VectorXd &state , double t){
    
    // Build derived state from current state
    Eigen::VectorXd derived_state = build_derived_state(state, t);

    // Save current time, state, and derived states into respective vectors 
    derived_states.push_back(derived_state);
    states.push_back( state );
    time.push_back( t );

    // Check for collision event
    if (std::hypot(state[0], state[1], state[2]) < central_body["radius"].asDouble()){
        throw std::runtime_error("Spacecraft impacted the central body surface.");
    }

    // Check for burn events
    if(t >= spacecraft["burns"][burn_counter]["time"].asDouble() && burn_counter < num_burns){

        std::cout << "Execuring burn at [" << spacecraft["burns"][burn_counter]["delta_v_rtn"][0].asDouble() <<
        ", " << spacecraft["burns"][burn_counter]["delta_v_rtn"][1].asDouble() << ", " <<
        spacecraft["burns"][burn_counter]["delta_v_rtn"][2].asDouble() << "] (rtn) km/s at t = " << 
        t << std::endl;

        // Convert burn in rtn frame to eci frame
        std::vector<double> delta_v_rtn_vec = {
        spacecraft["burns"][burn_counter]["delta_v_rtn"][0].asDouble(),
        spacecraft["burns"][burn_counter]["delta_v_rtn"][1].asDouble(),
        spacecraft["burns"][burn_counter]["delta_v_rtn"][2].asDouble()
        };

        Eigen::Vector3d delta_v_eci(Functions::rtn_to_eci_delta_v(state, delta_v_rtn_vec));
        
        // Add delta-v to current velocity state
        state[3] += delta_v_eci[0];
        state[4] += delta_v_eci[1];
        state[5] += delta_v_eci[2];

        burn_counter++;
    }
}

// Function that gets called on every simulation loop
void Simulator::ode_function(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double t){

    // Initialize derivitive vector
    dxdt = Eigen::VectorXd(x) * 0.0;

    // Define utility variables
    Eigen::Vector3d r_vec(x[0], x[1], x[2]);
    Eigen::Vector3d v_vec(x[3], x[4], x[5]);
    double r = r_vec.norm();

    // Calculate velocity of spacecraft relative to spinning earth
    Eigen::Vector3d omega_earth(0, 0, 7.292115e-5);
    Eigen::Vector3d v_rel = v_vec - omega_earth.cross(r_vec);
    

    // DRAG acceleration
    // Current total Unix time
    double current_unix = spacecraft["date"].asDouble() + t;
    
    // Calculate fractional DOY (0.0 = Jan 1 00:00:00)
    int doy = (current_unix - this->year_start_unix) / 86400;
    
    // Build NRLMSISE-00 inputs 
    std::vector<double> atmosphere_state = build_atmosphere_state(x, r, t);

    // TODO: Input all geomagnetic AP conditions
    std::array<double, 7> aps;
    aps.fill(atmosphere_state[5]);

    // Call NRLMSISE-00 class to output density, ignore if above alt limit
    double density = 0.0;

    if(atmosphere_state[1] < 999){

        density = nrlmsise00.density(doy, atmosphere_state[0],
                  atmosphere_state[1], atmosphere_state[2], atmosphere_state[3],
                atmosphere_state[4], atmosphere_state[5], aps); // density in kg/m^3
        }

    // Project drag force to negative relative velocity vector direction
    Eigen::Vector3d F_drag = -0.5*spacecraft["reference_area"].asDouble()*spacecraft["cd"].asDouble()*density*std::pow(v_rel.norm(), 2)*v_rel.normalized();
    
    std::cout << F_drag << std::endl;
    // Divide by satellite mass for total drag acceleration
    Eigen::Vector3d a_drag = F_drag / spacecraft["mass"].asDouble();

    // J2 acceleration
    double J2_coefficient = (1.5*central_body["J2"].asDouble()*mu*std::pow(central_body["equatorial_radius"].asDouble(), 2)) / std::pow(r, 5);

    double J2_x = J2_coefficient*x[0]*((5*std::pow(x[2], 2) / std::pow(r, 2)) - 1);
    double J2_y = J2_coefficient*x[1]*((5*std::pow(x[2], 2) / std::pow(r, 2)) - 1);
    double J2_z = J2_coefficient*x[2]*((5*std::pow(x[2], 2) / std::pow(r, 2)) - 3);

    // TODO: Add gravitational effect of Moon and Jupiter
    // TODO: Add solar radiation pressure

    // Add perturbing forces to two-body problem and integrate
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];

    dxdt[3] = (-(mu / std::pow(r, 3)) * x[0]) + a_drag[0] + J2_x;
    dxdt[4] = (-(mu / std::pow(r, 3)) * x[1]) + a_drag[1] + J2_y;
    dxdt[5] = (-(mu / std::pow(r, 3)) * x[2]) + a_drag[2] + J2_z;

}

// Function that calculates derived state values on each simulation loop
// TODO: Convert angular orbital parameters to degrees
Eigen::VectorXd Simulator::build_derived_state(Eigen::VectorXd state, double t){

    // DERIVED STATES TO CALCULATE
    // a,e,i,raan, omega, f, E, M, n, p, h, flight path angle
    Eigen::Vector3d r_vec(state[0], state[1], state[2]);
    Eigen::Vector3d v_vec(state[3], state[4], state[5]);

    double r = r_vec.norm();  
    double v = v_vec.norm();
    // Specific energy (vis-viva) equation 
    double Energy = (std::pow(v, 2) / 2) - (mu / r);

    // Semi-major axis
    double a = -mu / (2*Energy);
    
    // Mean motion
    double n = std::sqrt(mu / pow(a, 3));

    // Orbital Period
    double T = 2*M_PI/n;

    // Angular Momentum Vector
    Eigen::Vector3d h_vec = r_vec.cross(v_vec);
    double h = h_vec.norm();

    // Eccentricity Vector
    Eigen::Vector3d e_vec = ((1/mu) * (v_vec.cross(h_vec))) - (r_vec/r);
    double e = e_vec.norm();

    // Semi-latus rectrum
    double p = std::pow(h, 2) / mu;

    // Semi-minor axis
    double b = a*std::sqrt(1 - pow(e, 2));

    // Radii of apogee and perigee
    double ra = a*(1+e);
    double rp = a*(1-e);

    // True Anomaly
    double f = std::acos(((p / r) - 1) / e);

    if (r_vec.dot(v_vec) < 0) {
        f = 2.0 * M_PI - f;
    }
    // Eccentric Anomaly
    double E = 2*std::atan(std::sqrt((1-e) / (1+e))*std::tan(f/2));

    // Mean Anomaly
    double M = E - e*std::sin(E);

    // Flight Path Angle
    double gamma = std::atan((e*std::sin(f)) / (1+ e*std::cos(f)));

    // Inclination
    Eigen::Vector3d x(1, 0, 0);
    Eigen::Vector3d y(0, 1, 0);
    Eigen::Vector3d z(0, 0, 1);

    double i = std::acos((h_vec / h).dot(z));

    // Longitude of Ascending Node
    double laan = std::acos(x.dot(z.cross(h_vec / h)) / std::sin(i));

    if(laan < 0){
        laan += 2*M_PI;
    }

    // Argument of Periapsis
    Eigen::Vector3d node_vector = z.cross(h_vec);
    double omega = std::acos(node_vector.dot(e_vec) / (e*node_vector.norm())); 
    
    if (e_vec.z() < 0){
        omega = 2.0 * M_PI - omega;
    }

    // LLA
    std::vector<double> lla_calc = lla(state, r, t);
    double lat = lla_calc[0];
    double lon = lla_calc[1];
    double alt = lla_calc[2];

    Eigen::VectorXd derived_state(28);

    derived_state << v, r, Energy, a, n, T, h, h_vec[0], h_vec[1], h_vec[2], 
                 e, e_vec[0], e_vec[1], e_vec[2], p, ra, rp, 
                 b, f, E, M, gamma, i, laan, omega,
                 lat, lon, alt;
 
    return derived_state;

}


// Function that builds nrlmsise-00 inputs at every time step
std::vector<double> Simulator::build_atmosphere_state(Eigen::VectorXd state, double r, double t){

    /*
        void CNrlmsise00::gtd7d(const int doy, const double sec, const double& alt,
                 const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t)

    */

    double sec = spacecraft["initial_condition"]["seconds"].asInt() + t;

    // Calculate lla
    std::vector<double> lla_calc = lla(state, r, t);

    // Convert lat/lon to degrees
    double alt = lla_calc[2];
    double lat = lla_calc[0] * 180 / M_PI;
    double lon = lla_calc[1] * 180 / M_PI;


    // Get solar profile
    double f107a = spacecraft["initial_condition"]["f107a"].asDouble();
    double f107 = spacecraft["initial_condition"]["f107"].asDouble();
    double ap = spacecraft["initial_condition"]["ap"].asDouble();


    return { sec, alt, lat, lon, f107a, f107, ap } ;

} 

std::vector<double> Simulator::lla(Eigen::VectorXd state, double r, double t){

    
    double julian_days = ((spacecraft["date"].asDouble() + t) / 86400) + 2440587.5;
    double era = 2*M_PI*(0.7790572732640 + 1.00273781191135448*(julian_days - 2451545.0));

    era = std::fmod(era, 2*M_PI);

    // Calculate geocentric latitude and longitude
    double lat_geocentric = std::asin(state[2] / r);
    double lon = std::atan2(state[1], state[0]) - era;

    // Calculate geodetic latitude
    double lat = std::atan((std::pow(central_body["equatorial_radius"].asDouble(), 2) / std::pow(central_body["polar_radius"].asDouble(), 2))*std::tan(lat_geocentric));

    
    // Shift longitude to be clocked from International Date Line
    lon = std::fmod(lon, 2*M_PI);
    if(lon < 0){
        lon += 2*M_PI;
    }
    // WGS84 model for earth's oblate radius 
    double oblate_radius = std::sqrt((std::pow(std::pow(central_body["equatorial_radius"].asDouble(), 2)*std::cos(lat), 2) + std::pow(std::pow(central_body["polar_radius"].asDouble(), 2)*std::sin(lat), 2))
                        / (std::pow(central_body["equatorial_radius"].asDouble()*std::cos(lat), 2) + std::pow(central_body["polar_radius"].asDouble()*std::sin(lat), 2)));

    // Calculate altitude
    double alt = r - oblate_radius;

    return {lat, lon, alt};

}
// Function that writes states to output csv
void Simulator::write_output(std::vector<double>& time, 
                  std::vector<Eigen::VectorXd>& states, 
                  std::vector<Eigen::VectorXd>& derived_states
                  ) {
    
    // Filename
    const std::string& filename = spacecraft["name"].asString() + "_output.csv";

    // Hardcoded headers
    const std::vector<std::string> state_headers = {
        "ECI_X", "ECI_Y", "ECI_Z", "Vx", "Vy", "Vz"
    };
    
    const std::vector<std::string> derived_headers = {
        "v", "r", "E", "a", "n", "T", "h", "h_x", "h_y", "h_z", 
        "e", "e_x", "e_y", "e_z", "p", "ra", "rp", 
        "b", "f", "E_anom", "M", "gamma", "i", "laan", "omega",
        "lat", "lon", "alt"
    };
    
    // Create output directory if it doesn't exist
    std::string output_dir = "output";

    std::filesystem::create_directories(output_dir.c_str());


    // Construct full file path
    std::string filepath = output_dir + "/" + filename;
    
    std::ofstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filepath << std::endl;
        return;
    }
    
    // Set precision
    file << std::scientific << std::setprecision(10);
    
    // Check that all vectors have the same length

    std::cout << time.size() << std::endl;
    std::cout << states.size() << std::endl;
    std::cout << derived_states.size() << std::endl;


    if (time.size() != states.size() || time.size() != derived_states.size()) {
        std::cerr << "Error: time, states, and derived_states must have the same length" << std::endl;
        file.close();
        return;
    }
    
    if (time.empty()) {
        std::cerr << "Warning: No data to write" << std::endl;
        file.close();
        return;
    }
    
    // Get dimensions
    int num_states = states[0].size();
    int num_derived = derived_states[0].size();
    
    // Write header row
    file << "time";
    
    // State headers
    for (int i = 0; i < num_states; i++) {
        file << ",";
        if (i < state_headers.size()) {
            file << state_headers[i];
        } else {
            file << "state_" << i;
        }
    }
    
    // Derived state headers
    for (int i = 0; i < num_derived; i++) {
        file << ",";
        if (i < derived_headers.size()) {
            file << derived_headers[i];
        } else {
            file << "derived_" << i;
        }
    }
    file << "\n";
    
    // Write data rows
    for (size_t row = 0; row < time.size(); row++) {
        // Write time
        file << time[row];
        
        // Write state values
        for (int i = 0; i < states[row].size(); i++) {
            file << "," << states[row](i);
        }
        
        // Write derived state values
        for (int i = 0; i < derived_states[row].size(); i++) {
            file << "," << derived_states[row](i);
        }
        
        file << "\n";
    }
    
    file.close();
    std::cout << "Data written to " << filepath << std::endl;
}