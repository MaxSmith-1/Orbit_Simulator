#include <iostream>
#include <json/json.h>

#include <fstream>
#include <string>

#include <Eigen/Dense>
#include <vector>

// Header only class, just for the functions

class Functions{

    public:

    static Eigen::Vector3d rtn_to_eci_delta_v(Eigen::VectorXd state, const std::vector<double>& delta_v_rtn){
        

        // Calculate r hat vector
        Eigen::Vector3d r_vec(state[0], state[1], state[2]);
        Eigen::Vector3d r_hat = r_vec / r_vec.norm();
        
        // Calculate n hat vector
        Eigen::Vector3d v_vec(state[3], state[4], state[5]);
        Eigen::Vector3d h_vec = r_vec.cross(v_vec);
        Eigen::Vector3d n_hat = h_vec / h_vec.norm();

        // Calculate t hat vector
        Eigen::Vector3d t_hat = n_hat.cross(r_vec / r_vec.norm());

        // Assemble rtn to eci rotation matrix
        Eigen::Matrix3d RTN_to_ECI;
        RTN_to_ECI << r_hat[0], t_hat[0],   n_hat[0],
            r_hat[1],    t_hat[1], n_hat[1],
            r_hat[2],      t_hat[2],     n_hat[2];
        
        // Convert delta v vector to an eigen library vector
        Eigen::Vector3d delta_v = Eigen::Map<const Eigen::Vector3d>(delta_v_rtn.data());

        // Return rotation product
        return RTN_to_ECI * delta_v;

    }


};

