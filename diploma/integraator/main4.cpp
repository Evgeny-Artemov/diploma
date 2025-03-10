#include <iostream>
#include <vector>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <functional>
#include <iomanip>

using namespace std;
using namespace boost::numeric::odeint;

// Define a type for the system state
typedef std::vector<double> state_type;

// Global constants for unit conversion
const double R = 149598261.0; // Reference length (Earth's semi-major axis) in km
const double mu_dim = 132712440018.0; // Dimensional gravitational parameter
const double V = sqrt(mu_dim / R); // Reference velocity in km/s
const double T = R / V; // Reference time in seconds
const double m0 = 1000.0; // Initial spacecraft mass (kg)
const double P_dim = 0.22; // Thrust acceleration in m/s² (converted to km/s² below)
const double Isp = 2000.0; // Specific impulse in seconds

// Observer function to print results in dimensional units
void observer(const state_type &x, const double t) {
    // Convert to dimensional units for output
    double dim_t = t * T / 86400.0; // Time in days
    double dim_x = x[0] * R;
    double dim_y = x[1] * R;
    double dim_z = x[2] * R;
    double dim_vx = x[3] * V;
    double dim_vy = x[4] * V;
    double dim_vz = x[5] * V;
    double dim_m = m0 - x[6]; // Current mass (initial mass - consumed)
    
    // Calculate thrust status
    double lambda_v_norm = sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]);
    double P = (P_dim / 1000.0) * T * T / R; // Convert thrust from m/s² to non-dim units
    double W = Isp * 9.81 / V; // Convert Isp to non-dim units
    double switching_function = 1.0 - (P * lambda_v_norm) / (W * fabs(x[13] + 1e-10));
    string thrust_status = (switching_function < 0) ? "ON" : "OFF";
    
    cout << "Time = " << fixed << setprecision(2) << dim_t << " days" << endl;
    cout << "Position (km) = (" << scientific << setprecision(6) 
         << dim_x << ", " << dim_y << ", " << dim_z << ")" << endl;
    cout << "Velocity (km/s) = (" << fixed << setprecision(6) 
         << dim_vx << ", " << dim_vy << ", " << dim_vz << ")" << endl;
    cout << "Current Mass (kg) = " << fixed << setprecision(4) << dim_m << endl;
    cout << "Consumed Fuel (kg) = " << fixed << setprecision(4) << x[6] << endl;
    cout << "Thrust Status: " << thrust_status << endl;
    cout << "Switching function value: " << switching_function << endl;
    cout << "Lambda values: " << scientific << x[10] << ", " << x[11] << ", " << x[12] << ", " << x[13] << endl;
    cout << "----------------------------------------" << endl;
}

int main() {
    // Initialize the state vector with 14 elements
    state_type x(14);
    
    // Initialize position and velocity (non-dimensionalized)
    x[0] = -1.321562505626410E+08 / R; // Initial x position
    x[1] = 6.658171852356435E+07 / R;  // Initial y position
    x[2] = -2.596377182275057E+03 / R; // Initial z position
    x[3] = -1.389857140256131E+01 / V; // Initial Vx
    x[4] = -2.671292464300294E+01 / V; // Initial Vy
    x[5] = 1.718741410194013E-03 / V;  // Initial Vz
    x[6] = 0.0; // Initial consumed fuel mass (zero at start)
    
    // Initialize conjugate variables with non-zero values to enable optimal control
    // These values need to be tuned for your specific mission
    x[7] = 0.0;    // λx
    x[8] = 0.0;    // λy
    x[9] = 0.0;    // λz
    x[10] = -1.0;  // λVx - non-zero to enable thrust
    x[11] = -1.0;  // λVy - non-zero to enable thrust
    x[12] = 0.0;   // λVz
    x[13] = -0.1;  // λm - non-zero and negative to enable thrust
    
    cout << "Optimal Control Problem for Earth-Pandora Spacecraft Trajectory" << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Initial Position (km) = (" << scientific << setprecision(6) 
         << -1.321562505626410E+08 << ", " 
         << 6.658171852356435E+07 << ", " 
         << -2.596377182275057E+03 << ")" << endl;
    cout << "Initial Velocity (km/s) = (" << fixed << setprecision(6) 
         << -1.389857140256131E+01 << ", " 
         << -2.671292464300294E+01 << ", " 
         << 1.718741410194013E-03 << ")" << endl;
    cout << "Initial Mass (kg) = " << m0 << endl;
    cout << "Thrust Acceleration = " << P_dim << " m/s²" << endl;
    cout << "Specific Impulse = " << Isp << " seconds" << endl;
    cout << "--------------------------------------------------------" << endl;
    
    // Define the stepper type
    typedef runge_kutta_fehlberg78<state_type> stepper_type;

    // Define the system of differential equations
    auto system = [](const state_type &x, state_type &dxdt, double t) {
        // Non-dimensionalized parameters
        const double mu = 1.0; // Non-dimensionalized gravitational parameter
        const double P = (P_dim / 1000.0) * T * T / R; // Convert thrust from m/s² to non-dim units
        const double W = Isp * 9.81 / V; // Convert Isp to non-dim units
        
        // Calculate thrust direction (unit vector)
        double lambda_v_norm = sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]);
        
        // Initialize thrust control and thrust direction
        double delta = 0.0;
        double p0_x = 0.0, p0_y = 0.0, p0_z = 0.0;
        
        if (lambda_v_norm > 1e-10) {
            // Calculate thrust direction (opposite to velocity adjoint vector)
            p0_x = -x[10] / lambda_v_norm;
            p0_y = -x[11] / lambda_v_norm;
            p0_z = -x[12] / lambda_v_norm;
            
            // Apply Pontryagin's Maximum Principle with proper switching function
            double switching_function = 1.0 - (P * lambda_v_norm) / (W * fabs(x[13] + 1e-10));
            delta = (switching_function < 0) ? 1.0 : 0.0;
        }
        
        // Calculate radius vector magnitude
        double r_magnitude = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        double r_magnitude_cubed = r_magnitude * r_magnitude * r_magnitude;
        if (r_magnitude_cubed < 1e-10) r_magnitude_cubed = 1e-10;
        
        // Calculate current mass
        double current_mass = m0 - x[6];
        if (current_mass < 1e-10) current_mass = 1e-10;
        
        // Equations of motion
        dxdt[0] = x[3]; // dx/dt = Vx
        dxdt[1] = x[4]; // dy/dt = Vy
        dxdt[2] = x[5]; // dz/dt = Vz
        
        // Velocity equations with gravity and thrust
        dxdt[3] = -mu * x[0] / r_magnitude_cubed + delta * (P * m0 / current_mass) * p0_x;
        dxdt[4] = -mu * x[1] / r_magnitude_cubed + delta * (P * m0 / current_mass) * p0_y;
        dxdt[5] = -mu * x[2] / r_magnitude_cubed + delta * (P * m0 / current_mass) * p0_z;
        
        // Mass consumption equation - THIS IS THE KEY FIX
        dxdt[6] = delta * (P * m0 / (9.81 * Isp));
        
        // Conjugate variable equations for position
        double r_magnitude_5 = r_magnitude_cubed * r_magnitude * r_magnitude;
        dxdt[7] = mu * (3 * x[0] * (x[10]*x[0] + x[11]*x[1] + x[12]*x[2]) / r_magnitude_5 - x[10] / r_magnitude_cubed);
        dxdt[8] = mu * (3 * x[1] * (x[10]*x[0] + x[11]*x[1] + x[12]*x[2]) / r_magnitude_5 - x[11] / r_magnitude_cubed);
        dxdt[9] = mu * (3 * x[2] * (x[10]*x[0] + x[11]*x[1] + x[12]*x[2]) / r_magnitude_5 - x[12] / r_magnitude_cubed);
        
        // Conjugate variable equations for velocity
        dxdt[10] = -x[7];
        dxdt[11] = -x[8];
        dxdt[12] = -x[9];
        
        // Conjugate variable equation for mass - FIXED
        dxdt[13] = delta * P * m0 * m0 * (x[10]*p0_x + x[11]*p0_y + x[12]*p0_z) / (current_mass * current_mass);
    };

    // Create a controlled stepper with specified error tolerances
    auto stepper = make_controlled<stepper_type>(1.0e-10, 1.0e-10);

    // Perform integration with more frequent output
    integrate_adaptive(
        stepper,
        system,
        x,
        0.0,    // Start time
        0.5,    // End time (shorter period to see changes)
        0.05,   // Initial step size (smaller for more frequent updates)
        observer
    );

    // Print final summary
    cout << "\nFinal State Summary:" << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Final Position (km) = (" << scientific << setprecision(6) 
         << x[0]*R << ", " << x[1]*R << ", " << x[2]*R << ")" << endl;
    cout << "Final Velocity (km/s) = (" << fixed << setprecision(6) 
         << x[3]*V << ", " << x[4]*V << ", " << x[5]*V << ")" << endl;
    cout << "Final Mass (kg) = " << fixed << setprecision(4) << (m0 - x[6]) << endl;
    cout << "Total Fuel Consumed (kg) = " << fixed << setprecision(4) << x[6] << endl;

    return 0;
}