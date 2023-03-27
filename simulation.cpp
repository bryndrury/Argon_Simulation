// Assessment 3: n-body gravitational solver

// To avoid warnings tell the compiler to use a recent standard of C++:
// g++ -std=c++17 vector3d.cpp body.cpp compute_orbits.cpp -o compute_orbits
// ./compute_orbits sun_earth.csv test_case.csv 1e-4 110

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

#include <chrono>
#include <stdlib.h>

#include "vector3d.hpp"
#include "body.hpp"
// #include "import.hpp"

using std::cout, std::endl;

// *** Function declarations ***

// -----------------------------------------------------------
// ADD YOUR FUNCTION DECLARATIONS HERE

void boundary(vec &pos, vec &vel, const double L);
void update_acc(int N, std::vector<body> &system);
void vel_verlet(int N, int L, std::vector<body> &system, double dt);

// -----------------------------------------------------------
// Read input data from file
void read_init(std::string input_file, std::vector<body> &system);
// Read the components of a 3d vector from a line
void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z);
// Save the data to file
void save_data(std::ofstream& savefile, const std::vector<body> &system, double t);

int main(int argc, char* argv[])
{
    // Checking if number of arguments is equal to 4:
    if (argc != 7) {
        cout << "ERROR: need 5 arguments - particle_simulation <input_file> <output_file> <dt> <T> <Tsave> <L>" << endl;
        return EXIT_FAILURE;
    }
    // Process command line inputs:
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    double dt = atof(argv[3]); // Time step
    int T = atoi(argv[4]); // Total number of time steps
    int Tsave = atoi(argv[5]); // Number of steps between each save
    double L = atoi(argv[6]); // The box size

    std::vector<body> system; // Create an empty vector container for bodies
    read_init(input_file, system); // Read bodies from input file into system
    int N = system.size(); // Number of bodies in system

    cout << "--- Simulation of Particle Motion ---" << endl;
    cout << " number of particles N: " << N << endl;
    /*
    for(int p=0; p < N; p++)
    {
      cout << "- " << system[p].get_name() << endl; // Display names
    }
    */
    cout << "       time step dt: " << dt << endl;
    cout << "  number of steps T: " << T << endl;
    cout << "   save steps Tsave: " << Tsave << endl;

    std::ofstream savefile (output_file); // Open save file
    if (!savefile.is_open()) {
        cout << "Unable to open file: " << output_file << endl; // Exit if save file has not opened
        return EXIT_FAILURE;
    }
    savefile << std::setprecision(16); // Set the precision of the output to that of a double
    // Write a header for the save file
    savefile << dt << "," << T << "," << Tsave << "," << N << endl;
    for(int p=0; p < (N-1); p++)
    {
        savefile << system[p].get_name() << ",";
    }

    savefile << system[N-1].get_name() << endl;

    // -----------------------------------------------------------
    // ADD YOUR CODE HERE
    int count = Tsave;
    double time = 0;
    auto begin = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < T+1; t++)
    {
        if (count == Tsave)
        {
            save_data(savefile, system, time);
            count = 1;


           std::system("clear");
           cout << std::setprecision(2);

          auto now = std::chrono::high_resolution_clock::now();
          double est_finish = (std::chrono::duration_cast<std::chrono::seconds>(now-begin).count() / (100*time/T)) - std::chrono::duration_cast<std::chrono::seconds>(now-begin).count();
//               (std::chrono::duration_cast<std::chrono::seconds>(now-begin).count() / (time/T)) - std::chrono::duration_cast<std::chrono::seconds>(now-begin).count();

          cout << endl << endl << "- - - - - - - - - - - - - - - -" << endl;
          cout << "Simulation Progress: " << 100*t/T << "%" << endl;
          cout << std::setprecision(5);
          cout << "Estimated time to Completion: " << est_finish << " Seconds." << endl;
          cout << "- - - - - - - - - - - - - - - -" << endl << endl;
        }
        else {count++;}
        time += dt;

        vel_verlet(N, L, system, dt);
    }

    // -----------------------------------------------------------

    savefile.close();
    return EXIT_SUCCESS;
}

// *** Function implementations ***

// -----------------------------------------------------------
// ADD YOUR FUNCTION IMPLEMENTATIONS HERE
void update_acc(int N, std::vector<body> &system)
{
    double epsilon = 125.7 * 1.38e-23;
    double sigma = 0.3345e-9;

    for (int p = 0; p < N; p++)
    {
        vec acc;

        double m_inv = 1/ ( system[p].get_mass() );

        for (int j = 0; j < N; j++)
        {
            vec distance = ( system[j].get_pos() - system[p].get_pos() );

            double inv_rx = 1 / (distance).x();
            double inv_ry = 1 / (distance).y();
            double inv_rz = 1 / (distance).z();

            double a_x, a_y, a_z = 0;

//             cout << distance.x() << endl;

            if (distance.x() != 0)
            {
                a_x = -(4 * epsilon * m_inv) * ( (12 * pow(sigma, 12) * pow(inv_rx, 13)) - (6 * pow(sigma, 6) * pow(inv_rx, 7)) );
            }
            else a_x = 0;

            if (distance.y() != 0)
            {
                a_y = -(4 * epsilon * m_inv) * ( (12 * pow(sigma, 12) * pow(inv_ry, 13)) - (6 * pow(sigma, 6) * pow(inv_ry, 7)) );
            }
            else a_y = 0;

            if (distance.z() != 0)
            {
                a_z = -(4 * epsilon * m_inv) * ( (12 * pow(sigma, 12) * pow(inv_rz, 13)) - (6 * pow(sigma, 6) * pow(inv_rz, 7)) );
            }
            else a_z = 0;

//             cout << a_x << ", " << a_y << ", " << a_z << endl;

            vec a (a_x, a_y, a_z);

            acc += a;
        }
//         cout << " Acceleration: " << acc;
        system[p].set_acc(acc);
    }
}

void vel_verlet(int N, int L, std::vector<body> &system, double dt)
{
    std::vector<body> new_system = system;

    for (int p = 0; p < N; p++)
    {
        vec new_r;
        vec new_v;

        update_acc(N, new_system);

        new_r = system[p].get_pos() + ( system[p].get_vel() * dt ) + ( 0.5 * system[p].get_acc() * pow(dt, 2) );

        new_v = system[p].get_vel() + ( 0.5 * dt * ( new_system[p].get_acc() + system[p].get_acc() ) );

        boundary(new_r, new_v, L);

        new_system[p].set_pos(new_r);
        new_system[p].set_vel(new_v);
    }
    system = new_system;
}

void boundary(vec &pos, vec &vel, const double L)
{
//     cout << L << endl;
    double x = pos.x();
    double vx = vel.x();
    if(x > (L/2))
    {
        x = L - x;
        vx *= -1;
    }

    if(x < - (L/2))
    {
        x = - (L/2) - (x + (L/2));
        vx *= -1;
    }

    double y = pos.y();
    double vy = vel.y();
    if(y > (L/2))
    {
        y = L - y;
        vy *= -1;
    }

    if(y < - (L/2))
    {
        y = - (L/2) - (y + (L/2));
        vy *= -1;
    }

    double z = pos.z();
    double vz = vel.z();
    if(z > (L/2))
    {
        z = L - z;
        vz *= -1;
    }

    if(z < - (L/2))
    {
        z = - (L/2) - (z + (L/2));
        vz *= -1;
    }

    vec r (x,y,z);
    pos  = r;

    vec v (vx,vy,vz);
    vel = v;
}



// -----------------------------------------------------------
void read_init(std::string input_file, std::vector<body> &system)
{
    std::string line; // Declare a string to store each line
    std::string name; // String to store body name
    double m, x, y, z, vx, vy, vz; // Doubles to store vector components
    int line_cnt = 0; // Line counter

    // Declare and initialise an input file stream object
    std::ifstream data_file(input_file);

    while (getline(data_file, line)) // Read the file line by line
    {
        line_cnt++;
        std::stringstream data_line(line); // Create a string stream from the line
        switch (line_cnt)
        {
            case 1:
                name = line;
                break;
            case 2:
                m = std::stod(line); // Convert string line into double
                break;
            case 3:
                read_vector3d(data_line, x, y, z); // Read the 3 components of the vector on the line
                break;
            case 4:
                read_vector3d(data_line, vx, vy, vz); // Read the 3 components of the vector on the line
                break;
        }
        if (line_cnt==4) // Data for one body has been extracted
        {
            line_cnt = 0; // Reset line counter
            body b(name,m,vec(x,y,z),vec(vx,vy,vz)); // Package data into body
            system.push_back(b); // Add body to system
        }
    }
    // Close the file
    data_file.close();
}

void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z)
{
    std::string value; // Declare a string to store values in a line
    int val_cnt = 0; // Value counter along the line
    while (getline(data_line, value, ','))
    {
        val_cnt++;
        switch (val_cnt)
        {
            case 1:
                x = std::stod(value);
                break;
            case 2:
                y = std::stod(value);
                break;
            case 3:
                z = std::stod(value);
                break;
        }
    }
}

void save_data(std::ofstream& savefile, const std::vector<body> &system, double t)
{
    // Function for saving the simulation data to file.

    vec L; // Total angular momentum
    double E = 0.0; // Total energy
    for(int p = 0; p < system.size(); p++)
    {
        E += system[p].get_ke() + 0.5*system[p].get_gpe();
        L += system[p].get_L();
    }
    double Lmag = L.length(); // Magnitude of total angular momentum

    // Write a header for this time-step with time, total energy and total mag of L:
    savefile << t  << endl; //<< "," << E << "," << Lmag << endl;

    // Loop over the bodies:
    for(int p = 0; p < system.size(); p++)
    {
        // Output position and velocity for each body:
        savefile << system[p].get_pos() << "," << system[p].get_vel() << endl;
    }
}
