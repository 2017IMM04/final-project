//This program simulates independent excitons and charges, with no potential between charges.
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <random>


using std::vector;
using std::cout;
using std::endl;
using std::cin;


class Particle{
public:
    int type, x, y, z;
    Particle(){};
    Particle(int a, int b, int c, int d){
        type = a; x = b; y = c; z = d;
    }
    Particle(const Particle& other){
        type = other.type; x = other.x; y = other.y; z = other.z;
    }
    Particle& operator=(const Particle& other){
        type = other.type; x = other.x; y = other.y; z = other.z;
        return *this;
    }
};

class Event{
public:
    int type, par, x, y, z, fx, fy, fz;
    Event(){};
    Event(int a, int b, int c, int d, int e, int f, int g, int h){
        type = a; par = b; x = c; y = d; z = e; fx = f; fy = g; fz = h;
    }
    Event(const Event& other){
        type = other.type; par = other.par; x = other.x; y = other.y; z = other.z; fx = other.fx; fy = other.fy; fz = other.fz;
    }
    Event& operator=(const Event& other){
        type = other.type; par = other.par; x = other.x; y = other.y; z = other.z; fx = other.fx; fy = other.fy; fz = other.fz;
        return *this;
    }
};


//function declarations
inline int particle_at(vector<int>&, int, int, int);
inline double potential_at(vector<double>&, int, int, int);
inline int morphology_at(vector<int>&, int, int, int);
inline void update(vector<int>&, vector<double>&, int, int, int, int);

// spatial dimensions, global variables
int N_x = 60;
int N_y = 60;
int N_z = 30;
int max_steps = 100000;
int insert_charges = 60;

//global variables, constants
const double work_function = 0;//volts
const double dielectric = 3.5*8.8542e-12;
const double pi = 3.14159265359;
const double elementary_charge = 1.60217e-19;
const double boltzmann_constant = 1.38064852e-23;
const double T = 298;
const double lattice_constant = 3e-9;//meters
const int z_floor = 5;
const int z_ceiling = 35;
const double coulombs_law = elementary_charge/(4*pi*dielectric*lattice_constant);
const double beta_factor = elementary_charge/boltzmann_constant/T;

//data structures
vector<Particle> exciton_array;
vector<Particle> electron_array;
vector<Particle> hole_array;
vector<int> deletion_index;
vector<int> morphology_matrix;
vector<int> particle_type_matrix;
vector<double> potential_matrix;

//rate constants
double exciton_dissociation = 1e40;
double exciton_creation = 100 * N_x * N_y * 9;
double exciton_hop = 1e13;//value paper gave seems in error
double exciton_recombination = 2e9;
    

double charge_hop = 280e12;
double charge_hop_exponent = 2;
double charge_recombine = 1e7;
    
    int recombined = 0;
    int ex_combine = 0;
    int ex_hop = 0;
    int e_hop = 0;
    int h_hop = 0;
    
int excitons_created = 100000;

//main loop
int main(){
//SETUP
    //intialize random variables
    std::mt19937 generator(time(0));
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);
    std::uniform_int_distribution<int> x_y_distribution(0, N_x-1);
    //std::uniform_int_distribution<int> z_distribution(0, N_z-1);
    std::normal_distribution<double> z_distribution(15.0, 8.0);
    
    //linear vector to store data for fixed size 3d matrices; access by [x+y*N_x+z*N_x*N_y]
    morphology_matrix.resize(N_x*N_y*N_z, 0);
    particle_type_matrix.resize(N_x*N_y*N_z, 0);
    potential_matrix.resize(N_x*N_y*N_z, 0);
    int temp;
    int index = 0;
    //input morhpology
    while(cin >> temp){
        morphology_matrix[index]=temp;
        index ++;
    }
    //intialize potential
    for(int x = 0; x < N_x; x++){
        for(int y = 0; y < N_x; y++){
            for(int z = 0; z < N_z; z++){
                potential_matrix[x+y*N_x+z*N_x*N_y]= (z-(1.0+N_z)/2)/N_z*work_function;
            }
        }
    }//potential initialized to linear potential, with z=0 negative
        
        
    //variables to store global results
    int exciton_dissociations = 0;
    int exciton_recombinations = 0;
    int electrons_extracted = 0;
    int holes_extracted = 0;
    int recombinations = 0;
    double sum_exciton_lifetime = 0;
    double sum_exciton_path_length = 0;
    double sum_electron_lifetime = 0;
    double sum_electron_path_length = 0;
    double sum_hole_lifetime = 0;
    double sum_hole_path_length = 0;
    
        
    for(int n = 0; n < excitons_created; n++){//average results over n excitons
        
        //bookkeeping variables
        double smallest_tau = 1;
        double temp_tau = 0;
        Event smallest_event;
        int x,y,z;
        int ex,ey,ez;//electron coordinates
        int hx,hy,hz;//hole coordinates
        
        //exciton creation
        double exciton_lifetime = 0;
        double exciton_path_length = 0;
        
        x = x_y_distribution(generator);
        y = x_y_distribution(generator);
        z = floor(z_distribution(generator));
        if(z >= 30){z = 29;}
        if(z < 0){z = 0;}
        bool exciton_alive = 1;
        bool dissociate = 0;//default no dissociation
        while(exciton_alive){//while exciton still exists, do stuff
            int morph = morphology_at(morphology_matrix,x,y,z);
            smallest_tau = 1;//initialize smallest_tau
            for(int i = -3; i<4; i++){//for neighbors
                for(int j = -3; j<4; j++){
                    for(int k = -3; k<4; k++){
                        int r_square = i*i + j*j + k*k;
                        if(r_square <= 9 && r_square != 0){//if not self and is within cutoff
                            if(z+k >= 0 && z+k < N_z){//if new site not out of bounds of electrodes
                                if(r_square == 1){//if adjacent, consider dissociation
                                    if(morphology_at(morphology_matrix,x,y,z)*morphology_at(morphology_matrix,x+i,y+j,z+k) == -1){//at interface
                                        temp_tau = -1/exciton_dissociation*log(real_distribution(generator));
                                        if(temp_tau < smallest_tau){
                                            smallest_tau = temp_tau;
                                            if(morph == 1){
                                                smallest_event = Event(1,1,x,y,z,x+i,y+j,z+k);
                                            }
                                            else{
                                                smallest_event = Event(1,-1,x,y,z,x+i,y+j,z+k);
                                            }
                                        }
                                    }
                                }
                                //exciton hop
                                temp_tau = -1/(exciton_hop/(pow(r_square,3.0))/729)*log(real_distribution(generator));
                                if(temp_tau < smallest_tau){
                                    smallest_tau = temp_tau;
                                    smallest_event = Event(2,2,x,y,z,x+i,y+j,z+k);
                                }
                                //exciton recombination
                                temp_tau = -1/exciton_recombination*log(real_distribution(generator));
                                if(temp_tau < smallest_tau){
                                    smallest_tau = temp_tau;
                                    smallest_event = Event(3,2,x,y,z,0,0,0);
                                }
                            }
                        }
                    }
                }
            }
            switch(smallest_event.type){
                case 1://exciton dissociation
                    if(smallest_event.par == 1){
                        hx = smallest_event.x;hy = smallest_event.y; hz = smallest_event.z;
                        ex = smallest_event.fx; ey = smallest_event.fy; ez = smallest_event.fz;
                    }
                    else{
                        hx = smallest_event.fx;hy = smallest_event.fy; hz = smallest_event.fz;
                        ex = smallest_event.x; ey = smallest_event.y; ez = smallest_event.z;
                    }
                    exciton_lifetime += smallest_tau;
                    sum_exciton_lifetime += exciton_lifetime;
                    sum_exciton_path_length += exciton_path_length;
                    exciton_alive = 0;
                    dissociate = 1;
                    exciton_dissociations++;
                    break;
                case 2://exciton hop
                    exciton_path_length += sqrt(pow(x-smallest_event.fx,2)+pow(y-smallest_event.fy,2)+pow(z-smallest_event.fz,2));
                    exciton_lifetime += smallest_tau;
                    x = smallest_event.fx; y = smallest_event.fy; z = smallest_event.fz;
                    break;
                case 3://exciton recombination
                    exciton_lifetime += smallest_tau;
                    sum_exciton_lifetime += exciton_lifetime;
                    sum_exciton_path_length += exciton_path_length;
                    exciton_alive = 0;
                    exciton_recombinations++;
                    break;
            }
        }
        
        if(dissociate){
            double electron_lifetime = 0;
            double electron_path_length = 0;
            double hole_lifetime = 0;
            double hole_path_length = 0;
            int hops = 0;
            bool electron_alive = 1;
            bool hole_alive = 1;
            while(hops < 5000){
                smallest_tau = 1;
                if(ez == N_z-1 && electron_alive){//electron reached upper electrode
                    sum_electron_lifetime += electron_lifetime;
                    sum_electron_path_length += electron_path_length;
                    electrons_extracted++;
                    electron_alive = 0;
                }
                if(hz == 0 && hole_alive){//hole reached lower electrode
                    sum_hole_lifetime += hole_lifetime;
                    sum_hole_path_length += hole_path_length;
                    holes_extracted++;
                    hole_alive = 0;
                }
                
                if(!electron_alive && !hole_alive){
                    break;//both have been extracted or recombined
                }
                if(electron_alive){
                    for(int i = -3; i<4; i++){
                        for(int j = -3; j<4; j++){
                            for(int k = -3; k<4; k++){
                                int r_square = i*i + j*j + k*k;
                                if(r_square <= 9 && r_square != 0){//if not self and is within cutoff
                                    if(ez+k >= 0 && ez+k < N_z){//if not out of bounds of electrodes
                                        if(!(ex+i == hx && ey+j == hy && ez+k == hz) || !hole_alive){//that site is not the hole
                                            //electron hop
                                            if(morphology_at(morphology_matrix,ex+i,ey+j,ez+k) == -1){//material is electron conductor
                                                double energy_diff = -(potential_at(potential_matrix,ex+i,ey+j,ez+k)-potential_at(potential_matrix,ex,ey,ez));
                                                double factor = 0.0;
                                                if(energy_diff < 0){
                                                    factor = 1;
                                                }
                                                else{
                                                    factor = exp(-energy_diff*beta_factor);
                                                }
                                                temp_tau = -1/(charge_hop*factor*exp(-12*sqrt(r_square)))*log(real_distribution(generator));
                                                if(temp_tau < smallest_tau){
                                                    smallest_event = Event(4, -1, ex,ey,ez,ex+i,ey+j,ez+k);
                                                    smallest_tau = temp_tau;
                                                }
                                            }
                                            //electron recombination
                                            if(morphology_at(morphology_matrix, ex+i,ey+j,ez+k) == 1 && r_square == 1){
                                                temp_tau = -1/charge_recombine*log(real_distribution(generator));
                                                if(temp_tau < smallest_tau){
                                                    smallest_event = Event(6, -1, hx,hy,hz,hx+i,hy+j,hz+k);
                                                    smallest_tau = temp_tau;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(hole_alive){
                    for(int i = -3; i<4; i++){
                        for(int j = -3; j<4; j++){
                            for(int k = -3; k<4; k++){
                                int r_square = i*i + j*j + k*k;
                                if(r_square <= 9 && r_square != 0){//if not self and is within cutoff
                                    if(hz+k >= 0 && z+k < N_z){//if not out of bounds of electrodes
                                        if(!(hx+i == ex && hy+j == ey && hz+k == ez) || !electron_alive){//not electron there
                                            //hole hop
                                            if(morphology_at(morphology_matrix,hx+i,hy+j,hz+k) == 1){//material is hole conductor
                                                double energy_diff = (potential_at(potential_matrix,hx+i,hy+j,hz+k)-potential_at(potential_matrix,hx,hy,hz));
                                                double factor = 0.0;
                                                if(energy_diff < 0){
                                                    factor = 1;
                                                }
                                                else{
                                                    factor = exp(-energy_diff*beta_factor);
                                                }
                                                temp_tau = -1/(charge_hop*factor*exp(-12*sqrt(r_square)))*log(real_distribution(generator));
                                                if(temp_tau < smallest_tau){
                                                    smallest_event = Event(5, -1, hx,hy,hz,hx+i,hy+j,hz+k);
                                                    smallest_tau = temp_tau;
                                                }
                                            }
                                            //hole recombination
                                            if(morphology_at(morphology_matrix, hx+i,hy+j,hz+k) == -1 && r_square == 1){
                                                temp_tau = -1/charge_recombine*log(real_distribution(generator));
                                                if(temp_tau < smallest_tau){
                                                    smallest_event = Event(7, 1, hx,hy,hz,hx+i,hy+j,hz+k);
                                                    smallest_tau = temp_tau;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(electron_alive || hole_alive){//there will be events
                    switch(smallest_event.type){
                        case 4://electron hop
                            electron_path_length += sqrt(pow(ex-smallest_event.fx,2)+pow(ey-smallest_event.fy,2)+pow(ez-smallest_event.fz,2));
                            electron_lifetime += smallest_tau;
                            if(hole_alive){
                                hole_lifetime += smallest_tau;
                            }
                            ex = smallest_event.fx; ey = smallest_event.fy; ez = smallest_event.fz;
                            break;
                        case 5://hole hop
                            hole_path_length += sqrt(pow(hx-smallest_event.fx,2)+pow(hy-smallest_event.fy,2)+pow(hz-smallest_event.fz,2));
                            hole_lifetime += smallest_tau;
                            if(electron_alive){
                                electron_lifetime += smallest_tau;
                            }
                            hx = smallest_event.fx; hy = smallest_event.fy; hz = smallest_event.fz;
                            break;
                        case 6://electron recombination
                            electron_lifetime += smallest_tau;
                            if(hole_alive){hole_lifetime += smallest_tau;}
                            sum_electron_lifetime += electron_lifetime;
                            sum_electron_path_length += electron_path_length;
                            electron_alive = 0;
                            recombinations++;
                            break;
                        case 7://hole recombination
                            hole_lifetime += smallest_tau;
                            if(electron_alive){electron_lifetime += smallest_tau;}
                            sum_hole_lifetime += hole_lifetime;
                            sum_hole_path_length += hole_path_length;
                            hole_alive = 0;
                            recombinations++;
                            break;
                    }
                }
                hops++;
            }
        }
    }
    
    //OUTPUT
    cout << excitons_created << " excitons created" << endl;
    cout << exciton_dissociations << " excitons dissociated" << endl;
    cout << exciton_recombinations << " excitons recombined" << endl;
    cout << electrons_extracted << " electrons extracted" << endl;
    cout << holes_extracted << " holes extracted" << endl;
    cout << recombinations << " charges recombined" << endl;
    cout << sum_exciton_lifetime/excitons_created << " mean exciton lifetime" << endl;
    //cout << sum_electron_lifetime/(electrons_extracted+recombinations) << " mean electron lifetime" << endl;
    //cout << sum_hole_lifetime/(holes_extracted+recombinations) << " mean hole lifetime" << endl;
    cout << sum_exciton_path_length << " sum of exciton hop length" << endl;
    //cout << sum_electron_path_length << endl;
    //cout << sum_hole_path_length << endl;
    return 0;
}


//functions
inline int particle_at(vector<int>& particle_type_matrix, int x, int y, int z){
    return particle_type_matrix[(x%N_x+N_x)%N_x+(y%N_y+N_y)%N_y*N_x+z%N_z*N_x*N_y];
}

inline double potential_at(vector<double>& potential_matrix, int x, int y, int z){
    return potential_matrix[(x%N_x+N_x)%N_x+(y%N_y+N_y)%N_y*N_x+z%N_z*N_x*N_y];
}

inline int morphology_at(vector<int>& morphology_matrix, int x, int y, int z){
    return morphology_matrix[(x%N_x+N_x)%N_x+(y%N_y+N_y)%N_y*N_x+z%N_z*N_x*N_y];
}

inline void update(vector<int>& particle_type_matrix, vector<double>& potential_matrix, int type, int x, int y, int z){
    x = (x%N_x+N_x)%N_x; y = (y%N_y+N_y)%N_y; z = z%N_z;
    particle_type_matrix[x+y*N_x+z*N_x*N_y] = type;
    return;
}
