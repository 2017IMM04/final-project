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

class Bookkeep{
public:
    double smallest_tau;
    int exciton_index, electron_index, hole_index;
    int smallest_exciton, smallest_electron, smallest_hole;
    Event smallest_event;
    Bookkeep(){
        exciton_index = 0;
        electron_index = 0;
        hole_index = 0;
        smallest_exciton = 0;
        smallest_hole = 0;
        smallest_electron = 0;
        smallest_tau = 0;
        smallest_event = Event(0,0,0,0,0,0,0,0);
    }
};

//function declarations
inline int particle_at(vector<int>&, int, int, int);
inline double potential_at(vector<double>&, int, int, int);
inline int morphology_at(vector<int>&, int, int, int);
inline void update(vector<int>&, vector<double>&, int, int, int, int);

// spatial dimensions, global variables
int actual_x = 60;
int actual_y = 60;
int actual_z = 30;
int N_x = 60;
int N_y = 60;
int N_z = 40;
int max_steps = 80000000;

//global variables, constants
const double work_function = 0.5;//volts
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

int main(){
    //intialize random variables
    std::mt19937 generator(time(0));
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);
    std::uniform_int_distribution<int> x_y_distribution(0, N_x-1);
    std::uniform_int_distribution<int> z_distribution(0, actual_z-1);
    
    //rate constants
    double exciton_dissociation = 1e40;
    double exciton_creation = 100 * N_x * N_y * 9;
    double exciton_hop = 0.3e12;//value paper gave seems in error
    double exciton_recombination = 2e9;
    
    double charge_hop = 280e12;
    double charge_hop_exponent = 2;
    double charge_recombine = 1e10;
    
    double clock = 0;//clock time initialization
    double time = 1; //run for time
    
    double exciton_clock = 0;
    
    
    int extracted_electrons = 0;
    int extracted_holes = 0;
    int excitons_created = 0;
    int hops = 0;
    int recombined = 0;
    int ex_combine = 0;
    int ex_hop = 0;
    int e_hop = 0;
    int h_hop = 0;
    
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
    
    for(int x = 0; x < N_x; x++){
        for(int y = 0; y < N_x; y++){
            for(int z = 0; z < N_z; z++){
                potential_matrix[x+y*N_x+z*N_x*N_y]= (z-(1.0+N_z)/2)/actual_z*work_function;
            }
        }
    }//potential initialized to linear potential, with z=0 negative
    
    Bookkeep bookkeeper;
    
    int steps = 0;
    //main loop
    int flag =0;
    while(steps < max_steps){
        //infinite rate process
        //charge extraction
        
        deletion_index.clear();//vector to keep track of index of electrons that need to be deleted
        for(int i = electron_array.size()-1; i>=0;i--){
            int temp_x = electron_array[i].x;
            int temp_y = electron_array[i].y;
            int temp_z = electron_array[i].z;
            if(temp_z == z_ceiling-1){//reached electrode
                update(particle_type_matrix,potential_matrix,0,temp_x,temp_y,temp_z);
                electron_array.erase(electron_array.begin()+i);
                extracted_electrons ++;
            }
        }

        for(int i = hole_array.size()-1; i>=0;i--){
            int temp_x = hole_array[i].x;
            int temp_y = hole_array[i].y;
            int temp_z = hole_array[i].z;
            if(temp_z == z_floor){//reached electrode
                update(particle_type_matrix,potential_matrix,0,temp_x,temp_y,temp_z);
                hole_array.erase(hole_array.begin()+i);
                extracted_holes ++;
            }
        }
        
        double temp_tau = -1/exciton_creation*log(real_distribution(generator));
        bookkeeper.smallest_tau = temp_tau;
        int temp_x, temp_y, temp_z;
        do{
            temp_x = x_y_distribution(generator);
            temp_y = x_y_distribution(generator);
            temp_z = z_distribution(generator)+z_floor;
        }while(particle_type_matrix[temp_x+temp_y*N_x+temp_z*N_x*N_y] != 0);//if occupied, generate another
        bookkeeper.smallest_event = Event(0,2,temp_x,temp_y,temp_z,0,0,0);//set to exciton creation
        
        //exciton processes
        for(int n = 0; n < exciton_array.size(); n++){
            int x = exciton_array[n].x;
            int y = exciton_array[n].y;
            int z = exciton_array[n].z;
            int morph = morphology_at(morphology_matrix,x,y,z);
            double temp_tau = 0;
            for(int i = -3; i<4; i++){
                for(int j = -3; j<4; j++){
                    for(int k = -3; k<4; k++){
                        int r_square = i*i + j*j + k*k;
                        if(r_square <= 9 && r_square != 0){//if not self and is within cutoff
                            if(z+k >= z_floor && z+k < z_ceiling){//if new site not out of bounds of electrodes
                                if(particle_at(particle_type_matrix,x+i,y+j,z+k) == 0){//if empty site, consider dissociation and hop
                                    if(r_square == 1){//if adjacent, consider dissociation
                                        if(morphology_at(morphology_matrix,x,y,z)*morphology_at(morphology_matrix,x+i,y+j,z+k) == -1){//at interface
                                            temp_tau = -1/exciton_dissociation*log(real_distribution(generator));
                                            if(temp_tau < bookkeeper.smallest_tau){
                                                bookkeeper.smallest_tau = temp_tau;
                                                bookkeeper.smallest_exciton = n;
                                                if(morph == 1){
                                                    bookkeeper.smallest_event = Event(1,1,x,y,z,x+i,y+j,z+k);
                                                }
                                                else{
                                                    bookkeeper.smallest_event = Event(1,-1,x,y,z,x+i,y+j,z+k);
                                                }
                                            }
                                        }
                                    }
                                    //exciton hop
                                    temp_tau = -1/(exciton_hop/(pow(r_square,3.0))/729)*log(real_distribution(generator));
                                    if(temp_tau < bookkeeper.smallest_tau){
                                        bookkeeper.smallest_tau = temp_tau;
                                        bookkeeper.smallest_event = Event(2,2,x,y,z,x+i,y+j,z+k);
                                        bookkeeper.smallest_exciton = n;
                                    }
                                }
                                //exciton recombination
                                temp_tau = -1/exciton_recombination*log(real_distribution(generator));
                                if(temp_tau < bookkeeper.smallest_tau){
                                    bookkeeper.smallest_tau = temp_tau;
                                    bookkeeper.smallest_event = Event(3,2,x,y,z,0,0,0);
                                    bookkeeper.smallest_exciton = n;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //electrons processes
        for(int n = 0; n < electron_array.size(); n++){
            int x = electron_array[n].x;
            int y = electron_array[n].y;
            int z = electron_array[n].z;
            double temp_tau = 0;
            for(int i = -3; i<4; i++){
                for(int j = -3; j<4; j++){
                    for(int k = -3; k<4; k++){
                        int r_square = i*i + j*j + k*k;
                        if(r_square <= 9 && r_square != 0){//if not self and is within cutoff
                            if(z+k >= z_floor && z+k < z_ceiling){//if not out of bounds of electrodes
                                if(particle_at(particle_type_matrix,x+i,y+j,z+k) == 0){
                                    //electron hop
                                    if(morphology_at(morphology_matrix,x+i,y+j,z+k) == -1){//material is electron conductor
                                        double energy_diff = -(potential_at(potential_matrix,x+i,y+j,z+k)-potential_at(potential_matrix,x,y,z));
                                        double factor = 0.0;
                                        if(energy_diff < 0){
                                            factor = 1;
                                        }
                                        else{
                                            factor = exp(-energy_diff*beta_factor);
                                        }
                                        temp_tau = -1/(charge_hop*factor*exp(-12*sqrt(r_square)))*log(real_distribution(generator));
                                        if(temp_tau < bookkeeper.smallest_tau){
                                            bookkeeper.smallest_event = Event(4, -1, x,y,z,x+i,y+j,z+k);
                                            bookkeeper.smallest_electron = n;
                                            bookkeeper.smallest_tau = temp_tau;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //hole processes
        for(int n = 0; n < hole_array.size(); n++){
            int x = hole_array[n].x;
            int y = hole_array[n].y;
            int z = hole_array[n].z;
            double temp_tau = 0;
            for(int i = -3; i<4; i++){
                for(int j = -3; j<4; j++){
                    for(int k = -3; k<4; k++){
                        int r_square = i*i + j*j + k*k;
                        if(r_square <= 9 && r_square != 0){//if not self and is within cutoff
                            if(z+k >= z_floor && z+k < z_ceiling){//if not out of bounds of electrodes
                                if(particle_at(particle_type_matrix,x+i,y+j,z+k) == 0){
                                    //hole hop
                                    if(morphology_at(morphology_matrix,x+i,y+j,z+k) == 1){//material is hole conductor
                                        double energy_diff = (potential_at(potential_matrix,x+i,y+j,z+k)-potential_at(potential_matrix,x,y,z));
                                        double factor = 0.0;
                                        if(energy_diff < 0){
                                            factor = 1;
                                        }
                                        else{
                                            factor = exp(-energy_diff*beta_factor);
                                        }
                                        temp_tau = -1/(charge_hop*factor*exp(-12*sqrt(r_square)))*log(real_distribution(generator));
                                        if(temp_tau < bookkeeper.smallest_tau){
                                            bookkeeper.smallest_event = Event(5, -1, x,y,z,x+i,y+j,z+k);
                                            bookkeeper.smallest_hole = n;
                                            bookkeeper.smallest_tau = temp_tau;
                                        }
                                    }
                                }
                                //recombination
                                if(r_square == 1){//adjacent
                                    if(particle_at(particle_type_matrix,x+i,y+j,z+k) == -1){//if an electron
                                        temp_tau = -1/charge_recombine*log(real_distribution(generator));
                                        if(temp_tau < bookkeeper.smallest_tau){
                                            
                                            for(int m = 0; m < electron_array.size(); m++){
                                                if((electron_array[m].x-x-i)%N_x == 0 && (electron_array[m].y-y-j)%N_y == 0 && electron_array[m].z-z-k == 0){
                                                    bookkeeper.smallest_tau = temp_tau;
                                                    bookkeeper.smallest_event = Event(6, 1, x,y,z,x+i,y+j,z+k);
                                                    bookkeeper.smallest_hole = n;
                                                    bookkeeper.smallest_electron = m;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //run event with the smallest waiting time
        int x = bookkeeper.smallest_event.x;
        int y = bookkeeper.smallest_event.y;
        int z = bookkeeper.smallest_event.z;
        int fx = bookkeeper.smallest_event.fx;
        int fy = bookkeeper.smallest_event.fy;
        int fz = bookkeeper.smallest_event.fz;
        if(bookkeeper.smallest_event.type == 0){//exciton creation
            exciton_array.push_back(Particle(2, x,y,z));
            update(particle_type_matrix,potential_matrix,2,x,y,z);
            excitons_created ++;
        }
        else if(bookkeeper.smallest_event.type == 1){//exciton dissociation
            if(bookkeeper.smallest_event.par == 1){
                update(particle_type_matrix,potential_matrix, 1,x,y,z);
                update(particle_type_matrix,potential_matrix, -1,fx,fy,fz);
                hole_array.push_back(Particle(1,x,y,z));
                electron_array.push_back(Particle(-1,fx,fy,fz));
            }
            else{
                update(particle_type_matrix,potential_matrix, -1,x,y,z);
                update(particle_type_matrix,potential_matrix, 1,fx,fy,fz);
                electron_array.push_back(Particle(-1,x,y,z));
                hole_array.push_back(Particle(1,fx,fy,fz));
            }
            exciton_array.erase(exciton_array.begin()+bookkeeper.smallest_exciton);
        }
        else if(bookkeeper.smallest_event.type == 2){//exciton hop
            update(particle_type_matrix,potential_matrix, 0,x,y,z);
            update(particle_type_matrix,potential_matrix,2,fx,fy,fz);
            exciton_array.erase(exciton_array.begin()+bookkeeper.smallest_exciton);
            exciton_array.push_back(Particle(2,fx,fy,fz));
            ex_hop++;
        }
        else if(bookkeeper.smallest_event.type == 3){//exciton recombination
            update(particle_type_matrix,potential_matrix, 0,x,y,z);
            exciton_array.erase(exciton_array.begin()+bookkeeper.smallest_exciton);
            ex_combine++;
        }
        else if(bookkeeper.smallest_event.type == 4){//electron hop
            update(particle_type_matrix,potential_matrix, 0,x,y,z);
            update(particle_type_matrix,potential_matrix,-1,fx,fy,fz);
            electron_array.erase(electron_array.begin()+bookkeeper.smallest_electron);
            electron_array.push_back(Particle(-1,fx,fy,fz));
            e_hop++;
        }
        else if(bookkeeper.smallest_event.type == 5){//hole hop
            update(particle_type_matrix,potential_matrix, 0,x,y,z);
            update(particle_type_matrix,potential_matrix,1,fx,fy,fz);
            hole_array.erase(hole_array.begin()+bookkeeper.smallest_hole);
            hole_array.push_back(Particle(1,fx,fy,fz));
            h_hop++;
        }
        else if(bookkeeper.smallest_event.type == 6){//charge recombination
            update(particle_type_matrix,potential_matrix, 0,x,y,z);
            update(particle_type_matrix,potential_matrix, 0,fx,fy,fz);
            electron_array.erase(electron_array.begin()+bookkeeper.smallest_electron);
            hole_array.erase(hole_array.begin()+bookkeeper.smallest_hole);
            recombined ++;
        }
        steps ++;
        clock += bookkeeper.smallest_tau;
        
        if(steps%10000 == 0){
            cout << clock << "\t";
            cout << excitons_created << "\t";
            cout << extracted_holes << "\t";
            cout << extracted_electrons << "\t";
            cout << exciton_array.size() << "\t";
            cout << hole_array.size() << "\t";
            cout << electron_array.size() << "\t";
            cout << recombined << "\t";
            cout << ex_combine << "\n";
            //cout << ex_hop << "\t";
            //cout <<  e_hop << "\t";
            //cout << h_hop << "\n";
            }
    }
    
    /*for(int i = 0; i < exciton_array.size(); i++){
        cout << exciton_array[i].x << " " << exciton_array[i].y << " " << exciton_array[i].z << endl;
    }
    for(int i = 0; i < electron_array.size(); i++){
        cout << electron_array[i].x << " " << electron_array[i].y << " " << electron_array[i].z << endl;
    }
    for(int i = 0; i < exciton_array.size(); i++){
        cout << hole_array[i].x << " " << hole_array[i].y << " " << hole_array[i].z << endl;
    }*/
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
    //destroy potential caused by the original particle
    double final = 0.0;
    int par_type = particle_at(particle_type_matrix, x,y,z);//get original particle there
    if(abs(par_type) == 1){// if a charge
        for(int i = -3; i<4; i++){
            for(int j = -3; j<4; j++){
                for(int k = -3; k<4; k++){
                    double original = potential_at(potential_matrix,x+i,y+j,z+k);
                    int r_square = i*i + j*j + k*k;
                    if(r_square == 0){
                        final = original - par_type*coulombs_law;
                        potential_matrix[x+y*N_x+z*N_x*N_y] = final;
                    }
                    else if(r_square <= 9){
                        final = original - par_type*coulombs_law/sqrt(r_square);
                        potential_matrix[x+y*N_x+z*N_x*N_y] = final;
                    }
                }
            }
        }
    }
    //change particle
    particle_type_matrix[x+y*N_x+z*N_x*N_y] = type;
    if(abs(type) == 1){
        for(int i = -3; i<4; i++){
            for(int j = -3; j<4; j++){
                for(int k = -3; k<4; k++){
                    int r_square = i*i + j*j + k*k;
                    double original = potential_at(potential_matrix,x+i,y+j,z+k);
                    if(r_square == 0){
                        final = original + type*coulombs_law;
                        potential_matrix[x+y*N_x+z*N_x*N_y] = final;
                    }
                    else if(r_square <= 9){
                        final = original + type*coulombs_law/sqrt(r_square);
                        potential_matrix[x+y*N_x+z*N_x*N_y] = final;
                        
                    }
                }
            }
        }
    }
    return;
}
