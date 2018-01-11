#pragma warning(disable:4786)
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <random> // C++0x //

using namespace std;

uniform_real_distribution<double> runif(0.0,1.0);

const double min_step = 0.1;
const double plot_interval = 1.0;
const double t0 = 0.0;
const double Tmax = 100000.0;

struct CFG {
    vector<double> crvals;
    double srv;
    double dt;
    double r1;
    double r2;
    int rID;
};

void reactions( int Nr, vector<int> pnum, vector<double> pars, struct CFG& cfg ) {
    
    double r_rate[Nr];
    // parameters //
    double kon = pars[0];
    double koff = pars[1];
    double km = pars[2];
    double kp = pars[3];
    double kbind = pars[4];
    double kunbind = pars[5];
    double ka = pars[6];
    double dm = pars[7];
    double dp = pars[8];
    // variables //
    int LTRoff = pnum[0];
    int LTRon = pnum[1];
    int mRNA = pnum[2];
    int Tax = pnum[3];
    int LTRonTax = pnum[4];
    // reactions //
    r_rate[0] = kon*LTRoff;
    r_rate[1] = koff*LTRon;
    r_rate[2] = km*LTRon;
    r_rate[3] = kp*mRNA;
    r_rate[4] = kbind*LTRoff*Tax;
    r_rate[5] = kunbind*LTRonTax;
    r_rate[6] = ka*LTRonTax;
    r_rate[7] = dm*mRNA;
    r_rate[8] = dp*Tax;
    
    // routine //
    for( int i=0; i<Nr; i++ ) {
        cfg.srv += r_rate[i];
        cfg.crvals[i] = cfg.srv;
    }
}

void generate_random_numbers( mt19937& eng, struct CFG& cfg ) {
    cfg.r1 = runif(eng);
    cfg.r2 = runif(eng);
}


// /*
void update_pnum( vector<int>& nu, struct CFG& cfg ) {
    int rID = cfg.rID;
    if( rID==0 ) {
        // 1. Activating one LTR molecule from off to on #
        nu[0] += -1;
        nu[1] += 1;
    } else if( rID==1 ) {
        // 2. Deactivating one LTR molecule from on to off #
        nu[0] += 1;
        nu[1] += -1;
    } else if( rID==2 ) {
        /// 3. Transcription of RNA #
        nu[2] += 1;
    } else if( rID==3 ) {
        // 4. Translation of Tax #
        nu[3] += 1;
    } else if( rID==4 ) {
        // 5. Tax binding #
        nu[0] += -1;
        nu[3] += -1;
        nu[4] += 1;
    } else if( rID==5 ) {
        // 6. Tax unbinding #
        nu[0] += 1;
        nu[3] += 1;
        nu[4] += -1;
    } else if( rID==6 ) {
        // 7. activation #
        nu[2] += 1;
    } else if( rID==7 ) {
        // 8. mRNA decay #
        nu[2] += -1;
    } else if( rID==8 ) {
        // 9. Tax decay #
        nu[3] += -1;
    }
}
// */

void Gillespie_direct( int Nr, vector<int>& pnum, vector<double> pars, mt19937& eng, struct CFG& cfg ) {
    // generate two runif //
    generate_random_numbers(eng,cfg);
    cfg.srv = 0.0;
    // calculate reactions //
    reactions(Nr,pnum,pars,cfg);
    if( cfg.srv==0.0 ) {
        cfg.dt = min_step;
        cfg.rID = -1;
    } else {
        double val;
        double thr = cfg.r2*cfg.srv;
        cfg.rID = -1;
        int i = 1;
        val = cfg.crvals[0];
        if( val>=thr ) {
            cfg.rID = 0;
        } else {
            while( val<thr ) {
                val = cfg.crvals[i];
                i += 1;
            }
            cfg.rID = i-1;
        }
        cfg.dt = -log(cfg.r1)/cfg.srv;
    }
    // select the numbe of reactions to occur //
    update_pnum(pnum,cfg);
}

void shell_input( int Np, vector<double>& pars, int& rnd_seed, int argc, char **argv ) {
    
    int i;
    
    if( argc==Np+2 ) {
        rnd_seed = (int)atoi(argv[1]);
        for( i=0; i<Np; i++ ) {
            pars[i] = (double)atof(argv[i+2]);
        }
    } else if( argc==1 ) {
        cout << "Default parameter values are used." << endl;
        rnd_seed = 1;
        pars[0] = 0.0001; // kon
        pars[1] = 0.01; // koff
        pars[2] = 0.1; // km
        pars[3] = 10.0; // kp
        pars[4] = 0.001; // kbind
        pars[5] = 0.01; // kunbind
        pars[6] = 3.0; // ka
        pars[7] = 1.0; // dm
        pars[8] = 0.125; // dp
    } else {
        cout << "excess number of input (quit computation)." << endl;
        exit(1);
    }
}

int main( int argc, char** argv ) {
    
    int rnd_seed;
    //random_device generate_seed;
    //mt19937 eng(generate_seed());
    
    int Np = 9; // number of parameters
    vector<double> pars(Np);
    shell_input(Np,pars,rnd_seed,argc,argv);
    mt19937 eng(rnd_seed);

    string save_file_name = "Taxsim.txt";
    ofstream ofs( save_file_name );

    struct CFG cfg;
    int Nr = 9; // number of reactions
    for( int i=0; i<Nr; i++ ) {
         cfg.crvals.push_back(0.0);
    }
    double t = t0;
    double next_plot = plot_interval;
    
    int Ns = 5; /// number of species
    vector<int> pnum(Ns);
    pnum[0] = 0; // LTRoff
    pnum[1] = 1; // LTRon
    pnum[2] = 0; // mRNA
    pnum[3] = 0; // Tax
    pnum[4] = 0; // LTRonTax
    
    while( t<Tmax ) {
        Gillespie_direct(Nr,pnum,pars,eng,cfg);
        ofs << t << "\t" << pnum[3] << endl;
        t += cfg.dt;
    }
    return 0;

}
