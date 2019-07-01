#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{


public:

    int lx, ly, ns, orbs, IterMax, RandomSeed;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    int Total_Particles;
    double mus,mus_Cluster, Fill, pi, mu_old;
    double J_Hund;
    double U_onsite;
    double Lambda_SOC;
    string Auxilliary_fields_Seed_file_name;
    string Tag_output_files_with;
    double dw_dos, eta_dos,  w_min, w_max;


    double Disorder_Strength, RandomDisorderSeed;
    bool PBC;
    double Boltzman_constant;
    bool PNICTIDES_HOPPING;
    bool Read_Auxilliary_fields;
    double S_moment_max, S_moment_min;
    double L_moment_max, L_moment_min;
    double Rho_max, Rho_min;

    string File_OPs_in, File_OPs_out;

    Matrix<double> t2g_hopping_NN;
    Mat_1_doub Crystal_Field;


    bool Cooling_;
    bool ED_;

    bool Metropolis_Algorithm;
    bool Use_Saddle_point_Approx_on_density;
    bool Use_constant_local_Rho_den;

    char Dflag;

    bool Saving_Microscopic_States;
    int No_Of_Microscopic_States;

    double temp_max, beta_min;
    double temp_min, beta_max;
    double d_Temp;

    int Last_n_sweeps_for_measurement;
    int Measurement_after_each_m_sweeps;

    double Temperature, beta, Eav, maxmoment;
    double WindowSize, AccCount[2];

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){

    string PBC_string, Pnictides_Hopping_string;
    double ED_double, cooling_double;
    Boltzman_constant=1.0;
    double  Use_constant_local_Rho_den_double,Use_Saddle_point_Approx_on_density_double, Read_Auxilliary_fields_double ;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;


    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));
    TBC_mx = int(matchstring(inputfile_,"TwistedBoundaryCond_mx"));
    TBC_my = int(matchstring(inputfile_,"TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_,"TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_,"TBC_cellsY"));

    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;
    orbs=3;

    Fill = matchstring(inputfile_, "Fill");
    cout << "TotalNumberOfParticles = " << ns * Fill * 2.0 * orbs << endl;
    Total_Particles=ns * Fill * 2.0 * orbs;

    RandomSeed = matchstring(inputfile_,"RandomSeed");
    RandomDisorderSeed = matchstring(inputfile_,"RandomDisorderSeed");
    Disorder_Strength = matchstring(inputfile_,"Disorder_Strength");
    J_Hund = matchstring(inputfile_,"J_HUND");
    U_onsite = matchstring(inputfile_,"U_Onsite");
    Lambda_SOC = matchstring(inputfile_, "Lambda_SOC");
    S_moment_max = matchstring(inputfile_,"S_moment_max");
    S_moment_min = matchstring(inputfile_,"S_moment_min");
    L_moment_max = matchstring(inputfile_,"L_moment_max");
    L_moment_min = matchstring(inputfile_,"L_moment_min");
    Rho_max = matchstring(inputfile_,"Rho_max");
    Rho_min = matchstring(inputfile_,"Rho_min");
    dw_dos = matchstring(inputfile_,"dw_dos");
    eta_dos = matchstring(inputfile_,"eta_dos");
    w_min = matchstring(inputfile_,"w_min");
    w_max = matchstring(inputfile_,"w_max");


    Dflag = 'N';


    int SavingMicroscopicStates_int;
    SavingMicroscopicStates_int = int(matchstring(inputfile_, "SavingMicroscopicStates"));

    assert(SavingMicroscopicStates_int == 1 ||
           SavingMicroscopicStates_int == 0);
    if (SavingMicroscopicStates_int == 1)
    {
        Saving_Microscopic_States = true;
    }
    else
    {
        Saving_Microscopic_States = false;
    }

    No_Of_Microscopic_States = int(matchstring(inputfile_, "NoOfMicroscopicStates"));


    IterMax = int(matchstring(inputfile_, "MaxMCsweeps"));
    Last_n_sweeps_for_measurement = int(matchstring(inputfile_, "Last_n_sweeps_for_measurement"));
    Measurement_after_each_m_sweeps = int(matchstring(inputfile_, "Measurement_after_each_m_sweeps"));

    cooling_double = double(matchstring(inputfile_, "Cooling"));
    if (cooling_double == 1.0)
    {
        Cooling_ = true;

        temp_min = double(matchstring(inputfile_, "Temperature_min"));
        temp_max = double(matchstring(inputfile_, "Temperature_max"));
        d_Temp = double(matchstring(inputfile_, "dTemperature"));
        beta_max = double(Boltzman_constant / temp_min);
        beta_min = double(Boltzman_constant / temp_max);
    }
    else if (cooling_double == 0.0)
    {
        Cooling_ = false;

        Temperature = double(matchstring(inputfile_, "Temperature")); // temperature in kelvin
        beta = double(Boltzman_constant / Temperature);                         //Beta which is (T*k_b)^-1

        temp_min = Temperature;
        temp_max = Temperature;
        d_Temp = 10.0; //arbitrary positive number, not used in "not cooling"
    }
    else
    {
        cout << "ERROR: Cooling can be only 1 (true) or 0 (false)" << endl;
        assert(cooling_double == 0.0);
    }


    ED_double = double(matchstring(inputfile_, "Perform_ED"));
    if (ED_double == 1.0)
    {
        ED_ = true;
        lx_cluster = lx;
        ly_cluster = ly;
    }
    else if (ED_double == 0.0)
    {
        ED_ = false;
        lx_cluster = int(matchstring(inputfile_,"Cluster_lx"));
        ly_cluster = int(matchstring(inputfile_,"Cluster_ly"));

    }
    else
    {
        cout << "ERROR: Perform_ED can be only 1 (true) or 0 (false)" << endl;
        assert(ED_double == 0.0);
    }

    Read_Auxilliary_fields_double=double(matchstring(inputfile_,"Read_Auxilliary_fields"));
    if(Read_Auxilliary_fields_double==1.0){
        Read_Auxilliary_fields = true;
    }
    else if (Read_Auxilliary_fields_double == 0.0)
    {
        Read_Auxilliary_fields = false;
    }
    else
    {
        cout << "ERROR: Read_Auxilliary_fields can be only 1 (true) or 0 (false)" << endl;
        assert(Read_Auxilliary_fields_double == 0.0);
    }

    Auxilliary_fields_Seed_file_name=matchstring2(inputfile_,"Auxilliary_fields_Seed_file_name");
    Tag_output_files_with=matchstring2(inputfile_,"Tag_output_files_with");

    Use_Saddle_point_Approx_on_density_double=double(matchstring(inputfile_,"Use_Saddle_point_Approx_on_density"));
    if(Use_Saddle_point_Approx_on_density_double==1.0){
        Use_Saddle_point_Approx_on_density = true;
    }
    else if (Use_Saddle_point_Approx_on_density_double == 0.0)
    {
        Use_Saddle_point_Approx_on_density = false;
    }
    else
    {
        cout << "ERROR: Use_Saddle_point_Approx_on_density can be only 1 (true) or 0 (false)" << endl;
        assert(Use_Saddle_point_Approx_on_density_double == 0.0);
    }



    Use_constant_local_Rho_den_double=double(matchstring(inputfile_,"Use_constant_local_Rho_den"));
    if(Use_constant_local_Rho_den_double==1.0){
        Use_constant_local_Rho_den = true;
    }
    else if (Use_constant_local_Rho_den_double == 0.0)
    {
        Use_constant_local_Rho_den = false;
    }
    else
    {
        cout << "ERROR: Use_constant_local_Rho_den can be only 1 (true) or 0 (false)" << endl;
        assert(Use_constant_local_Rho_den_double == 0.0);
    }

    if( (Use_constant_local_Rho_den == Use_Saddle_point_Approx_on_density)
            &&
            Use_constant_local_Rho_den==true  ){
        cout <<"Use_constant_local_Rho_den and Use_Saddle_point_Approx_on_density cannot be true"<<endl;
        assert(!Use_constant_local_Rho_den);
    }


    PBC_string=matchstring2(inputfile_,"PBC");
    if(PBC_string=="true"){
        PBC=true;
    }
    else{
        PBC=false;
    }

    Pnictides_Hopping_string=matchstring2(inputfile_,"Pnictides_Hopping");
    if(Pnictides_Hopping_string=="true"){
        PNICTIDES_HOPPING=true;
    }
    else{
        PNICTIDES_HOPPING=false;
    }

    string Nearest_Neigh_Hopping_t2g_basis_row0;
    string Nearest_Neigh_Hopping_t2g_basis_row1;
    string Nearest_Neigh_Hopping_t2g_basis_row2;

    Nearest_Neigh_Hopping_t2g_basis_row0=matchstring2(inputfile_, "Nearest_Neigh_Hopping_t2g_basis_row0");
    Nearest_Neigh_Hopping_t2g_basis_row1=matchstring2(inputfile_, "Nearest_Neigh_Hopping_t2g_basis_row1");
    Nearest_Neigh_Hopping_t2g_basis_row2=matchstring2(inputfile_, "Nearest_Neigh_Hopping_t2g_basis_row2");

    stringstream t2g_row0_stream(Nearest_Neigh_Hopping_t2g_basis_row0);
    stringstream t2g_row1_stream(Nearest_Neigh_Hopping_t2g_basis_row1);
    stringstream t2g_row2_stream(Nearest_Neigh_Hopping_t2g_basis_row2);

    t2g_hopping_NN.resize(3,3);
    for(int n=0;n<3;n++){
        t2g_row0_stream >> t2g_hopping_NN(0,n);
        t2g_row1_stream >> t2g_hopping_NN(1,n);
        t2g_row2_stream >> t2g_hopping_NN(2,n);
    }


    string Crystal_Field_t2g;
    Crystal_Field_t2g=matchstring2(inputfile_, "Crystal_Field_t2g");

    stringstream Crystal_Field_t2g_stream(Crystal_Field_t2g);


    Crystal_Field.resize(3);
    for(int n=0;n<3;n++){
        Crystal_Field_t2g_stream >> Crystal_Field[n];
    }

    pi=4.00*atan(double(1.0));
    Eav=0.0;

    mus=0.25;
    mu_old=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



