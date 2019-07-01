#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;
#include "Matrix.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "MCEngine.h"
#include "random"


int main(int argc, char *argv[]) {
    if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }
    string ex_string_original =argv[0];
    cout<<"Executable used : "<<ex_string_original<<endl;
    string ex_string;
    ex_string=ex_string_original.substr (ex_string_original.length() - 3);

    string inputfile = argv[1];

    bool check_Non_Int=false;


    Parameters Parameters_;
    Parameters_.Initialize(inputfile);

    Coordinates Coordinates_(Parameters_.lx, Parameters_.ly);
    Coordinates CoordinatesCluster_(Parameters_.lx_cluster, Parameters_.ly_cluster);


    mt19937_64 Generator_(Parameters_.RandomSeed); //for random fields
    // mt19937_64 Generator2_(Parameters_.RandomDisorderSeed); //for random disorder


    MFParams MFParams_(Parameters_,Coordinates_,Generator_);
    Hamiltonian Hamiltonian_(Parameters_,Coordinates_,CoordinatesCluster_,MFParams_);
    Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);




    if(check_Non_Int==true){

        //Parameters_.J_HUND=0.0;
        // Hamiltonian_.InteractionsCreate();
        //  Hamiltonian_.Ham_.print();
        // Hamiltonian_.Check_up_down_symmetry();
        //Hamiltonian_.Check_Hermiticity();
        // Hamiltonian_.Diagonalize('V');
        // int temp=Parameters_.Total_Particles;
        //cout<<"mu for n=4 = "<<0.5*(Hamiltonian_.eigs_[temp-1] + Hamiltonian_.eigs_[temp])<<"   "<<
        // Hamiltonian_.eigs_[temp-1]<<"   "<<Hamiltonian_.eigs_[temp]<<endl;

        // double mu_temp = 0.5*(Hamiltonian_.eigs_[temp-1] + Hamiltonian_.eigs_[temp]);
        // Parameters_.mus = Hamiltonian_.chemicalpotential(mu_temp, Parameters_.Total_Particles);
        // cout<<"Energies near the Fermi level:"<<Hamiltonian_.eigs_[temp-1]<<"   "<<Hamiltonian_.eigs_[temp]<<endl;
        //  cout<<"Chemical potential = "<<Parameters_.mus<<endl;
        //  double Quantum_E=Hamiltonian_.E_QM();
        // double Classical_E=Hamiltonian_.GetCLEnergy();
        // cout<<setprecision(9);
        //  cout<<"Total_Energy = "<<Quantum_E+Classical_E<<endl;
        //double mu = chemicalpotential(0.5, temp);
        // Observables_.Get_Non_Interacting_dispersion();
        //Hamiltonian_.Ham_.print();
        //Observables_.Calculate_Akw();
        //Observables_.Calculate_Akw_at_w(mu);
        // Observables_.Calculate_Nw();
    }

    else{

        if(ex_string=="SPA"){
            cout<<setprecision(9);
            Observables_.Initialize();
            MCEngine MCEngine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            MCEngine_.RUN_MC();
        }
        else if (ex_string=="les"){
            if(!Parameters_.Read_Auxilliary_fields){
                cout<<"Read_Auxilliary_fields must be true"<<endl;
                assert(Parameters_.Read_Auxilliary_fields);
            }
            cout<<"Auxilliary fields are read from "<<Parameters_.Auxilliary_fields_Seed_file_name<<endl;






            for(int microstate_no=0;microstate_no<Parameters_.No_Of_Microscopic_States;microstate_no++){

                string DOFS_input = Parameters_.Auxilliary_fields_Seed_file_name + to_string(microstate_no) + ".txt";
                MFParams_.Read_classical_DOFs(DOFS_input);
                Parameters_.Dflag = 'V';
                Hamiltonian_.InteractionsCreate();
                Hamiltonian_.Diagonalize(Parameters_.Dflag);
                double initial_mu_guess;
                int n_states_occupied_zeroT;
                n_states_occupied_zeroT=Parameters_.ns*Parameters_.Fill*Parameters_.orbs*2.0;
                initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
                Parameters_.mus=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Fill);




                Observables_.Calculate_Single_Particle_Density_Matrix();
                Observables_.Calculate_two_point_correlations();
                Observables_.Calculate_2_point_Structure_factors();
                Observables_.Calculate_Nw_t2g();


                Observables_.Two_point_structure_factors_Average();
                Observables_.Nw_t2g_Average();


            }



            char temp_char[50];
            sprintf(temp_char,"%.4f",Parameters_.Temperature);




            //Average Nw_t2g------------------------------------------------//
            //---------------------------------------------------------------//
            string file_Nw_avg = "Thermal_Avg_Nw_t2g_Temp"+ string(temp_char) +".txt";
            ofstream file_Nw_avg_out(file_Nw_avg.c_str());
            file_Nw_avg_out<<"#(w-mu)     yz_up    xz_up      xy_up     ";
            file_Nw_avg_out<<"yz_dn    xz_dn      xy_dn   std_deviations....."<<endl;

            for(int i=0;i<Observables_.Thermal_avg_Nw_t2g[0].size();i++){
                file_Nw_avg_out<<Parameters_.w_min + (i*Parameters_.dw_dos);
                for(int type=0;type<6;type++){
                    file_Nw_avg_out<<"      "<< (Observables_.Thermal_avg_Nw_t2g[type][i])/(1.0*Parameters_.No_Of_Microscopic_States);
                }
                for(int type=0;type<6;type++){

                    file_Nw_avg_out<<"      "<<sqrt(
                                         (
                                             (   (Observables_.Thermal_avg_Nw_t2g_sqr[type][i])/(1.0*Parameters_.No_Of_Microscopic_States)    ) -
                                             ((Observables_.Thermal_avg_Nw_t2g[type][i]*Observables_.Thermal_avg_Nw_t2g[type][i] )/
                                              (Parameters_.No_Of_Microscopic_States*Parameters_.No_Of_Microscopic_States*1.0) )
                                             )
                                         );

                }
                file_Nw_avg_out<<endl;
            }
            //---------------------------------------------------------------//
            //---------------------------------------------------------------//



            //Avr_structure_factors--------------------------------------------//
            //---------------------------------------------------------------//
            string file_Sq_Lq_avg = "Thermal_Avg_Sq_Lq_Temp"+ string(temp_char) +".txt";
            ofstream file_Sq_Lq_avg_out(file_Sq_Lq_avg.c_str());
            file_Sq_Lq_avg_out<<"#qx   qy     S(q)   L(q)   std_dev(Sq)   std_dev(Lq)"<<endl;

            for(int qx=0;qx<Parameters_.lx;qx++){

                for(int qy=0;qy<Parameters_.ly;qy++){

                    file_Sq_Lq_avg_out<<qx<<"     "<<qy<<"      ";
                    file_Sq_Lq_avg_out<<(Observables_.SiSjQ_Mean(qx,qy)/(1.0*Parameters_.No_Of_Microscopic_States)).real()<<"       ";
                    file_Sq_Lq_avg_out<<(Observables_.LiLjQ_Mean(qx,qy)/(1.0*Parameters_.No_Of_Microscopic_States)).real()<<"       ";


                    file_Sq_Lq_avg_out<<
                                         sqrt(
                                             (
                                                 ((Observables_.SiSjQ_square_Mean(qx,qy))/(1.0*Parameters_.No_Of_Microscopic_States))
                                                 -
                                                 ((Observables_.SiSjQ_Mean(qx,qy)*Observables_.SiSjQ_Mean(qx,qy) )/
                                                  (Parameters_.No_Of_Microscopic_States*Parameters_.No_Of_Microscopic_States*1.0) )
                                                 ).real()
                                             )<<"       ";

                    file_Sq_Lq_avg_out<<
                                         sqrt(
                                             (
                                                 ((Observables_.LiLjQ_square_Mean(qx,qy))/(1.0*Parameters_.No_Of_Microscopic_States))
                                                 -
                                                 ((Observables_.LiLjQ_Mean(qx,qy)*Observables_.LiLjQ_Mean(qx,qy) )/
                                                  (Parameters_.No_Of_Microscopic_States*Parameters_.No_Of_Microscopic_States*1.0) )
                                                 ).real()
                                             )<<endl;

                }

                file_Sq_Lq_avg_out<<endl;
            }



            file_Sq_Lq_avg_out<<"# Avg_S2 and std.dev = "<<Observables_.Thermal_Avg_local_S2/(1.0*Parameters_.No_Of_Microscopic_States)<<"     ";
            file_Sq_Lq_avg_out<<sqrt(
                                    (
                                        ((Observables_.Thermal_Avg_local_S2_sqr)/(1.0*Parameters_.No_Of_Microscopic_States))
                                        -
                                        ((Observables_.Thermal_Avg_local_S2*Observables_.Thermal_Avg_local_S2 )/
                                         (Parameters_.No_Of_Microscopic_States*Parameters_.No_Of_Microscopic_States*1.0) )
                                        )
                                    )
                              <<endl;

            file_Sq_Lq_avg_out<<"# Avg_L2 and std.dev = "<<Observables_.Thermal_Avg_local_L2/(1.0*Parameters_.No_Of_Microscopic_States)<<"     ";
            file_Sq_Lq_avg_out<<sqrt(
                                    (
                                        ((Observables_.Thermal_Avg_local_L2_sqr)/(1.0*Parameters_.No_Of_Microscopic_States))
                                        -
                                        ((Observables_.Thermal_Avg_local_L2*Observables_.Thermal_Avg_local_L2 )/
                                         (Parameters_.No_Of_Microscopic_States*Parameters_.No_Of_Microscopic_States*1.0) )
                                        )
                                    )
                              <<endl;

            //---------------------------------------------------------------//
            //---------------------------------------------------------------//






        }

    }




    cout << "--------THE END--------" << endl;
} // main
