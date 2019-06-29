#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "tensor_type.h"

#ifndef MCENGINE_H
#define MCENGINE_H


class MCEngine{
public:
    MCEngine(Parameters& Parameters__, Coordinates& Coordinates__,
             MFParams& MFParams__, Hamiltonian& Hamiltonian__,
             Observables& Observables__)
        : Parameters_(Parameters__),Coordinates_(Coordinates__),
          MFParams_(MFParams__), Hamiltonian_(Hamiltonian__),
          Observables_(Observables__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns),
          orbs_(Parameters_.orbs), ED_(Parameters_.ED_)
    {

    }


    void RUN_MC();
    double Prob (double muu, double mu_new);
    double ProbCluster (double muu, double mu_new);
    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    Observables &Observables_;
    const int lx_, ly_, ns_, orbs_;
    bool ED_;

};

/*
 * ***********
 *  Functions in Class MCEngine ------
 *  ***********
*/

void MCEngine::RUN_MC(){

    complex<double> zero(0.0,0.0);
    double local_density_;
    int int_temp;
    ED_ = Parameters_.ED_;

    int MC_sweeps_used_for_Avg=Parameters_.Last_n_sweeps_for_measurement;
    int Gap_bw_sweeps = Parameters_.Measurement_after_each_m_sweeps;

    double Prev_Classical_E,Curr_Classical_E,P_new,P12, Prob_,muu_prev;
    double muu_prevCluster;
    double Prev_QuantE;
    double Curr_QuantE;
    double Curr_QuantECluster;
    double Prev_QuantECluster;
    int x,y,act;
    double saved_Params[8];
    //int Stheta_, Sphi_, S_moment_size_, Ltheta_, Lphi_,  L_moment_size_,  Rho_den_real_, Rho_den_imag_;
    enum Field_type {Stheta_, Sphi_, S_moment_size_, Ltheta_, Lphi_,  L_moment_size_,  Rho_den_real_, Rho_den_imag_};



    string File_Out_progress;
    string File_Out_theta_phi;

    double temp_=Parameters_.temp_max;

    double initial_mu_guess;
    int n_states_occupied_zeroT;



    //starting with a random guess

    while(temp_>=Parameters_.temp_min){

        cout << "Temperature = " << temp_<<" is being done"<<endl;
        Parameters_.Temperature=temp_;
        Parameters_.beta=double(1.0/Parameters_.Temperature);

        for(int ix=0;ix<lx_;ix++){
            for(int iy=0;iy<ly_;iy++){
                Observables_.SiSjQ_Mean_(ix,iy)=zero;
                Observables_.SiSjQ_square_Mean_(ix,iy)=zero;
                Observables_.LiLjQ_Mean_(ix,iy)=zero;
                Observables_.LiLjQ_square_Mean_(ix,iy)=zero;

                Observables_.Fields_LiLj_(ix,iy)=0.0;
                Observables_.Fields_LiLjQ_(ix,iy)=0.0;
                Observables_.Fields_LiLjQ_Mean_(ix,iy)=0.0;
                Observables_.Fields_LiLjQ_square_Mean_(ix,iy)=0.0;

                Observables_.Fields_SiSj_(ix,iy)=0.0;
                Observables_.Fields_SiSjQ_(ix,iy)=0.0;
                Observables_.Fields_SiSjQ_Mean_(ix,iy)=0.0;
                Observables_.Fields_SiSjQ_square_Mean_(ix,iy)=0.0;
            }
        }
        for(int i=0;i<ns_;i++){
            for(int j=0;j<ns_;j++){
                Observables_.SiSj_square_Mean_(i,j)=zero;
                Observables_.SiSj_Mean_(i,j)=zero;
                Observables_.LiLj_square_Mean_(i,j)=zero;
                Observables_.LiLj_Mean_(i,j)=zero;
            }
        }

        Observables_.AVG_Total_Energy=0.0;
        Observables_.AVG_Total_Energy_sqr=0.0;


        char temp_char[50];
        sprintf(temp_char,"%.4f",temp_);

        File_Out_theta_phi = "Auxilliary_Fields_Temp" + string(temp_char) + ".txt";
        ofstream File_Out_Theta_Phi(File_Out_theta_phi.c_str());

        File_Out_progress = "output_Temp" + string(temp_char) + ".txt";
        ofstream file_out_progress(File_Out_progress.c_str());

        file_out_progress<< "Total "<<Parameters_.IterMax<<" sweeps are performed."<<endl;
        file_out_progress<<"First "<<Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)<<
                           " sweeps are used for thermalization and every "<<Gap_bw_sweeps+1<<" in last "<<
                           Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg<<
                           " sweeps are used for measurement."<<endl;
        act=1;




        Parameters_.WindowSize = 0.1; //2f + 0.003f*beta0 ;
        Parameters_.Eav=0.0;
        Parameters_.Dflag='N'; // flag to calculate only Eigenvalue
        //std::string name="Output/Conf_" + to_string(ltemp) + ".dat";
        //Parameters_.beta = double(11604.0/ (Parameters_.temp +20.0) );
        //cout << "TEMP  " << Parameters_.temp << endl;


        file_out_progress<<"I_MC"<<setw(15)<<"S(0,Pi)"<<setw(15)<<"S(Pi,0)"<<setw(17)<<"S(0,0)"<<setw(17)<<"S(Pi,Pi)"<<setw(17)<<
                           "S(Pi/2,Pi/2)"<<setw(17)<<"< N_total(or N_cluster) >"<<setw(15)<<"E_CL"<<setw(15)<<"E_QM"<<setw(15)<<"mu (or mu_cluster)"<< endl;


        if(ED_){
        Prev_Classical_E = Hamiltonian_.GetCLEnergy();
        if(Parameters_.J_Hund !=0.0 || Parameters_.U_onsite!=0.0){
        Hamiltonian_.InteractionsCreate();
        Hamiltonian_.Diagonalize(Parameters_.Dflag);
        n_states_occupied_zeroT=Parameters_.ns*Parameters_.Fill*Parameters_.orbs*2.0;
        initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
        muu_prev=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Fill);
        Parameters_.mus = muu_prev;
        Prev_QuantE = Hamiltonian_.E_QM();
        Hamiltonian_.copy_eigs(1);
        }
        else{
            assert(Parameters_.J_Hund ==0.0 && Parameters_.U_onsite==0.0);
            Prev_QuantE =0.0;
            muu_prev=0.0;
        }

        cout<<"Initial Classical Energy[Full System] = "<<Prev_Classical_E<<endl;
        cout<<"Initial Quantum Energy[Full System] = "<<Prev_QuantE<<endl;
        cout<<"Initial mu_initial full system="<<muu_prev<<endl;

        //copying everything in cluster [cluster is system here.]
        muu_prevCluster=muu_prev;
        Hamiltonian_.HamCluster_ = Hamiltonian_.Ham_;
        Parameters_.mus_Cluster = muu_prev;
        Hamiltonian_.eigsCluster_ = Hamiltonian_.eigs_;
        Hamiltonian_.copy_eigs_Cluster(1);
        Prev_QuantECluster = Prev_QuantE;

        }

        int Confs_used=0;
        int measure_start=0;



        for(int count=0;count<Parameters_.IterMax;count++){

            for(int i=0;i<ns_;i++) {  // For each site


                //***Before change*************//

                if(!ED_){
                    //TCA is used
                    if(Parameters_.Use_Saddle_point_Approx_on_density){
                        cout<<"with TCA, Saddle point Approxiamation is not working at present,"<<endl;
                        cout<<" Either put local_density_field (Rho) to be constant or use ED"<<endl;
                        assert(!Parameters_.Use_Saddle_point_Approx_on_density);
                    }

                    Prev_Classical_E = Hamiltonian_.GetCLEnergy();

                    if(Parameters_.J_Hund !=0.0 || Parameters_.U_onsite!=0.0){
                        Hamiltonian_.InteractionsClusterCreate(i);
                        Hamiltonian_.DiagonalizeCluster(Parameters_.Dflag);
                        //n_states_occupied_zeroT=Parameters_.Fill*Hamiltonian_.eigsCluster_.size();
                        //initial_mu_guess=0.5*(Hamiltonian_.eigsCluster_[n_states_occupied_zeroT-1] + HamiltonianCluster_.eigs_[n_states_occupied_zeroT])
                        muu_prevCluster=Hamiltonian_.chemicalpotentialCluster(muu_prevCluster,Parameters_.Fill);
                        Parameters_.mus_Cluster = muu_prevCluster;
                        Prev_QuantECluster = Hamiltonian_.E_QMCluster();
                        Hamiltonian_.copy_eigs_Cluster(1);
                    }
                    else{
                        assert(Parameters_.J_Hund ==0.0 && Parameters_.U_onsite==0.0);
                        Parameters_.mus_Cluster=0.0;
                        Prev_QuantECluster=0.0;
                        muu_prevCluster=0.0;
                    }
                }
                else{
                    assert(ED_);
                    local_density_=0.0;
                    if(Parameters_.J_Hund !=0.0 || Parameters_.U_onsite!=0.0){
                    if(Parameters_.Use_Saddle_point_Approx_on_density){
                        Hamiltonian_.eigs_ = Hamiltonian_.eigsCluster_;
                        Hamiltonian_.Ham_ = Hamiltonian_.HamCluster_;
                        Parameters_.mus = Parameters_.mus_Cluster;
                        Observables_.Calculate_Single_Particle_Density_Matrix();
                        for(int orb_ind=0;orb_ind<3;orb_ind++){
                            for(int spin_ind=0;spin_ind<2;spin_ind++){
                                int_temp = Coordinates_.Nc_dof(i,orb_ind + 3*spin_ind);
                                local_density_ += Observables_.Local_density_[int_temp].real();
                            }
                        }
                    }
                }

                }

                //*******************************//

                x=Coordinates_.indx(i);
                y=Coordinates_.indy(i);


                saved_Params[Stheta_]=MFParams_.Stheta(x,y);
                saved_Params[Sphi_]=MFParams_.Sphi(x,y);

                saved_Params[S_moment_size_]=MFParams_.S_moment_size(x,y);

                saved_Params[Ltheta_]=MFParams_.Ltheta(x,y);
                saved_Params[Lphi_]=MFParams_.Lphi(x,y);

                saved_Params[L_moment_size_]=MFParams_.L_moment_size(x,y);

                saved_Params[Rho_den_real_]=MFParams_.Rho_den(x,y).real();
                saved_Params[Rho_den_imag_]=MFParams_.Rho_den(x,y).imag();


                MFParams_.FieldThrow(i, local_density_);
                Curr_Classical_E = Hamiltonian_.GetCLEnergy();

                if(Parameters_.J_Hund !=0.0 || Parameters_.U_onsite!=0.0){
                    Hamiltonian_.InteractionsClusterCreate(i);
                    Hamiltonian_.DiagonalizeCluster(Parameters_.Dflag);
                    Parameters_.mus_Cluster=Hamiltonian_.chemicalpotentialCluster(muu_prevCluster,Parameters_.Fill);
                    Curr_QuantECluster = Hamiltonian_.E_QMCluster();
                }
                else{
                    assert(Parameters_.J_Hund ==0.0 && Parameters_.U_onsite==0.0);
                    Parameters_.mus_Cluster=0.0;
                    Curr_QuantECluster=0.0;
                }

                //Ratio of Quantum partition functions
                /*P = [ Tr(exp(-beta(Hquant_new)))/Tr(exp(-beta(Hquant_old)))]*
                      [exp(-beta*E_classical(New)) / exp(-beta*E_classical(old))]
                     * [sin(STheta_i(New)) / sin(STheta_i(Old)) ]
                     * [sin(LTheta_i(New)) / sin(LTheta_i(Old)) ]*/
                /*exp(P12) = P
                  P12 = log (P)
                  */

                //same mu-refrence is used, otherwise engine does not work properly?
                if(Parameters_.J_Hund !=0.0 || Parameters_.U_onsite!=0.0){
                    P_new = ProbCluster(muu_prevCluster, Parameters_.mus_Cluster);
                   // P_new = ProbCluster(muu_prevCluster, muu_prevCluster);
                }
                else{
                    P_new = 0.0;
                }
                P12 = P_new - Parameters_.beta*((Curr_Classical_E)-(Prev_Classical_E));
                //cout<<P12<<endl;
                P12 += log ((sin(MFParams_.Stheta(x,y))/sin(saved_Params[Stheta_])));
                P12 += log ((sin(MFParams_.Ltheta(x,y))/sin(saved_Params[Ltheta_])));

                Prob_ = exp(P12);

                //---OR---
                //P12 = exp(-Parameters_.beta*((CurrE+Curr_QuantE)-(PrevE+Prev_QuantE)));
                //P12*= (sin(MFParams_.etheta(x,y))/sin(saved_Params[0]));

                //Heat bath algorithm [See page-129 of Prof. Elbio's Book]
                //Heat bath algorithm works for small changes i.e. when P12~1.0
                //Heat_Bath_Algo-----
                 // Prob_ =Prob_/(1.0+Prob_);
                //-----------------

                //Metropolis Algotithm--
                 Prob_=min(1.0,Prob_);
                //---------------------



                /*
       * VON NEUMANN's REJECTING METHOD:
       * Random number < P12 -----> ACCEPT
       * Random number > P12 -----> REJECT
       */
                //ACCEPTED
              if ( Prob_ > MFParams_.random() ) {
                    Parameters_.AccCount[0]++;
                    act=1;
                    if(ED_){
                        Prev_Classical_E=Curr_Classical_E;
                        Prev_QuantECluster = Curr_QuantECluster;
                        Hamiltonian_.copy_eigs_Cluster(1);
                        muu_prevCluster=Parameters_.mus_Cluster;
                    }
                }
                //REJECTED
                else{
                    Parameters_.AccCount[1]++;
                    act=0;
                    MFParams_.Stheta(x,y)=saved_Params[Stheta_];
                    MFParams_.Sphi(x,y)=saved_Params[Sphi_];
                    MFParams_.Ltheta(x,y)=saved_Params[Ltheta_];
                    MFParams_.Lphi(x,y)=saved_Params[Lphi_];
                    MFParams_.S_moment_size(x,y)=saved_Params[S_moment_size_];
                    MFParams_.L_moment_size(x,y)=saved_Params[L_moment_size_];
                    MFParams_.Rho_den(x,y).real(saved_Params[Rho_den_real_]);
                    MFParams_.Rho_den(x,y).imag(saved_Params[Rho_den_imag_]);
                }

            }// site loop


            if ( (count%10==0) ) {
                MFParams_.Adjust_MCWindow();
            }




            //No Measurement performed
            if(count < (Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)) ){
                if ( (count%10==0) ) {

                    Observables_.Fields_SiSj_LiLj();
                    file_out_progress << int(1.0*count) <<setw(20)<< Observables_.Fields_SiSjQ_(0,int(lx_/2)).real() <<setw(16)<< Observables_.Fields_SiSjQ_(int(lx_/2),0).real()
                                      <<setw(16)<<Observables_.Fields_SiSjQ_(0,0).real() <<setw(16)<< Observables_.Fields_SiSjQ_(int(lx_/2),int(lx_/2)).real() <<setw(16)<<
                                        Observables_.Fields_SiSjQ_(int(lx_/4),int(lx_/4)).real() <<setw(16)<< Hamiltonian_.ClusterDensity() <<setw(16)<< Curr_Classical_E
                                     <<setw(16)<< Curr_QuantECluster<<setw(15)<<Parameters_.mus_Cluster<< endl;
                }
            }
            //Measurement performed after every "Gap_bw_sweeps"
            else{

                if(measure_start==0){
                    measure_start++;
                    file_out_progress<<"----------Measurement is started----------"<<endl;
                    file_out_progress<<"I_MC      Avg{S(0,0)}    Avg{S(pi,pi)}    std.dev{S(0,0)}   std.dev{S(pi,pi)}   Avg{E_classical}  std.dev{E_classical}"<<endl;
                }
                int temp_count=count -
                        (Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg));
                int zero_or_not = temp_count % (Gap_bw_sweeps + 1);
                if( zero_or_not==0 ){

                    if((Parameters_.Saving_Microscopic_States==true) &&
                            (Confs_used<Parameters_.No_Of_Microscopic_States)
                            ){

                        char Confs_char[50];
                        sprintf(Confs_char,"%d",Confs_used);
                        string File_Out_fields_microState = "Fields_Temp" + string(temp_char) +
                                "MicroState" + string(Confs_char) +".txt";
                        ofstream File_Out_fields_MicroState(File_Out_fields_microState.c_str());


                        File_Out_fields_MicroState<<"#x"<<setw(15)<<"y"<<setw(15)<<"STheta(x,y)"<<setw(15)<<"SPhi(x,y)"<<
                                                    setw(15)<<"S_Moment_Size(x,y)"<<setw(15)<<"LTheta(x,y)"<<setw(15)<<"LPhi(x,y)"<<
                                                    setw(15)<<"L_Moment_Size(x,y)"<<setw(15)<<"Rho_den(x,y).real"<<
                                                    setw(15)<<"Rho_den(x,y).imag"<<endl;
                        for(int ix=0;ix<lx_;ix++){
                            for(int iy =0;iy<ly_;iy++){
                                File_Out_fields_MicroState<<ix<<setw(15)<<iy<<setw(15)<<MFParams_.Stheta(ix,iy)<<setw(15)<<MFParams_.Sphi(ix,iy)<<
                                                            setw(15)<<MFParams_.S_moment_size(ix,iy)<<setw(15)<<MFParams_.Ltheta(ix,iy)<<
                                                            setw(15)<<MFParams_.Lphi(ix,iy)<<setw(15)<<MFParams_.L_moment_size(ix,iy)<<
                                                            setw(15)<<MFParams_.Rho_den(ix,iy).real()<<setw(15)<<MFParams_.Rho_den(ix,iy).imag()<<endl;
                            }}

                    }

                    //IF TCA is used, Diagonalize the Full system for measurement
                    if(!ED_){
                        if(Parameters_.J_Hund !=0.0 || Parameters_.U_onsite!=0.0){
                        Parameters_.Dflag = 'V';
                        Hamiltonian_.InteractionsCreate();
                        Hamiltonian_.Diagonalize(Parameters_.Dflag);
                        Parameters_.mus=Hamiltonian_.chemicalpotential(muu_prevCluster,Parameters_.Fill);
                        Parameters_.Dflag = 'N';
                        }
                   }
                    else{
                        Hamiltonian_.eigs_ = Hamiltonian_.eigsCluster_;
                        Hamiltonian_.Ham_ = Hamiltonian_.HamCluster_;
                        Parameters_.mus = Parameters_.mus_Cluster;
                    }


                    Confs_used=Confs_used+1;
                    Observables_.Fields_SiSj_LiLj();
                    Observables_.Fields_Two_point_Average();
                   // Observables_.Calculate_Single_Particle_Density_Matrix();
                  //  Observables_.Calculate_two_point_correlations();
                   // Observables_.Calculate_2_point_Structure_factors();
                  //  Observables_.Two_point_Average();
                  //  Observables_.Two_point_structure_factors_Average();
                    Curr_QuantE = Hamiltonian_.E_QM();
                    Observables_.Total_Energy_Average( Curr_QuantE, Curr_Classical_E);

                    //double MC_steps_Avg_insitu = (1.0 + 1.0*(count - (Parameters_.IterMax - MC_steps_used_for_Avg)));

                    file_out_progress << int(1.0*count) <<setw(20)<<
                                         Observables_.Fields_SiSjQ_Mean_(0,0).real()/(Confs_used*1.0)
                                      <<setw(16)<<Observables_.Fields_SiSjQ_Mean_(int(lx_/2),int(lx_/2)).real()/(Confs_used*1.0)
                                     <<setw(16)<<


                                       sqrt(
                                           (( Observables_.Fields_SiSjQ_square_Mean_(0,0)/(Confs_used*1.0) ) -
                                            ((Observables_.Fields_SiSjQ_Mean_(0,0)*Observables_.Fields_SiSjQ_Mean_(0,0) )/(Confs_used*Confs_used*1.0) ) ).real()
                                           )

                                    <<setw(16)<<
                                      sqrt(
                                          (( Observables_.Fields_SiSjQ_square_Mean_(int(lx_/2),int(lx_/2))/(Confs_used*1.0) ) -
                                           ((Observables_.Fields_SiSjQ_Mean_(int(lx_/2),int(lx_/2))*Observables_.Fields_SiSjQ_Mean_(int(lx_/2),int(lx_/2)) )/(Confs_used*Confs_used*1.0) ) ).real()
                                          )
                                   <<setw(16)<<

                                     Observables_.AVG_Total_Energy/(Confs_used*1.0)
                                  <<setw(16)<<
                                    sqrt(  (Observables_.AVG_Total_Energy_sqr/(Confs_used*1.0)) -
                                           ((Observables_.AVG_Total_Energy*Observables_.AVG_Total_Energy)/(Confs_used*Confs_used*1.0))  )
                                 <<endl;

                }

            }


        }// count Loop [single sweep over lattice]
        file_out_progress << "Total "<<Confs_used<< " configurations were used were measurement"<<endl;

        temp_ = temp_ - Parameters_.d_Temp;

        //Writing Auxilliary fields for present temperature:
        //File_Out_Theta_Phi
        File_Out_Theta_Phi<<"# lx    ly    site   S_moment_size   Stheta    Sphi   L_moment_size  Ltheta   Lphi   rho_real    rho_imag"<<endl;
        for(int j=0;j<ly_;j++){
            for(int i=0;i<lx_;i++){

               File_Out_Theta_Phi<<i<<"   "<<j<<"   "<<i+(j*lx_)<<"   "<<MFParams_.S_moment_size(i,j)<<"   "<<
                                   MFParams_.Stheta(i,j)<<"   "<<MFParams_.Sphi(i,j)<<"   "<<MFParams_.L_moment_size(i,j)<<"   "
                                <<MFParams_.Ltheta(i,j)<<"   "<<MFParams_.Lphi(i,j)<<"   "
                               <<MFParams_.Rho_den(i,j).real()<<"   "<<MFParams_.Rho_den(i,j).imag()<<endl;
            }

        }

       // MFParams_.Read_classical_DOFs(File_Out_theta_phi);
    }//Temperature loop

} // ---------



double MCEngine::Prob(double muu, double mu_new){

    double P=0.0;
    double X,Y,X2;

    for(int i=0;i<2*orbs_*ns_;i++){
        X = Parameters_.beta*( (mu_new) - Hamiltonian_.eigs_[i]);
        Y = Parameters_.beta*( (muu) - Hamiltonian_.eigs_saved_[i]);
        //P += log(1 + exp(X)) - log(1 + exp(Y));

        if(X>5){
            P +=X;
        }
        else if(fabs(X)<0.001){
            P += log(2.0 + X);
        }
        else if(X<-5){
            P +=exp(X);
        }
        else{
            P +=log(1.0 + exp(X));
        }



        if(Y>5){
            P -=Y;
        }
        else if(fabs(Y)<0.001){
            P -= log(2.0 + Y);
        }
        else if(Y<-5){
            P -=exp(Y);
        }
        else{
            P -=log(1.0 + exp(Y));
        }




    }

    return P;

} // ---------



double MCEngine::ProbCluster(double muu, double mu_new){

    double P=0.0;
    double X,Y,X2;
    int ns = (Parameters_.lx_cluster)*(Parameters_.ly_cluster);

    for(int i=0;i<2*orbs_*ns;i++){
        X = Parameters_.beta*( (mu_new) - Hamiltonian_.eigsCluster_[i]);
        Y = Parameters_.beta*( (muu) - Hamiltonian_.eigsCluster_saved_[i]);
        //P += log(1 + exp(X)) - log(1 + exp(Y));

        if(X>5){
            P +=X;
        }
        else if(fabs(X)<0.001){
            P += log(2.0 + X);
        }
        else if(X<-5){
            P +=exp(X);
        }
        else{
            P +=log(1.0 + exp(X));
        }



        if(Y>5){
            P -=Y;
        }
        else if(fabs(Y)<0.001){
            P -= log(2.0 + Y);
        }
        else if(Y<-5){
            P -=exp(Y);
        }
        else{
            P -=log(1.0 + exp(Y));
        }




    }

    return P;

} // ---------





#endif // MCENGINE_H
