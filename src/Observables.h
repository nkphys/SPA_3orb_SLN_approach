#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"
#include "functions.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

//n, a, lda, ipvt, work, lwork, info
extern "C" void   zgetri_(int *,std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);

//zgetrf_ (&n, &n, &(_TEMP(0,0)), &n, &(ipvt[0]), &info);
extern "C" void   zgetrf_(int *,int *, std::complex<double> *, int *, int *, int *);


//zhetri (character UPLO, integer N, complex*16 dimension( lda, * ) A, integer LDA,
//integer IPIV, complex*16 dimension( * ) WORK, integer INFO)
extern "C" void   zhetri_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *);


//zhetrf(uplo,n,a,lda,ipiv,work,lwork,info)
extern "C" void   zhetrf_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);



class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {
        Initialize();
    }

    void Initialize();
    void Calculate_Akw_t2g();
    void Calculate_Nw_t2g();
    void Calculate_Single_Particle_Density_Matrix();
    complex<double> Two_particle_Den_Mat(int _alpha, int _beta, int _gamma, int _delta);
    void Calculate_two_point_correlations();
    void Calculate_2_point_Structure_factors();
    void Two_point_structure_factors_Average();
    void Two_point_Average();
    void Total_Energy_Average(double Curr_QuantE, double CurrE);
    complex<double> DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right);
    double Lorentzian(double x, double brd);
    void Fields_SiSj_LiLj();
    void Fields_Two_point_Average();
    void Nw_t2g_Average();


    //Quantum ObservablesXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    complex<double> SiSjQ(int i, int j);
    complex<double> SiSj(int i, int j);
    complex<double> LiLjQ(int i, int j);
    complex<double> LiLj(int i, int j);
    complex<double> SiSjQ_Mean(int i, int j);
    complex<double> SiSjQ_square_Mean(int i, int j);
    complex<double> SiSj_Mean(int i, int j);
    complex<double> SiSj_square_Mean(int i, int j);
    complex<double> LiLjQ_Mean(int i, int j);
    complex<double> LiLjQ_square_Mean(int i, int j);
    complex<double> LiLj_Mean(int i, int j);
    complex<double> LiLj_square_Mean(int i, int j);
    Mat_2_doub Nw_t2g;
    Mat_2_doub Thermal_avg_Nw_t2g;
    Mat_2_doub Thermal_avg_Nw_t2g_sqr;

    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<complex<double>> SiSj_, SiSj_Mean_, SiSj_square_Mean_;
    Matrix<complex<double>> LiLjQ_, LiLjQ_Mean_, LiLjQ_square_Mean_;
    Matrix<complex<double>> LiLj_, LiLj_Mean_, LiLj_square_Mean_;
    Mat_1_Complex_doub Local_density_, Local_density_Mean_, Local_density_square_Mean_;
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



    //Classical observablesXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    Matrix<complex<double>> Fields_SiSjQ_, Fields_SiSjQ_Mean_, Fields_SiSjQ_square_Mean_;
    Matrix<double> Fields_SiSj_, Fields_SiSj_Mean_, Fields_SiSj_square_Mean_;
    Matrix<complex<double>> Fields_LiLjQ_, Fields_LiLjQ_Mean_, Fields_LiLjQ_square_Mean_;
    Matrix<double> Fields_LiLj_, Fields_LiLj_Mean_, Fields_LiLj_square_Mean_;
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


    double AVG_Total_Energy, AVG_Total_Energy_sqr;
    double Avg_local_L2, Avg_local_S2;
    double Thermal_Avg_local_L2, Thermal_Avg_local_S2;
    double Thermal_Avg_local_L2_sqr, Thermal_Avg_local_S2_sqr;



    Mat_2_Complex_doub SP_Density_Matrix;

    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_;


};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/


void Observables::Initialize(){

    complex<double> zero(0.0,0.0);

    AVG_Total_Energy_sqr=0.0;
    AVG_Total_Energy=0.0;

    SiSj_.resize(ns_,ns_);
    SiSj_Mean_.resize(ns_,ns_);
    SiSj_square_Mean_.resize(ns_,ns_);
    LiLj_.resize(ns_,ns_);
    LiLj_Mean_.resize(ns_,ns_);
    LiLj_square_Mean_.resize(ns_,ns_);

    SiSjQ_Mean_.resize(lx_,ly_);
    SiSjQ_square_Mean_.resize(lx_,ly_);
    SiSjQ_.resize(lx_,ly_);
    LiLjQ_Mean_.resize(lx_,ly_);
    LiLjQ_square_Mean_.resize(lx_,ly_);
    LiLjQ_.resize(lx_,ly_);


    Fields_SiSj_.resize(lx_,ly_);
    Fields_SiSj_Mean_.resize(lx_,ly_);
    Fields_SiSj_square_Mean_.resize(lx_,ly_);
    Fields_LiLj_.resize(lx_,ly_);
    Fields_LiLj_Mean_.resize(lx_,ly_);
    Fields_LiLj_square_Mean_.resize(lx_,ly_);

    Fields_SiSjQ_Mean_.resize(lx_,ly_);
    Fields_SiSjQ_square_Mean_.resize(lx_,ly_);
    Fields_SiSjQ_.resize(lx_,ly_);
    Fields_LiLjQ_Mean_.resize(lx_,ly_);
    Fields_LiLjQ_square_Mean_.resize(lx_,ly_);
    Fields_LiLjQ_.resize(lx_,ly_);


    Local_density_.resize(6*ns_);
    Local_density_Mean_.resize(6*ns_);
    Local_density_square_Mean_.resize(6*ns_);

    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            SiSj_Mean_(i,j)=zero;
            SiSj_square_Mean_(i,j)=zero;
            LiLj_Mean_(i,j)=zero;
            LiLj_square_Mean_(i,j)=zero;
        }
    }

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            SiSjQ_Mean_(ix,iy)=zero;
            SiSjQ_square_Mean_(ix,iy)=zero;
            LiLjQ_Mean_(ix,iy)=zero;
            LiLjQ_square_Mean_(ix,iy)=zero;

            Fields_SiSj_Mean_(ix,iy)=0.0;
            Fields_SiSj_square_Mean_(ix,iy)=0.0;
            Fields_LiLj_Mean_(ix,iy)=0.0;
            Fields_LiLj_square_Mean_(ix,iy)=0.0;

            Fields_SiSjQ_Mean_(ix,iy)=zero;
            Fields_SiSjQ_square_Mean_(ix,iy)=zero;
            Fields_LiLjQ_Mean_(ix,iy)=zero;
            Fields_LiLjQ_square_Mean_(ix,iy)=zero;
        }
    }

    for(int i=0;i<6*ns_;i++){
        Local_density_[i]=zero;
        Local_density_Mean_[i]=zero;
        Local_density_square_Mean_[i]=zero;
    }

    Thermal_Avg_local_L2=0.0;
    Thermal_Avg_local_S2=0.0;
    Thermal_Avg_local_L2_sqr=0.0;
    Thermal_Avg_local_S2_sqr=0.0;

    Thermal_avg_Nw_t2g.resize(6);
    Thermal_avg_Nw_t2g_sqr.resize(6);

    int omega_index_size = int( (Parameters_.w_max - Parameters_.w_min)/(Parameters_.dw_dos) );
    for(int i=0;i<6;i++){
     Thermal_avg_Nw_t2g[i].resize(omega_index_size);
     Thermal_avg_Nw_t2g_sqr[i].resize(omega_index_size);
    }





} // ----------

complex<double> Observables::DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right){
    complex<double> temp_;
    temp_=zero_complex;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += conj(left[i])*right[i];
    }
    return temp_;

}

void Observables::Calculate_Akw_t2g(){


    //---------Read from input file-----------------------//
    string fileout="Akw_t2g.txt";
    double omega_min, omega_max, d_omega;
    /*
    double eta = 0.08;
    omega_min=Hamiltonian_.eigs_[0]-0.5-Parameters_.mus;omega_max=Hamiltonian_.eigs_[6*ns_ -1]+0.5-Parameters_.mus;d_omega=0.0005;
    */
    double eta = 0.001;
    omega_min=-1.6;omega_max=2.6;d_omega=0.03;
    //---------------------------------------------------//


    int UP_=0;
    int DN_=1;
    int YZ_, XZ_, XY_;
    YZ_=0; XZ_=1; XY_=2;

    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());

    int c1,c2;

    Mat_3_Complex_doub A_YZ, A_XZ,A_XY;
    A_YZ.resize(Parameters_.ns);
    A_XZ.resize(Parameters_.ns);
    A_XY.resize(Parameters_.ns);

    for (int i=0;i<Parameters_.ns;i++){
        A_YZ[i].resize(Parameters_.ns);
        A_XZ[i].resize(Parameters_.ns);
        A_XY[i].resize(Parameters_.ns);

        for(int j=0;j<Parameters_.ns;j++){
            A_YZ[i][j].resize(omega_index_max);
            A_XZ[i][j].resize(omega_index_max);
            A_XY[i][j].resize(omega_index_max);
        }
    }


    complex<double> Nup_check(0,0);
    complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            cout<<"Akw for "<<l<<"  "<<j<<" done"<<endl;
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
                A_YZ[j][l][omega_ind]=zero_complex;
                A_XZ[j][l][omega_ind]=zero_complex;
                A_XY[j][l][omega_ind]=zero_complex;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                    //c= l + or1*ns_ + ns_*orbs_*spin;
                    for(int spin=0;spin<2;spin++){
                        //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                        c1 = Coordinates_.Nc_dof(l,YZ_ + 3*spin);
                        c2 = Coordinates_.Nc_dof(j,YZ_ + 3*spin);
                        A_YZ[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                        c1 = Coordinates_.Nc_dof(l,XZ_ + 3*spin);
                        c2 = Coordinates_.Nc_dof(j,XZ_ + 3*spin);
                        A_XZ[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                        c1 = Coordinates_.Nc_dof(l,XY_ + 3*spin);
                        c2 = Coordinates_.Nc_dof(j,XY_ + 3*spin);
                        A_XY[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);
                    }
                }
            }
        }
    }

    cout << "Nup_check = "<<Nup_check<<endl;
    cout << "Ndn_check = "<<Ndn_check<<endl;

    complex<double> temp_YZ, temp_XZ, temp_XY;

    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    if(Parameters_.ly==1){
        //--------For 1D in x direction-----------------
        ky_i=0;
        for(kx_i=-1;kx_i<=(Parameters_.lx);kx_i++){
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------
    }
    else if(Parameters_.lx==1){
        //--------For 1D in y direction-----------------
        kx_i=0;
        for(ky_i=-1;ky_i<=(Parameters_.ly);ky_i++){
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------
    }
    else{
        assert(Parameters_.lx!=1);
        assert(Parameters_.ly!=1);

        //--------\Gamma to X-----------------
        ky_i=0;
        for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------

        //--------X to M-----------------
        kx_i=(Parameters_.lx/2);
        for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------

        //--------M to \Gamma[with one extra point,
        //                  because in gnuplor use "set pm3d corners2color c1"
        //                  ]-----------------
        kx_i=(Parameters_.lx/2) - 1;
        ky_i=(Parameters_.lx/2) - 1;
        for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
            temp_pair.first = kx_i;
            temp_pair.second = kx_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


    }


    file_Akw_out<< "#k_point   kx_i   ky_i   (ky_i*Parameters_.lx) + kx_i    omega_min + (d_omega*omega_ind)   omega_ind    Akw_YZ.real()    Akw_XZ.real()    Akw_XY.real()    _YZ.imag()    _XZ.imag()    _XY.imag()"<<endl;


    double k22_offset=1*PI;
    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp_YZ=zero_complex;
            temp_XZ=zero_complex;
            temp_XY=zero_complex;


            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp_YZ += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_YZ[j][l][omega_ind];


                    temp_XZ += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_XZ[j][l][omega_ind];

                    temp_XY += one_complex*
                            exp(iota_complex*((kx+k22_offset)*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              (ky+k22_offset)*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_XY[j][l][omega_ind];

                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "
                        <<temp_YZ.real()<<"    "<<temp_XZ.real()<<"    "<<temp_XY.real()<<"    "
                       <<temp_YZ.imag()<<"    "<<temp_XZ.imag()<<"    "<<temp_XY.imag()<<"    "<<endl;

        }
        file_Akw_out<<endl;
    }



}

void Observables::Calculate_Nw_t2g(){


    Nw_t2g.resize(6);

    //---------Read from input file-----------------------//
    double omega_min, omega_max, d_omega;
    double eta = Parameters_.eta_dos;
    omega_min= Parameters_.w_min;
    omega_max= Parameters_.w_max;
    d_omega=Parameters_.dw_dos;
    //---------------------------------------------------//

    int c1;
    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );
    double temp_val;
    for(int i=0;i<6;i++){
    Nw_t2g[i].resize(omega_index_max);
    }

    //---------------------------------------------------------------------------------//
    //************************************Nw_jm****************************************//
    //---------------------------------------------------------------------------------//
    /*
    ofstream file_Nw_out(fileout.c_str());


    file_Nw_out<<"#(w-mu)    jm_3by2_m3by2     jm_3by2_3by2     ";
    file_Nw_out<<"jm_3by2_m1by2       jm_3by2_1by2     jm_1by2_m1by2   jm_1by2_1by2"<<endl;

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        file_Nw_out<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int state_type=0;state_type<6;state_type++){
            temp_val=0.0;

            for(int site=0;site<Coordinates_.ns_;site++){
                c1=Coordinates_.Nc_dof_(site,state_type);

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    temp_val +=  (conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                                  Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                }
            }
            file_Nw_out<<temp_val<<"      ";
        }
        file_Nw_out<<endl;
    }

    file_Nw_out<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

    */
    //---------------------------------------------------------------------------------//
    //*********************************************************************************//



    //---------------------------------------------------------------------------------//
    //************************************Nw_jm****************************************//
    //---------------------------------------------------------------------------------//

/*
    string fileout_t2g="Nw_t2g_" + Parameters_.Tag_output_files_with +   ".txt";
    ofstream file_Nw_out_t2g(fileout_t2g.c_str());
    file_Nw_out_t2g<<"#(w-mu)    yz_up    xz_up      xy_up     ";
    file_Nw_out_t2g<<"yz_dn    xz_dn      xy_dn"<<endl;
    */

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
       // file_Nw_out_t2g<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int t2g_type=0;t2g_type<6;t2g_type++){
            temp_val=0.0;

            for(int site=0;site<Coordinates_.ns_;site++){
                c1=Coordinates_.Nc_dof_(site,t2g_type);

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    temp_val += ( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                                  Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                }

            }

            Nw_t2g[t2g_type][omega_ind]=temp_val;
           // file_Nw_out_t2g<<temp_val<<"      ";
        }

       // file_Nw_out_t2g<<endl;
    }

   // file_Nw_out_t2g<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

    //---------------------------------------------------------------------------------//
    //*********************************************************************************//






    /*
    int n_chosen=Parameters_.Total_Particles - 1;
    file_Nw_out_t2g<<"#Gap = "<<Hamiltonian_.eigs_[n_chosen+1] - Hamiltonian_.eigs_[n_chosen]<<endl;

    string fileout_Eigen="Eigen_spectrum_"  + Parameters_.Tag_output_files_with + ".txt";
    ofstream file_Eigen_out(fileout_Eigen.c_str());

    file_Eigen_out<<"#n  E[n]   Fermi[n]"<<endl;
    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
        file_Eigen_out<<n<<"\t"<<Hamiltonian_.eigs_[n]<<"\t"<<(1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)) <<endl;
    }
    */


}

double Observables::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}

complex<double> Observables::SiSjQ(int i,int j){return SiSjQ_(i,j);}
complex<double> Observables::SiSj(int i,int j){return SiSj_(i,j);}

complex<double> Observables::LiLjQ(int i,int j){return LiLjQ_(i,j);}
complex<double> Observables::LiLj(int i,int j){return LiLj_(i,j);}

complex<double> Observables::SiSjQ_Mean(int i,int j){return SiSjQ_Mean_(i,j);}
complex<double> Observables::SiSjQ_square_Mean(int i,int j){return SiSjQ_square_Mean_(i,j);}

complex<double> Observables::LiLjQ_Mean(int i,int j){return LiLjQ_Mean_(i,j);}
complex<double> Observables::LiLjQ_square_Mean(int i,int j){return LiLjQ_square_Mean_(i,j);}

complex<double> Observables::SiSj_Mean(int i,int j){return SiSj_Mean_(i,j);}
complex<double> Observables::SiSj_square_Mean(int i,int j){return SiSj_square_Mean_(i,j);}

complex<double> Observables::LiLj_Mean(int i,int j){return LiLj_Mean_(i,j);}
complex<double> Observables::LiLj_square_Mean(int i,int j){return LiLj_square_Mean_(i,j);}


void Observables::Fields_SiSj_LiLj()
{

    double Cos_ij, Sin_ij, ei, ai, phase;
    double ei_L, ai_L;

    int site_, site_p, ax, ay;

    Mat_1_doub sx_, sz_, sy_;
    sx_.resize(ns_);
    sy_.resize(ns_);
    sz_.resize(ns_);

    Mat_1_doub Lx_, Lz_, Ly_;
    Lx_.resize(ns_);
    Ly_.resize(ns_);
    Lz_.resize(ns_);

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            site_ = Coordinates_.Nc(i, j);
            ei = MFParams_.Stheta(i, j);
            ai = MFParams_.Sphi(i, j);
            sx_[site_] = MFParams_.S_moment_size(i, j)*cos(ai) * sin(ei); //*MFParams_.S_moment_size(i, j)* ;
            sy_[site_] = MFParams_.S_moment_size(i, j)*sin(ai) * sin(ei);
            sz_[site_] = MFParams_.S_moment_size(i, j)*cos(ei);

            ei_L = MFParams_.Ltheta(i, j);
            ai_L = MFParams_.Lphi(i, j);
            Lx_[site_] = MFParams_.L_moment_size(i, j) *cos(ai_L) * sin(ei_L);//MFParams_.L_moment_size(i, j) *
            Ly_[site_] = MFParams_.L_moment_size(i, j) *sin(ai_L) * sin(ei_L);
            Lz_[site_] = MFParams_.L_moment_size(i, j) *cos(ei_L);
        }
    }

    for (int xr = 0; xr < lx_; xr++)
    {
        for (int yr = 0; yr < ly_; yr++)
        {
            Fields_SiSj_(xr, yr) = 0.0;
            Fields_LiLj_(xr, yr) = 0.0;
            for (int i = 0; i < lx_; i++)
            {
                for (int j = 0; j < ly_; j++)
                {
                    site_ = Coordinates_.Nc(i, j);
                    ax = (i + xr) % lx_;
                    ay = (j + yr) % ly_;
                    site_p = Coordinates_.Nc(ax, ay);
                    Fields_SiSj_(xr, yr) += sx_[site_] * sx_[site_p];
                    Fields_SiSj_(xr, yr) += sy_[site_] * sy_[site_p];
                    Fields_SiSj_(xr, yr) += sz_[site_] * sz_[site_p];

                    Fields_LiLj_(xr, yr) += Lx_[site_] * Lx_[site_p];
                    Fields_LiLj_(xr, yr) += Ly_[site_] * Ly_[site_p];
                    Fields_LiLj_(xr, yr) += Lz_[site_] * Lz_[site_p];
                }
            }
            Fields_SiSj_(xr, yr) *= double(1.0 / (lx_ * ly_));
            Fields_LiLj_(xr, yr) *= double(1.0 / (lx_ * ly_));
            //cout << xr << " "<< yr<< " "<<  SiSj_(xr,yr) << endl;
        }
    }

    for (int qx = 0; qx < lx_; qx++)
    {
        for (int qy = 0; qy < ly_; qy++)
        {
            Fields_SiSjQ_(qx, qy) = zero_complex;
            Fields_LiLjQ_(qx, qy) = zero_complex;
            for (int xr = 0; xr < lx_; xr++)
            {
                for (int yr = 0; yr < ly_; yr++)
                {
                    phase = 2.0 * Parameters_.pi * (double(qx * xr) / double(lx_) + double(qy * yr) / double(ly_));
                    Cos_ij = cos(phase);
                    Sin_ij = sin(phase);
                    Fields_SiSjQ_(qx, qy) += Fields_SiSj_(xr, yr) * complex<double>(Cos_ij, Sin_ij);

                    Fields_LiLjQ_(qx, qy) += Fields_LiLj_(xr, yr) * complex<double>(Cos_ij, Sin_ij);
                }
            }
            Fields_SiSjQ_(qx, qy) *= double(1.0 / (lx_ * ly_));
            Fields_LiLjQ_(qx, qy) *= double(1.0 / (lx_ * ly_));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

} // ----------


void Observables::Fields_Two_point_Average(){

    for(int i=0; i<lx_; i++) {
        for(int j=0; j<ly_; j++) {
            Fields_SiSj_Mean_(i,j) += Fields_SiSj_(i,j);
            Fields_SiSj_square_Mean_(i,j) += ( Fields_SiSj_(i,j)*Fields_SiSj_(i,j) )   ;

            Fields_LiLj_Mean_(i,j) += Fields_LiLj_(i,j);
            Fields_LiLj_square_Mean_(i,j) += ( Fields_LiLj_(i,j)*Fields_LiLj_(i,j) )   ;

            Fields_SiSjQ_Mean_(i,j) += Fields_SiSjQ_(i,j);
            Fields_SiSjQ_square_Mean_(i,j) += ( Fields_SiSjQ_(i,j)*Fields_SiSjQ_(i,j) )   ;

            Fields_LiLjQ_Mean_(i,j) += Fields_LiLjQ_(i,j);
            Fields_LiLjQ_square_Mean_(i,j) += ( Fields_LiLjQ_(i,j)*Fields_LiLjQ_(i,j) )   ;

            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }
} // ----------


void Observables::Calculate_Single_Particle_Density_Matrix(){

    /*
      NOTE:
      SP_Density_Matrix[alpha][beta] = <c_{alpha^{daggger}} c_{beta}>
     */


    SP_Density_Matrix.resize(ns_*6);
    for(int i=0;i<ns_*6;i++){
        SP_Density_Matrix[i].resize(ns_*6);
    }

    for(int alpha_=0;alpha_<ns_*6;alpha_++){
        for(int beta_=0;beta_<ns_*6;beta_++){
            SP_Density_Matrix[alpha_][beta_] = zero_complex;
            for(int n=0;n<ns_*6;n++){
                SP_Density_Matrix[alpha_][beta_] += conj(Hamiltonian_.Ham_(alpha_,n))*Hamiltonian_.Ham_(beta_,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));
            }
        }
    }

    for(int i=0;i<6*ns_;i++){
        Local_density_[i]=SP_Density_Matrix[i][i] ;
    }



}

complex<double> Observables::Two_particle_Den_Mat(int _alpha, int _beta, int _gamma, int _delta){

    complex<double> temp;
    complex<double> delta_gamma_beta;

    if(_gamma == _beta){
        delta_gamma_beta=one_complex;
    }
    else{
        assert(_gamma != _beta);
        delta_gamma_beta=zero_complex;
    }

    temp = (SP_Density_Matrix[_alpha][_beta]*SP_Density_Matrix[_gamma][_delta])
            +
            (SP_Density_Matrix[_alpha][_delta]*(delta_gamma_beta - SP_Density_Matrix[_gamma][_beta]));
    return temp;
}

void Observables::Calculate_two_point_correlations(){

    int YZ_,XZ_,XY_;
    YZ_=0; XZ_=1; XY_=2;

    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    int _Sz, _Sx,_Sy, _Lz, _Lx,_Ly;
    int _Jz_eff, _Jx_eff,_Jy_eff;
    int _M_Total_z, _M_Total_x,_M_Total_y;

    Mat_1_string opr_type;
    opr_type.clear();
    opr_type.push_back("Sz");_Sz=0;
    opr_type.push_back("Sx");_Sx=1;
    opr_type.push_back("Sy");_Sy=2;
    opr_type.push_back("Lz");_Lz=3;
    opr_type.push_back("Lx");_Lx=4;
    opr_type.push_back("Ly");_Ly=5;
    opr_type.push_back("Jz_eff");_Jz_eff=6;
    opr_type.push_back("Jx_eff");_Jx_eff=7;
    opr_type.push_back("Jy_eff");_Jy_eff=8;
    opr_type.push_back("M_Total_z");_M_Total_z=9;
    opr_type.push_back("M_Total_x");_M_Total_x=10;
    opr_type.push_back("M_Total_y");_M_Total_y=11;


    Mat_3_Complex_doub Oprs_;
    Oprs_.resize(opr_type.size());
    for(int opr_no_=0;opr_no_<opr_type.size();opr_no_++){
        Oprs_[opr_no_].clear();
        Oprs_[opr_no_].resize(6);
        for(int i=0;i<6;i++){
            Oprs_[opr_no_][i].resize(6);
        }
    }




    //Opr. =\sum_{orb1,orb2,spin1,spin2}Temp_Mat{orb1,spin1;orb2,spin2}c*_{orb1,spin1}c_{orb2,spin2}
    //    for(int orb1=0;orb1<3;orb1++){
    //        for(int spin1=0;spin1<2;spin1++){
    //            //orb1 + 3*spin1
    //            for(int orb2=0;orb2<3;orb2++){
    //                for(int spin2=0;spin2<2;spin2++){
    //                }
    //            }
    //        }
    //    }

    int opr_no;

    //Sz---
    assert(opr_type[_Sz]=="Sz");
    for(int orb1=0;orb1<3;orb1++){
        for(int spin1=0;spin1<2;spin1++){
            Oprs_[_Sz][orb1 + 3*spin1][orb1 + 3*spin1]=one_complex*(0.5*(1.0-(2.0*spin1)));
        }}


    //Sx---
    assert(opr_type[_Sx]=="Sx");
    for(int orb1=0;orb1<3;orb1++){
        Oprs_[_Sx][orb1 + 3*0][orb1 + 3*1]=one_complex*(0.5);
        Oprs_[_Sx][orb1 + 3*1][orb1 + 3*0]=one_complex*(0.5);
    }

    //Sy---
    assert(opr_type[_Sy]=="Sy");
    for(int orb1=0;orb1<3;orb1++){
        Oprs_[_Sy][orb1 + 3*0][orb1 + 3*1]=(-1.0*iota_complex)*(0.5);
        Oprs_[_Sy][orb1 + 3*1][orb1 + 3*0]=(-1.0*iota_complex)*(-0.5);
    }

    //Lz---
    assert(opr_type[_Lz]=="Lz");
    for(int spin1=0;spin1<2;spin1++){
        Oprs_[_Lz][YZ_ + 3*spin1][XZ_ + 3*spin1]=iota_complex;
        Oprs_[_Lz][XZ_ + 3*spin1][YZ_ + 3*spin1]=iota_complex*(-1.0);
    }

    //Lx---
    opr_no=_Lx;
    assert(opr_type[opr_no]=="Lx");
    for(int spin1=0;spin1<2;spin1++){
        Oprs_[opr_no][XZ_ + 3*spin1][XY_ + 3*spin1]=iota_complex;
        Oprs_[opr_no][XY_ + 3*spin1][XZ_ + 3*spin1]=iota_complex*(-1.0);
    }

    //Ly---
    opr_no=_Ly;
    assert(opr_type[opr_no]=="Ly");
    for(int spin1=0;spin1<2;spin1++){
        Oprs_[opr_no][XY_ + 3*spin1][YZ_ + 3*spin1]=iota_complex;
        Oprs_[opr_no][YZ_ + 3*spin1][XY_ + 3*spin1]=iota_complex*(-1.0);
    }

    //J_eff = S-L----------
    assert(opr_type[_Jz_eff]=="Jz_eff");
    assert(opr_type[_Jx_eff]=="Jx_eff");
    assert(opr_type[_Jy_eff]=="Jy_eff");
    assert(opr_type[_M_Total_z]=="M_Total_z");
    assert(opr_type[_M_Total_x]=="M_Total_x");
    assert(opr_type[_M_Total_y]=="M_Total_y");
    for(int ir=0;ir<6;ir++){
        for(int ic=0;ic<6;ic++){
            Oprs_[_Jz_eff][ir][ic]=Oprs_[_Sz][ir][ic] - Oprs_[_Lz][ir][ic];
            Oprs_[_Jx_eff][ir][ic]=Oprs_[_Sx][ir][ic] - Oprs_[_Lx][ir][ic];
            Oprs_[_Jy_eff][ir][ic]=Oprs_[_Sy][ir][ic] - Oprs_[_Ly][ir][ic];
            Oprs_[_M_Total_z][ir][ic]=2.0*Oprs_[_Sz][ir][ic] + Oprs_[_Lz][ir][ic];
            Oprs_[_M_Total_x][ir][ic]=2.0*Oprs_[_Sx][ir][ic] + Oprs_[_Lx][ir][ic];
            Oprs_[_M_Total_y][ir][ic]=2.0*Oprs_[_Sy][ir][ic] + Oprs_[_Ly][ir][ic];
        }
    }




    //  string corrs_out = "corrs_" + Parameters_.Tag_output_files_with + ".txt";
    //  ofstream file_corrs_out(corrs_out.c_str());
    //  file_corrs_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)     SS[site_i][site_j]     LL[site_i][site_j]     JeffJeff[site_i][site_j]     MM[site_i][site_j]"<<endl;

    int i1,i2,j1,j2;

    complex<double> temp_val_SS, temp_val_LL, temp_val_JeffJeff, temp_val_MM;
    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){

            temp_val_SS = zero_complex;
            temp_val_LL = zero_complex;
            temp_val_JeffJeff = zero_complex;
            temp_val_MM = zero_complex;

            for(int i_row=0;i_row<6;i_row++){
                for(int i_col=0;i_col<6;i_col++){

                    for(int j_row=0;j_row<6;j_row++){
                        for(int j_col=0;j_col<6;j_col++){

                            //S.S
                            for(int comp=_Sz;comp<_Sy+1;comp++){
                                if( (Oprs_[comp][i_row][i_col] !=zero_complex)
                                        &&
                                        (Oprs_[comp][j_row][j_col] !=zero_complex)   ){
                                    i1=Coordinates_.Nc_dof(i,i_row);
                                    i2=Coordinates_.Nc_dof(i,i_col);
                                    j1=Coordinates_.Nc_dof(j,j_row);
                                    j2=Coordinates_.Nc_dof(j,j_col);
                                    temp_val_SS += Oprs_[comp][i_row][i_col]*Oprs_[comp][j_row][j_col]*
                                            Two_particle_Den_Mat(i1,i2,j1,j2);
                                }
                            }

                            //L.L
                            for(int comp=_Lz;comp<_Ly+1;comp++){
                                if( (Oprs_[comp][i_row][i_col] !=zero_complex)
                                        &&
                                        (Oprs_[comp][j_row][j_col] !=zero_complex)   ){
                                    i1=Coordinates_.Nc_dof(i,i_row);
                                    i2=Coordinates_.Nc_dof(i,i_col);
                                    j1=Coordinates_.Nc_dof(j,j_row);
                                    j2=Coordinates_.Nc_dof(j,j_col);
                                    temp_val_LL += Oprs_[comp][i_row][i_col]*Oprs_[comp][j_row][j_col]*
                                            Two_particle_Den_Mat(i1,i2,j1,j2);
                                }
                            }

                            /*
                            //J_eff.J_eff
                            for(int comp=_Jz_eff;comp<_Jy_eff+1;comp++){
                                if( (Oprs_[comp][i_row][i_col] !=zero_complex)
                                        &&
                                        (Oprs_[comp][j_row][j_col] !=zero_complex)   ){
                                    i1=Coordinates_.Nc_dof(i,i_row);
                                    i2=Coordinates_.Nc_dof(i,i_col);
                                    j1=Coordinates_.Nc_dof(j,j_row);
                                    j2=Coordinates_.Nc_dof(j,j_col);
                                    temp_val_JeffJeff += Oprs_[comp][i_row][i_col]*Oprs_[comp][j_row][j_col]*
                                            Two_particle_Den_Mat(i1,i2,j1,j2);
                                }
                            }

                            //M.M
                            for(int comp=_M_Total_z;comp<_M_Total_y+1;comp++){
                                if( (Oprs_[comp][i_row][i_col] !=zero_complex)
                                        &&
                                        (Oprs_[comp][j_row][j_col] !=zero_complex)   ){
                                    i1=Coordinates_.Nc_dof(i,i_row);
                                    i2=Coordinates_.Nc_dof(i,i_col);
                                    j1=Coordinates_.Nc_dof(j,j_row);
                                    j2=Coordinates_.Nc_dof(j,j_col);
                                    temp_val_MM += Oprs_[comp][i_row][i_col]*Oprs_[comp][j_row][j_col]*
                                            Two_particle_Den_Mat(i1,i2,j1,j2);
                                }
                            }
                            */

                        }
                    }
                }
            }
            SiSj_(i,j)=temp_val_SS;
            LiLj_(i,j)=temp_val_LL;
            //  file_corrs_out<<i<<setw(15)<<Coordinates_.indx(i)<<setw(15)<<Coordinates_.indy(i)<<setw(15)<<j<<setw(15)<<Coordinates_.indx(j)<<setw(15)<<Coordinates_.indy(j)<<
            //                  setw(15)<<temp_val_SS.real()<<//"\t"<<temp_val_SS.imag()<<
            //                  setw(15)<<temp_val_LL.real()<<//"\t"<<temp_val_LL.imag()<<
            //                  setw(15)<<temp_val_JeffJeff.real()<<//"\t"<<temp_val_JeffJeff.imag()<<
            //                  setw(15)<<temp_val_MM.real()<<endl;//"\t"<<temp_val_MM.imag()<<endl;

        }
    }


    //------------------------------//
}

void Observables::Calculate_2_point_Structure_factors(){


    int ix, iy, jx, jy;
    int xr, yr;
    double phase, Cos_ij, Sin_ij;

   /*
    string Sq_Lq_out = "Sq_Lq_" + Parameters_.Tag_output_files_with + ".txt";
    ofstream file_Sq_Lq_out(Sq_Lq_out.c_str());
    file_Sq_Lq_out<<"#qx      qy     S(qx,qy)     L(qx,qy)"<<endl;
   */

    for(int qx=0; qx<lx_; qx++) {
        for(int qy=0; qy<ly_; qy++) {
            SiSjQ_(qx,qy)=complex<double>(0.0,0.0);
            LiLjQ_(qx,qy)=complex<double>(0.0,0.0);
            for(int i=0;i<ns_;i++){
                for(int j=0;j<ns_;j++){

                    ix = Coordinates_.indx(i);iy = Coordinates_.indy(i);
                    jx = Coordinates_.indx(j);jy = Coordinates_.indy(j);
                    xr=jx-ix;
                    yr=jy-iy;

                    phase=2.0*Parameters_.pi*(double(qx*xr)/double(lx_)+double(qy*yr)/double(ly_));
                    Cos_ij = cos(phase);
                    Sin_ij = sin(phase);
                    SiSjQ_(qx,qy) += SiSj_(i,j)*complex<double>(Cos_ij,Sin_ij);
                    LiLjQ_(qx,qy) += LiLj_(i,j)*complex<double>(Cos_ij,Sin_ij);

                }
            }
            SiSjQ_(qx,qy)*= double(1.0/(lx_*ly_*lx_*ly_));
            LiLjQ_(qx,qy)*= double(1.0/(lx_*ly_*lx_*ly_));

            //file_Sq_Lq_out << qx << "     "<< qy<< "     "<<  SiSjQ_(qx,qy).real() << "      "<<LiLjQ_(qx,qy).real()<<endl;
        }
        //file_Sq_Lq_out <<endl;
    }


    Avg_local_L2=0.0;
    Avg_local_S2=0.0;
    for(int i=0;i<ns_;i++){
       Avg_local_L2 +=  LiLj_(i,i).real();
       Avg_local_S2 +=  SiSj_(i,i).real();
    }
    Avg_local_L2 = Avg_local_L2/(1.0*ns_);
    Avg_local_S2 = Avg_local_S2/(1.0*ns_);

    /*
    file_Sq_Lq_out<<"# Avg_local_S2 = "<<Avg_local_S2<<endl;
    file_Sq_Lq_out<<"# Avg_local_L2 = "<<Avg_local_L2<<endl;
    */

}

void Observables::Two_point_Average(){

    for(int i=0; i<ns_; i++) {
        for(int j=0; j<ns_; j++) {
            SiSj_Mean_(i,j) += SiSj_(i,j);
            SiSj_square_Mean_(i,j) += ( SiSj_(i,j)*SiSj_(i,j) )   ;

            LiLj_Mean_(i,j) += LiLj_(i,j);
            LiLj_square_Mean_(i,j) += ( LiLj_(i,j)*LiLj_(i,j) )   ;
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }
} // ----------


void Observables::Two_point_structure_factors_Average(){

    for(int qx=0; qx<lx_; qx++) {
        for(int qy=0; qy<ly_; qy++) {
            SiSjQ_Mean_(qx,qy) += SiSjQ_(qx,qy);
            SiSjQ_square_Mean_(qx,qy) += ( SiSjQ_(qx,qy)*SiSjQ_(qx,qy) )   ;

            LiLjQ_Mean_(qx,qy) += LiLjQ_(qx,qy);
            LiLjQ_square_Mean_(qx,qy) += ( LiLjQ_(qx,qy)*LiLjQ_(qx,qy) )   ;
        }
    }

    Thermal_Avg_local_S2  +=  Avg_local_S2;
    Thermal_Avg_local_L2  +=  Avg_local_L2;
    Thermal_Avg_local_S2_sqr  +=  Avg_local_S2*Avg_local_S2;
    Thermal_Avg_local_L2_sqr  +=  Avg_local_L2*Avg_local_L2;

} // ----------

void  Observables::Nw_t2g_Average(){

    int omega_index_size = int( (Parameters_.w_max - Parameters_.w_min)/(Parameters_.dw_dos) );


    for(int type=0;type<6;type++){
        for(int i=0;i<omega_index_size;i++){
            Thermal_avg_Nw_t2g[type][i] += Nw_t2g[type][i];
            Thermal_avg_Nw_t2g_sqr[type][i] += Nw_t2g[type][i]*Nw_t2g[type][i];
        }
    }


}

void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE){

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE)*(Curr_QuantE + CurrE);
}


#endif // OBSERVABLES_H
