#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{
public:


    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator_(Generator__)
    {
        //setupm_arr();
        initialize();
    }

    // Define Fields
    Matrix<double> Stheta, Sphi, S_moment_size;
    Matrix<double> Ltheta, Lphi, L_moment_size;
    Matrix<complex<double>> Rho_den;


    Matrix<double> Disorder;

    double random();
    void FieldThrow(int site, double density_site);
    void FieldThrow(int site, double density_site, int Field_type);
    void initialize();
    void Adjust_MCWindow();
    void Calculate_Fields_Avg();
    void Read_classical_DOFs(string filename);


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator_;
    int lx_,ly_,ns_;
    uniform_real_distribution<double> dis_;
    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};


void MFParams::Adjust_MCWindow(){
    double ratio;
    ratio=Parameters_.AccCount[0]/(Parameters_.AccCount[0]+Parameters_.AccCount[1]);
    //cout<<"ratio= "<< ratio << "temp= "<<Parameters_.temp << endl;
    Parameters_.AccCount[0]=0;
    Parameters_.AccCount[1]=0;
    Parameters_.WindowSize *= abs(1.0 + 1.0*(ratio-0.5));

    if(Parameters_.WindowSize < 0.1){
        Parameters_.WindowSize =0.2;
    }
    //Parameters_.WindowSize = 1.0;
    cout << "Ratio: " << ratio << "  window size:  "<<Parameters_.WindowSize<< endl;
    return;
} // ----------


void MFParams::FieldThrow(int site, double density_at_site){
    int a,b;

    double Pi=Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    double temp_rho_real, temp_rho_imag;

    a=Coordinates_.indx(site);
    b=Coordinates_.indy(site);

    //Angles of S
    Sphi(a,b) += 2*Pi*(random()-0.5)*MC_Window;
    if( Sphi(a,b) < 0.0) {Sphi(a,b) += 2.0*Pi; }
    if( Sphi(a,b) >=2.0*Pi) {Sphi(a,b) -= 2.0*Pi;}

    Stheta(a,b) +=  Pi*(random()-0.5)*MC_Window;
    if ( Stheta(a,b) < 0.0 ) {
        Stheta(a,b) = - Stheta(a,b);
        Sphi(a,b) = fmod( Sphi(a,b)+Pi, 2.0*Pi );
    }
    if ( Stheta(a,b) > Pi ) {
        Stheta(a,b) -= 2.0*Pi ;
        Sphi(a,b) = fmod( Sphi(a,b) + Pi, 2.0*Pi );
    }

    //Angles of L
    Lphi(a,b) += 2*Pi*(random()-0.5)*MC_Window;
    if( Lphi(a,b) < 0.0) {Lphi(a,b) += 2.0*Pi; }
    if( Lphi(a,b) >=2.0*Pi) {Lphi(a,b) -= 2.0*Pi;}

    Ltheta(a,b) +=  Pi*(random()-0.5)*MC_Window;
    if ( Ltheta(a,b) < 0.0 ) {
        Ltheta(a,b) = - Ltheta(a,b);
        Lphi(a,b) = fmod( Lphi(a,b)+Pi, 2.0*Pi );
    }
    if ( Ltheta(a,b) > Pi ) {
        Ltheta(a,b) -= 2.0*Pi ;
        Lphi(a,b) = fmod( Lphi(a,b) + Pi, 2.0*Pi );
    }

    //Moment sizes
    S_moment_size(a,b) += (Parameters_.S_moment_max - Parameters_.S_moment_min)*
            (random()-0.5)*MC_Window;
    L_moment_size(a,b) += (Parameters_.L_moment_max - Parameters_.L_moment_min)*
            (random()-0.5)*MC_Window;

    //Rho
    if(Parameters_.Use_Saddle_point_Approx_on_density){
        assert(!Parameters_.Use_constant_local_Rho_den);
        temp_rho_real = 0.0;
        temp_rho_imag = 2.0*density_at_site;
    }
    if(!Parameters_.Use_Saddle_point_Approx_on_density && !Parameters_.Use_constant_local_Rho_den){
        temp_rho_real = 0.0;//Rho_den(a,b).real() + (Parameters_.Rho_max - Parameters_.Rho_min)*
        //(random()-0.5)*MC_Window;
        temp_rho_imag = Rho_den(a,b).imag() + (Parameters_.Rho_max - Parameters_.Rho_min)*
                (random()-0.5)*MC_Window;
    }
    if(Parameters_.Use_constant_local_Rho_den){
        assert(!Parameters_.Use_Saddle_point_Approx_on_density);
        temp_rho_real = 0.0;
        temp_rho_imag = -2.0*Parameters_.Fill*6.0;

    }
    Rho_den(a,b) = complex<double> (temp_rho_real,temp_rho_imag);

} // ----------


void MFParams::FieldThrow(int site, double density_at_site, int Field_type){
    int a,b;

    double Pi=Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    double temp_rho_real, temp_rho_imag;

    a=Coordinates_.indx(site);
    b=Coordinates_.indy(site);

    //Angles of S
    if(Field_type==0){
        Sphi(a,b) += 2*Pi*(random()-0.5)*MC_Window;
        if( Sphi(a,b) < 0.0) {Sphi(a,b) += 2.0*Pi; }
        if( Sphi(a,b) >=2.0*Pi) {Sphi(a,b) -= 2.0*Pi;}

        Stheta(a,b) +=  Pi*(random()-0.5)*MC_Window;
        if ( Stheta(a,b) < 0.0 ) {
            Stheta(a,b) = - Stheta(a,b);
            Sphi(a,b) = fmod( Sphi(a,b)+Pi, 2.0*Pi );
        }
        if ( Stheta(a,b) > Pi ) {
            Stheta(a,b) -= 2.0*Pi ;
            Sphi(a,b) = fmod( Sphi(a,b) + Pi, 2.0*Pi );
        }
    }

    //Angles of L
    if(Field_type==1){
        Lphi(a,b) += 2*Pi*(random()-0.5)*MC_Window;
        if( Lphi(a,b) < 0.0) {Lphi(a,b) += 2.0*Pi; }
        if( Lphi(a,b) >=2.0*Pi) {Lphi(a,b) -= 2.0*Pi;}

        Ltheta(a,b) +=  Pi*(random()-0.5)*MC_Window;
        if ( Ltheta(a,b) < 0.0 ) {
            Ltheta(a,b) = - Ltheta(a,b);
            Lphi(a,b) = fmod( Lphi(a,b)+Pi, 2.0*Pi );
        }
        if ( Ltheta(a,b) > Pi ) {
            Ltheta(a,b) -= 2.0*Pi ;
            Lphi(a,b) = fmod( Lphi(a,b) + Pi, 2.0*Pi );
        }
    }

    //Moment sizes
    if(Field_type==2){
        S_moment_size(a,b) += (Parameters_.S_moment_max - Parameters_.S_moment_min)*
                (random()-0.5)*MC_Window;
    }
    if(Field_type==3){
        L_moment_size(a,b) += (Parameters_.L_moment_max - Parameters_.L_moment_min)*
                (random()-0.5)*MC_Window;
    }

    //Rho
    if(Field_type==4){
        if(Parameters_.Use_Saddle_point_Approx_on_density){
            assert(!Parameters_.Use_constant_local_Rho_den);
            temp_rho_real = 0.0;
            temp_rho_imag = 2.0*density_at_site;
        }
        if(!Parameters_.Use_Saddle_point_Approx_on_density && !Parameters_.Use_constant_local_Rho_den){
            temp_rho_real = 0.0;//Rho_den(a,b).real() + (Parameters_.Rho_max - Parameters_.Rho_min)*
            //(random()-0.5)*MC_Window;
            temp_rho_imag = Rho_den(a,b).imag() + (Parameters_.Rho_max - Parameters_.Rho_min)*
                    (random()-0.5)*MC_Window;
        }
        if(Parameters_.Use_constant_local_Rho_den){
            assert(!Parameters_.Use_Saddle_point_Approx_on_density);
            temp_rho_real = 0.0;
            temp_rho_imag = -2.0*Parameters_.Fill*6.0;

        }
        Rho_den(a,b) = complex<double> (temp_rho_real,temp_rho_imag);
    }

} // ----------

double MFParams::random(){

    /*
    double random_double;
    random_double=(rand()%RAND_MAX);
    random_double=random_double/RAND_MAX;

    return random_double;
    */

    return dis_(Generator_);

}

void MFParams::initialize(){


    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;

    ns_ = lx_*ly_;
    // srand(Parameters_.RandomSeed);


    Stheta.resize(lx_,ly_);
    Sphi.resize(lx_,ly_);
    S_moment_size.resize(lx_,ly_);

    Ltheta.resize(lx_,ly_);
    Lphi.resize(lx_,ly_);
    L_moment_size.resize(lx_,ly_);

    Rho_den.resize(lx_,ly_);

    Disorder.resize(lx_,ly_);
    Disorder.fill(0.0);


    if(!Parameters_.Read_Auxilliary_fields){
        for(int j=0;j<ly_;j++){
            for(int i=0;i<lx_;i++){
                //ephi(i,j)=(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                //etheta(i,j)=0.5*Parameters_.pi + grnd()*0.2;

                //q=(pi,pi)
                // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                // etheta(i,j)=0.5*(pow(-1.0,j+i)  + 1.0 )*PI ;//+ grnd()*0.2;

                //q=(0,pi)
                //ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                //etheta(i,j)=0.5*(pow(-1.0,j)  + 1.0 )*PI; //+ grnd()*0.2;

                //q=(0,0)
                // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                // etheta(i,j)=0.0; //+ grnd()*0.2;

                //RANDOM

                Sphi(i,j)=2.0*random()*PI;
                Stheta(i,j)=random()*PI;
                S_moment_size(i,j)=Parameters_.S_moment_min + (random()*
                                                               (Parameters_.S_moment_max - Parameters_.S_moment_min) );

                Lphi(i,j)=2.0*random()*PI;
                Ltheta(i,j)=random()*PI;
                L_moment_size(i,j)=Parameters_.L_moment_min + (random()*
                                                               (Parameters_.L_moment_max - Parameters_.L_moment_min) );

                if(Parameters_.Use_constant_local_Rho_den == true){
                    Rho_den(i,j) = complex<double>(0.0,-2.0*Parameters_.Fill*6.0);
                }
                else{
                    Rho_den(i,j).real(Parameters_.Rho_min + (random()*
                                                             (Parameters_.Rho_max - Parameters_.Rho_min) ));

                    Rho_den(i,j).imag(Parameters_.Rho_min + (random()*
                                                             (Parameters_.Rho_max - Parameters_.Rho_min) ));
                }
            }
        }
    }
    else{
        assert(Parameters_.Read_Auxilliary_fields);
        Read_classical_DOFs(Parameters_.Auxilliary_fields_Seed_file_name);
    }


} // ----------

void MFParams::Calculate_Fields_Avg(){

    //NOT_USED

} // ----------


void MFParams::Read_classical_DOFs(string filename){

    string tmp_str;
    int temp_lx, temp_ly, temp_site;
    double temp_rho_real, temp_rho_imag;

    ifstream fl_in(filename.c_str());
    getline(fl_in, tmp_str);

    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){

            fl_in>>temp_lx>>temp_ly>>temp_site>>S_moment_size(i,j)>>Stheta(i,j)>>Sphi(i,j)>>
                                                                    L_moment_size(i,j)>>Ltheta(i,j)>>Lphi(i,j)>>temp_rho_real>>temp_rho_imag;
            Rho_den(i,j)=complex<double>(temp_rho_real,temp_rho_imag);

            assert(temp_lx==i);
            assert(temp_ly==j);
            assert( temp_site==i+(j*lx_) );
        }
    }


} // ----------


#endif
