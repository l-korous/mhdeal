#include "complete_elliptic_integrals.h"
#include "initialCondition.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialCondition<equationsType, dim>::InitialCondition(Parameters<dim>& parameters) : Function<dim>(Equations<equationsType, dim>::n_components), parameters(parameters)
{
};


template <EquationsType equationsType, int dim>
void InitialCondition<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points, std::vector<Vector<double> > & result) const
{
}

template <EquationsType equationsType, int dim>
Parameters<dim>& InitialCondition<equationsType, dim>::getParams() const
{
  return parameters;
}


 /***************************************************************************
          MHD Blast initial condition
  ***************************************************************************/
template <EquationsType equationsType, int dim>
MHDBlastIC<equationsType, dim>::MHDBlastIC(Parameters<dim>& parameters) : 
                                           InitialCondition<equationsType,dim>(parameters)
{
};

 
template <EquationsType equationsType, int dim>
void MHDBlastIC<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points,
                                                  std::vector<Vector<double> > &result) const
{
  // The result is a two-dimensional array, first dimension is for the integration point, second for the component (density, momentum-x, ...)
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 1.;
    result[i][1] = 0.;
    result[i][2] = 0.;
    result[i][3] = 0.;
    result[i][4] = 1.0 / std::sqrt(2.);
    result[i][5] = 1.0 / std::sqrt(2.);
    result[i][6] = 0.;
    if (points[i].norm() < 0.1)
        result[i][7] = 10. / (this->getParams().gas_gamma - 1.0) + 0.5;
    else
        result[i][7] = 0.1 / (this->getParams().gas_gamma - 1.0) + 0.5;
  }
}


 /***************************************************************************
          Titov & Demoulin initial condition
  ***************************************************************************/
template <EquationsType equationsType, int dim>
TitovDemoulinIC<equationsType, dim>::TitovDemoulinIC(Parameters<dim> &parameters) : 
               InitialCondition<equationsType,dim>(parameters)
{
  // plasma beta
    beta=0.05;
    Lg=0.0;          

    invLg=0.0;
    if(Lg > 0.0) invLg=1.0/Lg;
    
    //======================== TD-model specific parameters

    // Torus winding number
    N_t=-3.0;

    // Torus major radius
    R_t=4.0;
  
    // Submerging of torus main axis in units of R_t
    d2R_t=2.0/R_t;

    // Distance of magnetic charges from x=0 plane in units of R_t
    L2R_t=2.0/R_t;
      
    
    //======================= Calculate dependent TD model parameters

    // Normalised magnetic charge corresponding to global equilibrium (Eq. 6)
    q_mag=0.25*fabs(N_t)*(log(8.0*R_t)-1.25)
      *(1+L2R_t*L2R_t)*sqrt(1+L2R_t*L2R_t)/L2R_t;

    // Sign of winding: corresponds to sign of I_O in TD paper
    iSgn=(N_t>=0)?1.0:-1.0;

    // "Helicity" factor inside tho loop (used later in B_theta_internal calcs)
    heliFactor=2.0*(N_t*N_t)/(R_t*R_t);


    //======================= Parameters for plasma

    // The coronal/prominence temperature ratio and its inverse value.
    Tc2Tp=1.0;
    //  Tp2Tc=1.0/Tc2Tp;

    //======================= TCPR definition

    // density jump half-width...
    t_rho=0.12;
    
    // ...and its inverse value
    densGrad=1.0/t_rho;
}

template <EquationsType equationsType, int dim>
void TitovDemoulinIC<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points,
                                            std::vector<Vector<double> >   &value_list) const
{
  for(unsigned int p=0; p<points.size(); ++p)
    this->point_value(points[p], value_list[p]);
}

 /***************************************************************************
          Calculate the field according to TD paper (A&A 351, 707, 1999)
          Fill the structure with gravity-stratified plasma. 
  ***************************************************************************/
template <EquationsType equationsType, int dim>
void TitovDemoulinIC<equationsType, dim>::point_value(const Point<dim> &p, 
                                                      Vector<double> &result) const
{
    //========== Calculate the vector potential for I_t-generated toroidal field 
    double xx,yy,zz;
    Point<dim> &ca = this->getParams().corner_a;
    Point<dim> &cb = this->getParams().corner_b;
    
    xx=p[0]-(ca[0]+cb[0])*0.5;
    yy=p[1]-(ca[1]+cb[1])*0.5;
    zz=p[2]-ca[2];

    Vector<double> theta0(3);
    double r_maj,r_min;
    const double dd=1e-8;  // precision of numerical derivatives
    const double idd=0.5/dd;
    double df[6][3];
    double P[6][3];
    double dP[3][3];
    
    // calculate vector potentil in 6 close points (+dx,-dx,+dy,-dy,+dz,-dz)
    for(unsigned int i=0;i<6;++i){
      df[i][0]=df[i][1]=df[i][2]=0.0;
      df[i][int(i/2)]=double(1-2*int(i%2))*dd;
    }
    
    for(unsigned int i=0;i<6;++i){
      double x=xx+df[i][0];
      double y=yy+df[i][1];
      double z=zz+df[i][2];
      
      // Distances from torus major and minor axes
      r_maj=sqrt(y*y+(z+d2R_t*R_t)*(z+d2R_t*R_t));
      r_min=sqrt(x*x+(r_maj-R_t)*(r_maj-R_t));

      // Unit vector of toroidal coordinate theta
      theta0[0]=0.0;
      theta0[1]=-(z+d2R_t*R_t)/r_maj;
      theta0[2]=y/r_maj;

      // Common radial factor for A_tor
      double rFactor=fabs(N_t)*sqrt(1.0/(R_t*r_maj));

      // Argument of elliptical integral
      double kr=2.0*sqrt(r_maj*R_t/((r_maj+R_t)*(r_maj+R_t)+x*x));

      //---- Sew-up internal and external solutions
      if(r_min>1.0){ //---- external region 
          // Elliptical integrals
          double Ek,Kk;
          Complete_Elliptic_Integrals_Modulus(kr, Kk, Ek);

          double Ak=((2.0-kr*kr)*Kk-2.0*Ek)/kr;

          P[i][0]=(rFactor*Ak)*theta0[0];
          P[i][1]=(rFactor*Ak)*theta0[1];
          P[i][2]=(rFactor*Ak)*theta0[2];
          
      }else{ //---- inside the torus

          // ka=kr at r_min=1 (=torus surface) 
          double ka=2.0*sqrt(r_maj*R_t/(4.0*r_maj*R_t+1.0));

          double Ek,Kk;
          Complete_Elliptic_Integrals_Modulus(ka, Kk, Ek);

          double Ak=((2.0-ka*ka)*Kk-2.0*Ek)/ka;
          double Ak_prime=((2.0-ka*ka)*Ek-2.0*(1.0-ka*ka)*Kk)/
                          (ka*ka*(1.0-ka*ka));

          double cf=(rFactor*(Ak+Ak_prime*(kr-ka)));
          P[i][0]=cf*theta0[0];
          P[i][1]=cf*theta0[1];
          P[i][2]=cf*theta0[2];
      }
    }
    
    // calculate derivatives of vector potential
    for(unsigned int i=0;i<3;++i){
      dP[i][0]=(P[2*i][0]-P[2*i+1][0])*idd;
      dP[i][1]=(P[2*i][1]-P[2*i+1][1])*idd;
      dP[i][2]=(P[2*i][2]-P[2*i+1][2])*idd;
    }
    
    //====================== Calculate the full state field

    double pressure;
    
    // Distances from torus major and minor axes
    r_maj=sqrt(yy*yy+(zz+d2R_t*R_t)*(zz+d2R_t*R_t));
    r_min=sqrt(xx*xx+(r_maj-R_t)*(r_maj-R_t));

    // Unit vector of toroidal coordinate theta
    theta0[0]=0.0;
    theta0[1]=-(zz+d2R_t*R_t)/r_maj;
    theta0[2]=yy/r_maj;

    // Radius vectors originating in magnetic charges
    Vector<double> r_plus(3);//(xx-L2R_t*R_t,yy,zz+d2R_t*R_t);
    r_plus[0]=xx-L2R_t*R_t;
    r_plus[1]=yy;
    r_plus[2]=zz+d2R_t*R_t;
    Vector<double> r_minus(3);//(xx+L2R_t*R_t,yy,zz+d2R_t*R_t);
    r_minus[0]=xx+L2R_t*R_t;
    r_minus[1]=yy;
    r_minus[2]=zz+d2R_t*R_t;

    double rp=r_plus.l2_norm();
    double rm=r_minus.l2_norm();

    //========== Calculate the magnetic field in TD equilibrium by parts

    //--- Q-generated part of field

    double cf1=q_mag/(rp*rp*rp);
    double cf2=q_mag/(rm*rm*rm);
    Vector<double> B_loc(3);//=cf1*r_plus-cf2*r_minus;
    B_loc.sadd(0.0,cf1,r_plus,-cf2,r_minus);
 
    // add vector potential part B = curl A
    B_loc[0]+=dP[1][2]-dP[2][1];
    B_loc[1]+=dP[2][0]-dP[0][2];
    B_loc[2]+=dP[0][1]-dP[1][0];
    
    /*
      barta@asu.cas.cz
      10/02/2012

      With the following density-scaling the velocities are normalised 
      to V_A outside the flux-rope; L_g is the coronal height-scale.

      ----------

      The density (and consequent temperature) jump rho~H(r_min-1) 
      replaced by smoother rho~tgh(r_min-1) profile (TPCR-like).
    */

    double rho_0=0.5*(1.0-Tc2Tp)*tanh(densGrad*(r_min-1.0))+0.5*(1+Tc2Tp);
    
    if(r_min>1.0){ // external region

      result[0]=rho_0*exp(-zz*invLg);         // mass density outside
      
      B_loc.sadd(1.0,iSgn*R_t/r_maj,theta0);

      pressure=beta*exp(-zz*invLg);
    }else{ // inside the torus

      result[0]=rho_0*exp(-zz*Tc2Tp*invLg);   // mass density in the loop

      B_loc.sadd(1.0,(iSgn*(sqrt(1.0+heliFactor*(1.0-r_min*r_min))+
                    R_t/r_maj-1.0)),theta0);

      pressure=beta*exp(-zz*invLg);
    }
    
    result[1]=0.0;                    
    result[2]=0.0;                          // momentum density
    result[3]=0.0;
    
    result[4]=B_loc[0];
    result[5]=B_loc[1];                  // magnetic field
    result[6]=B_loc[2];
    // energy density
    result[7]=pressure/(this->getParams().gas_gamma-1.0)+
              result[4]*result[4]+result[5]*result[5]+result[6]*result[6];
    
//     v[8]=0.0;                    
//     v[9]=0.0;                          // current density
//     v[10]=0.0;
//     
//     v[11]=0.0; // eta
}

template class InitialCondition<EquationsTypeMhd, 3>;
template class MHDBlastIC<EquationsTypeMhd, 3>;
template class TitovDemoulinIC<EquationsTypeMhd, 3>;
