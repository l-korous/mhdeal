#include "common.h"


d EquationImplementation::matrixVolValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, Point<DIM> quadPoint)
{
  d result = 0.;
  // Time derivative.
  if (comp_i == comp_j)
    result += u_val * v_val;

  return result;
}


d EquationImplementation::matrixBoundaryEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}


d EquationImplementation::matrixInternalEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N, 
  Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}



d EquationImplementation::rhsVolValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad, Point<DIM> quadPoint)
{
  d result = 0.;

  // Time derivative.
  result += Un_val[comp_i] * v_val;

  return result;
}


d EquationImplementation::rhsBoundaryEdgeValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}


d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad,
  d v_valN, vec Un_valN, dimVec v_gradN,
  vecDimVec Un_gradN, bool v_N, Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}

void EquationImplementation::Jacobians(FullMatrix<double> *J,
                                    std::vector<Vector<double> > lv,
                                    const unsigned int qp)
{
  double v[11],iRh,iRh2,Uk,p,gmmo,Ec1,Ec2,Ec3,E1,E2,E3;
  
  // using shorter notation for old solution
  // order of the variables is following: rho, v(3), B(3), U, J(3)
  for(unsigned int i=0;i<11;i++)
    v[i]=lv[qp](i);
  
  iRh=1.0/v[0];
  iRh2=iRh*iRh;
  Uk=iRh*(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
  p=gamma*(v[7]-v[4]*v[4]+v[5]*v[5]+v[6]*v[6]-Uk);
  gmmo=gamma-1.0;
  Ec1=(v[3]*v[5]-v[2]*v[6])*iRh;
  Ec1=(v[1]*v[6]-v[3]*v[4])*iRh;
  Ec1=(v[2]*v[4]-v[1]*v[5])*iRh;
  E1=Ec1+eta*v[8];
  E2=Ec2+eta*v[9];
  E3=Ec3+eta*v[10];
  
  J[0](0,0)=0;
  J[0](1,0)=1;
  J[0](2,0)=0;
  J[0](3,0)=0;
  J[0](4,0)=0;
  J[0](5,0)=0;
  J[0](6,0)=0;
  J[0](7,0)=0;
  J[0](8,0)=0;
  J[0](9,0)=0;
  J[0](10,0)=0;

  J[0](0,1)=-(v[1]*v[1]*iRh2)+(gmmo*Uk)*.5*iRh;
  J[0](1,1)=(2*v[1])*iRh-(gmmo*v[1])*iRh;
  J[0](2,1)=-((gmmo*v[2])*iRh);
  J[0](3,1)=-((gmmo*v[3])*iRh);
  J[0](4,1)=-2*v[4]+0.5*(2*v[4]-2*gmmo*v[4]);
  J[0](5,1)=0.5*(2*v[5]-2*gmmo*v[5]);
  J[0](6,1)=0.5*(2*v[6]-2*gmmo*v[6]);
  J[0](7,1)=0.5*gmmo;
  J[0](8,1)=0;
  J[0](9,1)=0;
  J[0](10,1)=0;

  J[0](0,2)=-((v[1]*v[2])*iRh2);
  J[0](1,2)=v[2]*iRh;
  J[0](2,2)=v[1]*iRh;
  J[0](3,2)=0;
  J[0](4,2)=-v[5];
  J[0](5,2)=-v[4];
  J[0](6,2)=0;
  J[0](7,2)=0;
  J[0](8,2)=0;
  J[0](9,2)=0;
  J[0](10,2)=0;

  J[0](0,3)=-((v[1]*v[3])*iRh2);
  J[0](1,3)=v[3]*iRh;
  J[0](2,3)=0;
  J[0](3,3)=v[1]*iRh;
  J[0](4,3)=-v[6];
  J[0](5,3)=0;
  J[0](6,3)=-v[4];
  J[0](7,3)=0;
  J[0](8,3)=0;
  J[0](9,3)=0;
  J[0](10,3)=0;

  J[0](0,4)=0;
  J[0](1,4)=0;
  J[0](2,4)=0;
  J[0](3,4)=0;
  J[0](4,4)=0;
  J[0](5,4)=0;
  J[0](6,4)=0;
  J[0](7,4)=0;
  J[0](8,4)=0;
  J[0](9,4)=0;
  J[0](10,4)=0;

  J[0](0,5)=Ec3*iRh;
  J[0](1,5)=v[5]*iRh;
  J[0](2,5)=-(v[4]*iRh);
  J[0](3,5)=0;
  J[0](4,5)=-(v[2]*iRh);
  J[0](5,5)=v[1]*iRh;
  J[0](6,5)=0;
  J[0](7,5)=0;
  J[0](8,5)=0;
  J[0](9,5)=0;
  J[0](10,5)=-eta;

  J[0](0,6)=-Ec2*iRh;
  J[0](1,6)=v[6]*iRh;
  J[0](2,6)=0;
  J[0](3,6)=-(v[4]*iRh);
  J[0](4,6)=-(v[3]*iRh);
  J[0](5,6)=0;
  J[0](6,6)=v[1]*iRh;
  J[0](7,6)=0;
  J[0](8,6)=0;
  J[0](9,6)=eta;
  J[0](10,6)=0;

  J[0](0,7)=2*iRh*(v[5]*Ec3-v[6]*Ec2)+v[1]*gmmo*Uk*iRh2-(v[1]*(Uk+p))*iRh2;
  J[0](1,7)=2*(v[5]*v[5]*iRh+v[6]*v[6]*iRh)+(v[1]*((2*v[1])*iRh-(2*gamma*v[1])*iRh))*iRh+(Uk+p)*iRh;
  J[0](2,7)=-(2*v[4]*v[5]*iRh)-(v[1]*2*gmmo*v[2])*iRh2;
  J[0](3,7)=-(2*v[4]*v[6]*iRh)-(v[1]*2*gmmo*v[3])*iRh2;
  J[0](4,7)=-(2*gamma*v[4]*v[1]*iRh)+2*(-((v[5]*v[2])*iRh)-(v[6]*v[3])*iRh);
  J[0](5,7)=-(2*gamma*v[5]*v[1]*iRh)+2*((v[5]*v[1])*iRh-E3);
  J[0](6,7)=-(2*gamma*v[6]*v[1]*iRh)+2*((v[6]*v[1])*iRh+E2);
  J[0](7,7)=(gamma*v[1])*iRh;
  J[0](8,7)=0;
  J[0](9,7)=2*eta*v[6];
  J[0](10,7)=-2*eta*v[5];

  J[0](0,8)=0;
  J[0](1,8)=0;
  J[0](2,8)=0;
  J[0](3,8)=0;
  J[0](4,8)=0;
  J[0](5,8)=0;
  J[0](6,8)=0;
  J[0](7,8)=0;
  J[0](8,8)=0;
  J[0](9,8)=0;
  J[0](10,8)=0;

  J[0](0,9)=0;
  J[0](1,9)=0;
  J[0](2,9)=0;
  J[0](3,9)=0;
  J[0](4,9)=0;
  J[0](5,9)=0;
  J[0](6,9)=1;
  J[0](7,9)=0;
  J[0](8,9)=0;
  J[0](9,9)=0;
  J[0](10,9)=0;

  J[0](0,10)=0;
  J[0](1,10)=0;
  J[0](2,10)=0;
  J[0](3,10)=0;
  J[0](4,10)=0;
  J[0](5,10)=-1;
  J[0](6,10)=0;
  J[0](7,10)=0;
  J[0](8,10)=0;
  J[0](9,10)=0;
  J[0](10,10)=0;

  
  
  J[1](0,0)=0;
  J[1](1,0)=0;
  J[1](2,0)=1;
  J[1](3,0)=0;
  J[1](4,0)=0;
  J[1](5,0)=0;
  J[1](6,0)=0;
  J[1](7,0)=0;
  J[1](8,0)=0;
  J[1](9,0)=0;
  J[1](10,0)=0;

  J[1](0,1)=-((v[1]*v[2])*iRh2);
  J[1](1,1)=v[2]*iRh;
  J[1](2,1)=v[1]*iRh;
  J[1](3,1)=0;
  J[1](4,1)=-v[5];
  J[1](5,1)=-v[4];
  J[1](6,1)=0;
  J[1](7,1)=0;
  J[1](8,1)=0;
  J[1](9,1)=0;
  J[1](10,1)=0;

  J[1](0,2)=-(v[2]*v[2]*iRh2)+(gmmo*Uk)*.5*iRh;
  J[1](1,2)=-((gmmo*v[1])*iRh);
  J[1](2,2)=(2*v[2])*iRh-(gmmo*v[2])*iRh;
  J[1](3,2)=-((gmmo*v[3])*iRh);
  J[1](4,2)=0.5*(2*v[4]-2*gmmo*v[4]);
  J[1](5,2)=-2*v[5]+0.5*(2*v[5]-2*gmmo*v[5]);
  J[1](6,2)=0.5*(2*v[6]-2*gmmo*v[6]);
  J[1](7,2)=0.5*gmmo;
  J[1](8,2)=0;
  J[1](9,2)=0;
  J[1](10,2)=0;

  J[1](0,3)=-((v[2]*v[3])*iRh2);
  J[1](1,3)=0;
  J[1](2,3)=v[3]*iRh;
  J[1](3,3)=v[2]*iRh;
  J[1](4,3)=0;
  J[1](5,3)=-v[6];
  J[1](6,3)=-v[5];
  J[1](7,3)=0;
  J[1](8,3)=0;
  J[1](9,3)=0;
  J[1](10,3)=0;

  J[1](0,4)=-(Ec3*iRh);
  J[1](1,4)=-(v[5]*iRh);
  J[1](2,4)=v[4]*iRh;
  J[1](3,4)=0;
  J[1](4,4)=v[2]*iRh;
  J[1](5,4)=-(v[1]*iRh);
  J[1](6,4)=0;
  J[1](7,4)=0;
  J[1](8,4)=0;
  J[1](9,4)=0;
  J[1](10,4)=eta;

  J[1](0,5)=0;
  J[1](1,5)=0;
  J[1](2,5)=0;
  J[1](3,5)=0;
  J[1](4,5)=0;
  J[1](5,5)=0;
  J[1](6,5)=0;
  J[1](7,5)=0;
  J[1](8,5)=0;
  J[1](9,5)=0;
  J[1](10,5)=0;

  J[1](0,6)=Ec1*iRh;
  J[1](1,6)=0;
  J[1](2,6)=v[6]*iRh;
  J[1](3,6)=-(v[5]*iRh);
  J[1](4,6)=0;
  J[1](5,6)=-(v[3]*iRh);
  J[1](6,6)=v[2]*iRh;
  J[1](7,6)=0;
  J[1](8,6)=-eta;
  J[1](9,6)=0;
  J[1](10,6)=0;

  J[1](0,7)=2*iRh*(-(v[4]*Ec3)+v[6]*Ec1)+v[2]*gmmo*Uk*iRh2-(v[2]*(Uk+p))*iRh2;
  J[1](1,7)=-(2*v[4]*v[5]*iRh)-2*gmmo*v[1]*v[2]*iRh2;
  J[1](2,7)=2*(v[4]*v[4]*iRh+v[6]*v[6]*iRh)+(v[2]*((2*v[2])*iRh-(2*gamma*v[2])*iRh))*iRh+(Uk+p)*iRh;
  J[1](3,7)=-(2*v[5]*v[6]*iRh)-2*gmmo*v[2]*v[3]*iRh2;
  J[1](4,7)=-(2*gamma*v[4]*v[2]*iRh)+2*((v[4]*v[2])*iRh+E3);
  J[1](5,7)=-(2*gamma*v[5]*v[2]*iRh)+2*(-((v[4]*v[1])*iRh)-(v[6]*v[3])*iRh);
  J[1](6,7)=-(2*gamma*v[6]*v[2]*iRh)+2*((v[6]*v[2])*iRh-E1);
  J[1](7,7)=(gamma*v[2])*iRh;
  J[1](8,7)=-2*eta*v[6];
  J[1](9,7)=0;
  J[1](10,7)=2*eta*v[4];

  J[1](0,8)=0;
  J[1](1,8)=0;
  J[1](2,8)=0;
  J[1](3,8)=0;
  J[1](4,8)=0;
  J[1](5,8)=0;
  J[1](6,8)=-1;
  J[1](7,8)=0;
  J[1](8,8)=0;
  J[1](9,8)=0;
  J[1](10,8)=0;

  J[1](0,9)=0;
  J[1](1,9)=0;
  J[1](2,9)=0;
  J[1](3,9)=0;
  J[1](4,9)=0;
  J[1](5,9)=0;
  J[1](6,9)=0;
  J[1](7,9)=0;
  J[1](8,9)=0;
  J[1](9,9)=0;
  J[1](10,9)=0;

  J[1](0,10)=0;
  J[1](1,10)=0;
  J[1](2,10)=0;
  J[1](3,10)=0;
  J[1](4,10)=1;
  J[1](5,10)=0;
  J[1](6,10)=0;
  J[1](7,10)=0;
  J[1](8,10)=0;
  J[1](9,10)=0;
  J[1](10,10)=0;

  
  
  J[2](0,0)=0;
  J[2](1,0)=0;
  J[2](2,0)=0;
  J[2](3,0)=1;
  J[2](4,0)=0;
  J[2](5,0)=0;
  J[2](6,0)=0;
  J[2](7,0)=0;
  J[2](8,0)=0;
  J[2](9,0)=0;
  J[2](10,0)=0;

  J[2](0,1)=-((v[1]*v[3])*iRh2);
  J[2](1,1)=v[3]*iRh;
  J[2](2,1)=0;
  J[2](3,1)=v[1]*iRh;
  J[2](4,1)=-v[6];
  J[2](5,1)=0;
  J[2](6,1)=-v[4];
  J[2](7,1)=0;
  J[2](8,1)=0;
  J[2](9,1)=0;
  J[2](10,1)=0;

  J[2](0,2)=-((v[2]*v[3])*iRh2);
  J[2](1,2)=0;
  J[2](2,2)=v[3]*iRh;
  J[2](3,2)=v[2]*iRh;
  J[2](4,2)=0;
  J[2](5,2)=-v[6];
  J[2](6,2)=-v[5];
  J[2](7,2)=0;
  J[2](8,2)=0;
  J[2](9,2)=0;
  J[2](10,2)=0;

  J[2](0,3)=-(v[3]*v[3]*iRh2)+(gmmo*Uk)*.5*iRh;
  J[2](1,3)=-((gmmo*v[1])*iRh);
  J[2](2,3)=-((gmmo*v[2])*iRh);
  J[2](3,3)=(2*v[3])*iRh-(gmmo*v[3])*iRh;
  J[2](4,3)=0.5*(2*v[4]-2*gmmo*v[4]);
  J[2](5,3)=0.5*(2*v[5]-2*gmmo*v[5]);
  J[2](6,3)=-2*v[6]+0.5*(2*v[6]-2*gmmo*v[6]);
  J[2](7,3)=0.5*gmmo;
  J[2](8,3)=0;
  J[2](9,3)=0;
  J[2](10,3)=0;

  J[2](0,4)=Ec2*iRh;
  J[2](1,4)=-(v[6]*iRh);
  J[2](2,4)=0;
  J[2](3,4)=v[4]*iRh;
  J[2](4,4)=v[3]*iRh;
  J[2](5,4)=0;
  J[2](6,4)=-(v[1]*iRh);
  J[2](7,4)=0;
  J[2](8,4)=0;
  J[2](9,4)=-eta;
  J[2](10,4)=0;

  J[2](0,5)=-(Ec1*iRh);
  J[2](1,5)=0;
  J[2](2,5)=-(v[6]*iRh);
  J[2](3,5)=v[5]*iRh;
  J[2](4,5)=0;
  J[2](5,5)=v[3]*iRh;
  J[2](6,5)=-(v[2]*iRh);
  J[2](7,5)=0;
  J[2](8,5)=eta;
  J[2](9,5)=0;
  J[2](10,5)=0;

  J[2](0,6)=0;
  J[2](1,6)=0;
  J[2](2,6)=0;
  J[2](3,6)=0;
  J[2](4,6)=0;
  J[2](5,6)=0;
  J[2](6,6)=0;
  J[2](7,6)=0;
  J[2](8,6)=0;
  J[2](9,6)=0;
  J[2](10,6)=0;

  J[2](0,7)=2*iRh*(v[4]*Ec2-v[5]*Ec1)+v[3]*gmmo*Uk*iRh2-(v[3]*(Uk+p))*iRh2;
  J[2](1,7)=-((2*v[4]*v[6])*iRh)-2*gmmo*v[1]*v[3])*iRh2;
  J[2](2,7)=-((2*v[5]*v[6])*iRh)-2*gmmo*v[2]*v[3])*iRh2;
  J[2](3,7)=2*(v[4]*v[4]*iRh+v[5]*v[5]*iRh)+(v[3]*((2*v[3])*iRh-(2*gamma*v[3])*iRh))*iRh+(Uk+p)*iRh;
  J[2](4,7)=-(2*gamma*v[4]*v[3]*iRh)+2*((v[4]*v[3])*iRh-E2);
  J[2](5,7)=-(2*gamma*v[5]*v[3]*iRh)+2*((v[5]*v[3])*iRh+E1);
  J[2](6,7)=2*(-((v[4]*v[1])*iRh)-(v[5]*v[2])*iRh)-(2*gamma*v[6]*v[3])*iRh;
  J[2](7,7)=(gamma*v[3])*iRh;
  J[2](8,7)=2*eta*v[5];
  J[2](9,7)=-2*eta*v[4];
  J[2](10,7)=0;

  J[2](0,8)=0;
  J[2](1,8)=0;
  J[2](2,8)=0;
  J[2](3,8)=0;
  J[2](4,8)=0;
  J[2](5,8)=1;
  J[2](6,8)=0;
  J[2](7,8)=0;
  J[2](8,8)=0;
  J[2](9,8)=0;
  J[2](10,8)=0;

  J[2](0,9)=0;
  J[2](1,9)=0;
  J[2](2,9)=0;
  J[2](3,9)=0;
  J[2](4,9)=-1;
  J[2](5,9)=0;
  J[2](6,9)=0;
  J[2](7,9)=0;
  J[2](8,9)=0;
  J[2](9,9)=0;
  J[2](10,9)=0;

  J[2](0,10)=0;
  J[2](1,10)=0;
  J[2](2,10)=0;
  J[2](3,10)=0;
  J[2](4,10)=0;
  J[2](5,10)=0;
  J[2](6,10)=0;
  J[2](7,10)=0;
  J[2](8,10)=0;
  J[2](9,10)=0;
  J[2](10,10)=0;

}
