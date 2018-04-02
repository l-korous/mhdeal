#include "numericalFlux.h"

template <EquationsType equationsType, int dim>
void NumFlux<equationsType, dim>::Q(std::array<double, Equations<equationsType, dim>::n_components> &result, const std::array<double, Equations<equationsType, dim>::n_components> &W, const Tensor<1, dim> &normal)
{
  double forResult[8];
  double b = asin(normal[2]);
  double cb = cos(b);
  double a;
  if (std::abs(normal[0]) < NEGLIGIBLE)
  {
    if (std::abs(cb) < NEGLIGIBLE)
      a = 0.;
    else
      a = asin(normal[1] / cb);
  }
  else
    a = acos(normal[0] / cb);

  double sa = sin(a);
  double sb = normal[2];
  double ca = cos(a);

  forResult[0] = W[0];

  forResult[1]
    = (ca * cb * W[1])
    + (sa * cb * W[2])
    + (sb * W[3]);

  forResult[2]
    = (-sa * W[1])
    + (ca * W[2]);

  forResult[3]
    = (-ca * sb * W[1])
    - (sa * sb * W[2])
    + (cb * W[3]);

  forResult[5]
    = (ca * cb * W[5])
    + (sa * cb * W[6])
    + (sb * W[7]);

  forResult[6]
    = (-sa * W[5])
    + (ca * W[6]);

  forResult[7]
    = (-ca * sb * W[5])
    - (sa * sb * W[6])
    + (cb * W[7]);

  forResult[4] = W[4];

  for (unsigned int d = 0; d < 8; d++)
    result[d] = forResult[d];
}

template <EquationsType equationsType, int dim>
void NumFlux<equationsType, dim>::Q_inv(std::array<double, Equations<equationsType, dim>::n_components> &result, std::array<double, Equations<equationsType, dim>::n_components> &W, const Tensor<1, dim> &normal)
{
  double forResult[8];
  double b = asin(normal[2]);
  double cb = cos(b);
  double a;
  if (std::abs(normal[0]) < NEGLIGIBLE)
  {
    if (std::abs(cb) < NEGLIGIBLE)
      a = 0.;
    else
      a = asin(normal[1] / cb);
  }
  else
    a = acos(normal[0] / cb);

  double sa = sin(a);
  double sb = normal[2];
  double ca = cos(a);

  forResult[0] = W[0];

  forResult[1]
    = (ca * cb * W[1])
    - (sa * W[2])
    - (ca * sb * W[3]);

  forResult[1 + 1]
    = (sa * cb * W[1])
    + (ca * W[1 + 1])
    - (sa * sb * W[3]);

  forResult[3]
    = sb * W[1]
    + cb * W[3];

  forResult[5]
    = (ca * cb * W[5])
    - (sa * W[6])
    - (ca * sb * W[7]);

  forResult[6]
    = (sa * cb * W[5])
    + (ca * W[6])
    - (sa * sb * W[7]);

  forResult[7]
    = sb * W[5]
    + cb * W[7];

  forResult[4] = W[4];

  for (unsigned int d = 0; d < 8; d++)
    result[d] = forResult[d];
}

template <EquationsType equationsType, int dim>
void NumFluxLaxFriedrich<equationsType, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const std::array<double, Equations<equationsType, dim>::n_components> &Wplus_, 
  const std::array<double, Equations<equationsType, dim>::n_components> &Wminus_, std::array<double, Equations<equationsType, dim>::n_components> &normal_flux, double& max_speed) const
{
  double Fl[Equations<equationsType, dim>::n_components], Fr[Equations<equationsType, dim>::n_components], hl[2], hr[2];
  double Udl[Equations<equationsType, dim>::n_components], Udr[Equations<equationsType, dim>::n_components], Ul[Equations<equationsType, dim>::n_components], Ur[Equations<equationsType, dim>::n_components];
  double sp[5], vbstl, vbstr, Bsgnl, Bsgnr, invsumd;

  std::array<double, Equations<equationsType, dim>::n_components> ul, ur;
  Q(ul, Wplus_, normal);
  Q(ur, Wminus_, normal);

  double B = 0.5*(ul[5] + ur[5]);
  double B2 = B*B;

  // Densities, energies.
  hl[0] = 1.0 / ul[0];
  double Ukl = 0.5 * hl[0] * (ul[1] * ul[1] + ul[2] * ul[2] + ul[3] * ul[3]);
  double Uml = 0.5 * (ul[5] * ul[5] + ul[6] * ul[6] + ul[7] * ul[7]);
  hl[1] = (parameters.gas_gamma - 1) * (ul[4] - Ukl - Uml);

  hr[0] = 1.0 / ur[0];
  double Ukr = 0.5 * hr[0] * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]);
  double Umr = 0.5 * (ur[5] * ur[5] + ur[6] * ur[6] + ur[7] * ur[7]);
  hr[1] = (parameters.gas_gamma - 1) * (ur[4] - Ukr - Umr);

  // sound speed
  double al2 = parameters.gas_gamma * hl[1] * hl[0];
  double ar2 = parameters.gas_gamma * hr[1] * hr[0];

  // fast magnetoacoustic speed
  double cl = (al2 + (2. * Uml * hl[0]));
  cl = sqrt(0.5 * (cl + sqrt((cl * cl) - (4.0 * al2 * ul[5] * ul[5] * hl[0]))));

  double cr = ar2 + (2. * Umr * hr[0]);
  cr = sqrt(0.5 * (cr + sqrt((cr * cr) - (4.0 * ar2 * ur[5] * ur[5] * hr[0]))));

  // total pressure
  double ptl = hl[1] + Uml;
  double ptr = hr[1] + Umr;

  // maximum of fast magnetoacoustic speeds L/R
  double cm = (cl > cr) ? cl : cr;
  if (ul[1] * hl[0] <= ur[1] * hr[0]) {
    sp[0] = ul[1] * hl[0] - cm;
    sp[4] = ur[1] * hr[0] + cm;
  }
  else {
    sp[0] = ur[1] * hr[0] - cm;
    sp[4] = ul[1] * hl[0] + cm;
  }

  max_speed = std::max(max_speed, (std::max(std::abs(sp[0]), std::abs(sp[4]))));

  std::array<std::array <std::array<double, Equations<equationsType, dim>::n_components>::value_type, 3>, Equations<equationsType, dim>::n_components > iflux, oflux;

  Equations<equationsType, dim>::compute_flux_matrix(Wplus_, iflux, this->parameters);
  Equations<equationsType, dim>::compute_flux_matrix(Wminus_, oflux, this->parameters);

  for (unsigned int di = 0; di < Equations<equationsType, dim>::n_components; ++di)
  {
    normal_flux[di] = 0.;
    for (unsigned int d = 0; d < dim; ++d)
      normal_flux[di] += 0.5 * (iflux[di][d] + oflux[di][d]) * normal[d];

    normal_flux[di] += this->parameters.lax_friedrich_stabilization_value * (Wplus_[di] - Wminus_[di]);
  }
}

template <EquationsType equationsType, int dim>
void NumFluxHLLD<equationsType, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const std::array<double, Equations<equationsType, dim>::n_components> &Wplus_,
  const std::array<double, Equations<equationsType, dim>::n_components> &Wminus_, std::array<double, Equations<equationsType, dim>::n_components> &normal_flux, double& max_speed) const
{
  double Fl[Equations<equationsType, dim>::n_components], Fr[Equations<equationsType, dim>::n_components], hl[2], hr[2];
  double Udl[Equations<equationsType, dim>::n_components], Udr[Equations<equationsType, dim>::n_components], Ul[Equations<equationsType, dim>::n_components], Ur[Equations<equationsType, dim>::n_components];
  double sp[5], vbstl, vbstr, Bsgnl, Bsgnr, invsumd;

  std::array<double, Equations<equationsType, dim>::n_components> ul, ur;
  Q(ul, Wplus_, normal);
  Q(ur, Wminus_, normal);

  double B = 0.5*(ul[5] + ur[5]);
  double B2 = B*B;

  // Densities, energies.
  hl[0] = 1.0 / ul[0];
  double Ukl = 0.5 * hl[0] * (ul[1] * ul[1] + ul[2] * ul[2] + ul[3] * ul[3]);
  double Uml = 0.5 * (ul[5] * ul[5] + ul[6] * ul[6] + ul[7] * ul[7]);
  hl[1] = (parameters.gas_gamma - 1) * (ul[4] - Ukl - Uml);

  hr[0] = 1.0 / ur[0];
  double Ukr = 0.5 * hr[0] * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]);
  double Umr = 0.5 * (ur[5] * ur[5] + ur[6] * ur[6] + ur[7] * ur[7]);
  hr[1] = (parameters.gas_gamma - 1) * (ur[4] - Ukr - Umr);

  // sound speed
  double al2 = parameters.gas_gamma * hl[1] * hl[0];
  double ar2 = parameters.gas_gamma * hr[1] * hr[0];

  // fast magnetoacoustic speed
  double cl = (al2 + (2. * Uml * hl[0]));
  cl = sqrt(0.5 * (cl + sqrt((cl * cl) - (4.0 * al2 * ul[5] * ul[5] * hl[0]))));

  double cr = ar2 + (2. * Umr * hr[0]);
  cr = sqrt(0.5 * (cr + sqrt((cr * cr) - (4.0 * ar2 * ur[5] * ur[5] * hr[0]))));

  // total pressure
  double ptl = hl[1] + Uml;
  double ptr = hr[1] + Umr;

  // maximum of fast magnetoacoustic speeds L/R
  double cm = (cl > cr) ? cl : cr;
  if (ul[1] * hl[0] <= ur[1] * hr[0]) {
    sp[0] = ul[1] * hl[0] - cm;
    sp[4] = ur[1] * hr[0] + cm;
  }
  else {
    sp[0] = ur[1] * hr[0] - cm;
    sp[4] = ul[1] * hl[0] + cm;
  }

  max_speed = std::max(max_speed, (std::max(std::abs(sp[0]), std::abs(sp[4]))));

  // Calculate left flux
  double E2 = hl[0] * (ul[1] * ul[7] - ul[3] * ul[5]);
  double E3 = hl[0] * (ul[2] * ul[5] - ul[1] * ul[6]);

  Fl[0] = ul[1];
  Fl[1] = hl[0] * ul[1] * ul[1] - ul[5] * ul[5] + ptl;
  Fl[2] = hl[0] * ul[1] * ul[2] - ul[5] * ul[6];
  Fl[3] = hl[0] * ul[1] * ul[3] - ul[5] * ul[7];
  Fl[5] = 0.0;
  Fl[6] = -E3;
  Fl[7] = E2;
  Fl[4] = hl[0] * ul[1] * (hl[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Ukl) + (E2 * ul[7] - E3 * ul[6]);

  // Calculate right flux
  E2 = hr[0] * (ur[1] * ur[7] - ur[3] * ur[5]);
  E3 = hr[0] * (ur[2] * ur[5] - ur[1] * ur[6]);

  Fr[0] = ur[1];
  Fr[1] = hr[0] * ur[1] * ur[1] - ur[5] * ur[5] + ptr;
  Fr[2] = hr[0] * ur[1] * ur[2] - ur[5] * ur[6];
  Fr[3] = hr[0] * ur[1] * ur[3] - ur[5] * ur[7];
  Fr[5] = 0.0;
  Fr[6] = -E3;
  Fr[7] = E2;
  Fr[4] = hr[0] * ur[1] * (hr[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Ukr) + (E2*ur[7] - E3*ur[6]);

  // Upwind flux in the case of supersonic flow
  if (sp[0] >= 0.0) {
    for (int j = 0; j < 8; j++)
      normal_flux[j] = Fl[j];
    Q_inv(normal_flux, normal_flux, normal);
    return;
  }
  if (sp[4] <= 0.0) {  // use F_R
    for (int j = 0; j < 8; j++)
      normal_flux[j] = Fr[j];
    Q_inv(normal_flux, normal_flux, normal);
    return;
  }

  // Determine Alfven and middle speeds
  double Sl = sp[0] - ul[1] * hl[0];
  double Sr = sp[4] - ur[1] * hr[0];
  sp[2] = (ul[1] * Sl - ur[1] * Sr - ptl + ptr) / (ul[0] * Sl - ur[0] * Sr);
  double sml = sp[0] - sp[2];
  double smr = sp[4] - sp[2];

  Ul[0] = ul[0] * Sl / sml;  // Density
  Ur[0] = ur[0] * Sr / smr;

  Ul[5] = Udl[5] = ul[5];
  Ur[5] = Udr[5] = ur[5];
  Ul[5] = Ur[5] = Udl[5] = Udr[5] = B;
  double srdl = sqrt(Ul[0]);
  double srdr = sqrt(Ur[0]);

  sp[1] = sp[2] - fabs(Ul[5]) / srdl;  // Sl*
  sp[3] = sp[2] + fabs(Ur[5]) / srdr;  // Sr*

  double ptst = ptl + ul[0] * Sl*(Sl - sml);
  double ptstr = ptr + ur[0] * Sr*(Sr - smr);

  // F*_L
  Ul[1] = Ul[0] * sp[2];

  cl = ul[0] * Sl*sml - Ul[5] * Ul[5];
  if (fabs(cl) < NEGLIGIBLE*ptst) {
    Ul[2] = Ul[0] * ul[2] * hl[0];
    Ul[3] = Ul[0] * ul[3] * hl[0];

    Ul[6] = ul[6];
    Ul[7] = ul[7];
  }
  else {
    cl = 1.0 / cl;
    cm = Ul[5] * (Sl - sml)*cl;
    Ul[2] = Ul[0] * (ul[2] * hl[0] - ul[6] * cm);
    Ul[3] = Ul[0] * (ul[3] * hl[0] - ul[7] * cm);
    cm = (ul[0] * Sl*Sl - Ul[5] * Ul[5])*cl;
    Ul[6] = ul[6] * cm;
    Ul[7] = ul[7] * cm;
  }
  vbstl = (Ul[1] * Ul[5] + Ul[2] * Ul[6] + Ul[3] * Ul[7]) / Ul[0];

  Ul[4] = (Sl*ul[4] - ptl*ul[1] * hl[0] + ptst*sp[2] + Ul[5] *
    ((ul[1] * ul[5] + ul[2] * ul[6] + ul[3] * ul[7])*hl[0] - vbstl)) / sml;

  // F*_R
  Ur[1] = Ur[0] * sp[2];
  cl = ur[0] * Sr*smr - Ur[5] * Ur[5];
  if (fabs(cl) < NEGLIGIBLE*ptstr) {
    Ur[2] = Ur[0] * ur[2] * hr[0];
    Ur[3] = Ur[0] * ur[3] * hr[0];

    Ur[6] = ur[6];
    Ur[7] = ur[7];
  }
  else {
    cl = 1.0 / cl;
    cm = Ur[5] * (Sr - smr)*cl;
    Ur[2] = Ur[0] * (ur[2] * hr[0] - ur[6] * cm);
    Ur[3] = Ur[0] * (ur[3] * hr[0] - ur[7] * cm);
    cm = (ur[0] * Sr*Sr - Ur[5] * Ur[5])*cl;
    Ur[6] = ur[6] * cm;
    Ur[7] = ur[7] * cm;
  }
  vbstr = (Ur[1] * Ur[5] + Ur[2] * Ur[6] + Ur[3] * Ur[7]) / Ur[0];

  Ur[4] = (Sr*ur[4] - ptr*ur[1] * hr[0] + ptstr*sp[2] + Ur[5] *
    ((ur[1] * Ur[5] + ur[2] * ur[6] + ur[3] * ur[7])*hr[0] - vbstr)) / smr;

  if (sp[1] >= 0.0) {
    for (int j = 0; j < 8; j++)
      normal_flux[j] = Fl[j] + sp[0] * (Ul[j] - ul[j]);
    Q_inv(normal_flux, normal_flux, normal);
    return;
  }
  if (sp[3] <= 0.0 && sp[2] < 0.0) {
    for (int j = 0; j < 8; j++)
      normal_flux[j] = Fr[j] + sp[4] * (Ur[j] - ur[j]);
    Q_inv(normal_flux, normal_flux, normal);
    return;
  }

  // F**_L and F**_R
  if (B2 < NEGLIGIBLE*(ptst + ptstr)) {
    for (int j = 0; j < 8; j++) {
      Udl[j] = Ul[j];
      Udr[j] = Ur[j];
    }
  }
  else {
    invsumd = 1.0 / (srdl + srdr);
    Bsgnl = (Ul[5] > 0.0) ? 1.0 : -1.0;
    Bsgnr = (Ur[5] > 0.0) ? 1.0 : -1.0;

    Udl[0] = Ul[0];
    Udr[0] = Ur[0];

    Udl[1] = Ul[1];
    Udr[1] = Ur[1];

    cm = invsumd*(srdl*Ul[2] / Ul[0] + srdr*Ur[2] / Ur[0]);
    cl = invsumd*(Ur[6] - Ul[6]);
    Udl[2] = Ul[0] * (cm + Bsgnl*cl);
    Udr[2] = Ur[0] * (cm + Bsgnr*cl);

    cm = invsumd*(srdl*Ul[3] / Ul[0] + srdr*Ur[3] / Ur[0]);
    cl = invsumd*(Ur[7] - Ul[7]);
    Udl[3] = Ul[0] * (cm + Bsgnl*cl);
    Udr[3] = Ur[0] * (cm + Bsgnr*cl);

    cm = invsumd*(srdl*Ur[6] + srdr*Ul[6]);
    cl = invsumd*srdl*srdr*(Ur[2] / Ur[0] - Ul[2] / Ul[0]);
    Udl[6] = cm + Bsgnl*cl;
    Udr[6] = cm + Bsgnr*cl;

    cm = invsumd*(srdl*Ur[7] + srdr*Ul[7]);
    cl = invsumd*srdl*srdr*(Ur[3] / Ur[0] - Ul[3] / Ul[0]);
    Udl[7] = cm + Bsgnl*cl;
    Udr[7] = cm + Bsgnr*cl;

    Udl[4] = Ul[4] - srdl*Bsgnl*(vbstl - sp[2] * Ul[5] - (Udl[2] * Udl[6] + Udl[3] * Udl[7]) / Udl[0]);
    Udr[4] = Ur[4] + srdr*Bsgnr*(vbstr - sp[2] * Ur[5] - (Udr[2] * Udr[6] + Udr[3] * Udr[7]) / Udr[0]);
  }

  if (sp[2] >= 0.0) {
    cm = sp[1] - sp[0];
    for (int j = 0; j < 8; j++)
      normal_flux[j] = Fl[j] + sp[1] * Udl[j] - sp[0] * ul[j] - cm*Ul[j];
    Q_inv(normal_flux, normal_flux, normal);
    return;
  }
  else
  {
    cm = sp[3] - sp[4];
    for (int j = 0; j < 8; j++)
      normal_flux[j] = Fr[j] + sp[3] * Udr[j] - sp[4] * ur[j] - cm*Ur[j];
    Q_inv(normal_flux, normal_flux, normal);
    return;
  }
}

template class NumFlux<EquationsTypeMhd, 3>;
template class NumFluxLaxFriedrich<EquationsTypeMhd, 3>;
template class NumFluxHLLD<EquationsTypeMhd, 3>;
