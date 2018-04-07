#include "numericalFlux.h"

template <EquationsType equationsType, int dim>
void NumFlux<equationsType, dim>::Q(component_vector &result, const component_vector &W, const Tensor<1, dim> &normal)
{
  double forResult[n_components];
  for (unsigned int d = 0; d < n_components; d++)
    forResult[d] = W[d];
  if (normal[0] > 0.5) { // nothing
  }
  else if (normal[0] < -0.5) {
    forResult[1] = -W[1];
    forResult[5] = -W[5];
    forResult[2] = -W[2];
    forResult[6] = -W[6];
  }
  else if (normal[1] > 0.5) {
    forResult[1] = W[2];
    forResult[2] = -W[1];
    forResult[5] = W[6];
    forResult[6] = -W[5];
  }
  else if (normal[1] < -0.5) {
    forResult[1] = -W[2];
    forResult[2] = W[1];
    forResult[5] = -W[6];
    forResult[6] = W[5];
  }
  else if (normal[2] > 0.5) {
    forResult[1] = W[3];
    forResult[3] = -W[1];
    forResult[5] = W[7];
    forResult[7] = -W[5];
  }
  else if (normal[2] < -0.5) {
    forResult[1] = -W[3];
    forResult[3] = W[1];
    forResult[5] = -W[7];
    forResult[7] = W[5];
  }

  for (unsigned int d = 0; d < n_components; d++)
    result[d] = forResult[d];
}

template <EquationsType equationsType, int dim>
void NumFlux<equationsType, dim>::Q_inv(component_vector &result, component_vector &W, const Tensor<1, dim> &normal)
{
  double forResult[n_components];
  for (unsigned int d = 0; d < n_components; d++)
    forResult[d] = W[d];
  if (normal[0] > 0.5) { // nothing
  }
  else if (normal[0] < -0.5) {
    forResult[1] = -W[1];
    forResult[5] = -W[5];
    forResult[2] = -W[2];
    forResult[6] = -W[6];
  }
  else if (normal[1] > 0.5) {
    forResult[1] = -W[2];
    forResult[2] = W[1];
    forResult[5] = -W[6];
    forResult[6] = W[5];
  }
  else if (normal[1] < -0.5) {
    forResult[1] = W[2];
    forResult[2] = -W[1];
    forResult[5] = W[6];
    forResult[6] = -W[5];
  }
  else if (normal[2] > 0.5) {
    forResult[1] = -W[3];
    forResult[3] = W[1];
    forResult[5] = -W[7];
    forResult[7] = W[5];
  }
  else if (normal[2] < -0.5) {
    forResult[1] = W[3];
    forResult[3] = -W[1];
    forResult[5] = W[7];
    forResult[7] = -W[5];
  }
  for (unsigned int d = 0; d < n_components; d++)
    result[d] = forResult[d];
}

template <EquationsType equationsType, int dim>
void NumFluxLaxFriedrich<equationsType, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const component_vector &Wplus_,
  const component_vector &Wminus_, component_vector &normal_flux, double& max_speed) const
{
  double hl[2], hr[2], spd[5];

  component_vector ul, ur;
  Q(ul, Wplus_, normal);
  Q(ur, Wminus_, normal);

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

  // maximum of fast magnetoacoustic speeds L/R
  double cm = (cl > cr) ? cl : cr;
  if (ul[1] * hl[0] <= ur[1] * hr[0]) {
    spd[0] = ul[1] * hl[0] - cm;
    spd[4] = ur[1] * hr[0] + cm;
  }
  else {
    spd[0] = ur[1] * hr[0] - cm;
    spd[4] = ul[1] * hl[0] + cm;
  }

  max_speed = std::max(max_speed, (std::max(std::abs(spd[0]), std::abs(spd[4]))));

  std::array<std::array <component_vector::value_type, 3>, n_components > iflux, oflux;

  Equations<equationsType, dim>::compute_flux_matrix(Wplus_, iflux, this->parameters);
  Equations<equationsType, dim>::compute_flux_matrix(Wminus_, oflux, this->parameters);

  for (unsigned int di = 0; di < n_components; ++di)
  {
    normal_flux[di] = 0.;
    for (unsigned int d = 0; d < dim; ++d)
      normal_flux[di] += 0.5 * (iflux[di][d] + oflux[di][d]) * normal[d];

    normal_flux[di] += this->parameters.lax_friedrich_stabilization_value * (Wplus_[di] - Wminus_[di]);
  }
}

template <EquationsType equationsType, int dim>
void NumFluxHLLD<equationsType, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const component_vector &Wplus_,
  const component_vector &Wminus_, component_vector &normal_flux, double& max_speed) const
{
  component_vector flux_lf;
  if (parameters.debug)
  {
    NumFluxLaxFriedrich<equationsType, dim> lf(this->parameters);
    lf.numerical_normal_flux(normal, Wplus_, Wminus_, flux_lf, max_speed);
  }

  double Fl[n_components], Fr[n_components], hl[2], hr[2];
  double Uldst[n_components], Urdst[n_components], Ulst[n_components], Urst[n_components];
  double spd[5], vbstl, vbstr, Bsgnl, Bsgnr, invsumd;

  component_vector ul, ur;
  Q(ul, Wplus_, normal);
  Q(ur, Wminus_, normal);

  double B = 0.5*(ul[5] + ur[5]);
  double B2 = B*B;

  // Densities, energies.
  hl[0] = 1.0 / ul[0];
  double Ukl = 0.5 * hl[0] * (ul[1] * ul[1] + ul[2] * ul[2] + ul[3] * ul[3]);
  double Uml = 0.5 * (ul[5] * ul[5] + ul[6] * ul[6] + ul[7] * ul[7]);
  hl[1] = (parameters.gas_gamma - 1.) * (ul[4] - Ukl - Uml);

  hr[0] = 1.0 / ur[0];
  double Ukr = 0.5 * hr[0] * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]);
  double Umr = 0.5 * (ur[5] * ur[5] + ur[6] * ur[6] + ur[7] * ur[7]);
  hr[1] = (parameters.gas_gamma - 1.) * (ur[4] - Ukr - Umr);

  // sound speed
  double al2 = parameters.gas_gamma * hl[1] * hl[0];
  double ar2 = parameters.gas_gamma * hr[1] * hr[0];

  // fast magnetoacoustic speed
  double cl = al2 + (2. * Uml * hl[0]);
  cl = sqrt(0.5 * (cl + sqrt((cl * cl) - (4.0 * al2 * ul[5] * ul[5] * hl[0]))));

  double cr = ar2 + (2. * Umr * hr[0]);
  cr = sqrt(0.5 * (cr + sqrt((cr * cr) - (4.0 * ar2 * ur[5] * ur[5] * hr[0]))));

  // total pressure
  double ptl = hl[1] + Uml;
  double ptr = hr[1] + Umr;

  // maximum of fast magnetoacoustic speeds L/R
  double cm = (cl > cr) ? cl : cr;
  if (ul[1] * hl[0] <= ur[1] * hr[0]) {
    spd[0] = ul[1] * hl[0] - cm;
    spd[4] = ur[1] * hr[0] + cm;
  }
  else {
    spd[0] = ur[1] * hr[0] - cm;
    spd[4] = ul[1] * hl[0] + cm;
  }

  max_speed = std::max(max_speed, (std::max(std::abs(spd[0]), std::abs(spd[4]))));

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
  Fl[4] = (ul[4] + ptl) * (ul[1] / ul[0]) - (ul[5] * (ul[1] * ul[5] + ul[2] * ul[6] + ul[3] * ul[7]));

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
  Fr[4] = (ur[4] + ptr) * (ur[1] / ur[0]) - (ur[5] * (ur[1] * ur[5] + ur[2] * ur[6] + ur[3] * ur[7]));

  // Upwind flux in the case of supersonic flow
  if (spd[0] >= 0.0) {
    for (int j = 0; j < n_components; j++)
      normal_flux[j] = Fl[j];
    Q_inv(normal_flux, normal_flux, normal);
    if (parameters.debug)
      for (int j = 0; j < n_components; j++)
      {
        if ((std::abs(flux_lf[j]) > 1e-10) && (std::abs(normal_flux[j]) > 1e-10))
          if (std::abs(flux_lf[j] - normal_flux[j]) > std::min(std::abs(flux_lf[j]), std::abs(normal_flux[j])))
          {
            std::cout << "n: " << normal << ", component: " << j << ", L-F: " << flux_lf[j] << ", result: " << normal_flux[j] << std::endl;
          }
      }
    return;
  }
  if (spd[4] <= 0.0) {
    for (int j = 0; j < n_components; j++)
      normal_flux[j] = Fr[j];
    Q_inv(normal_flux, normal_flux, normal);
    if (parameters.debug)
      for (int j = 0; j < n_components; j++)
      {
        if ((std::abs(flux_lf[j]) > 1e-10) && (std::abs(normal_flux[j]) > 1e-10))
          if (std::abs(flux_lf[j] - normal_flux[j]) > std::min(std::abs(flux_lf[j]), std::abs(normal_flux[j])))
          {
            std::cout << "n: " << normal << ", component: " << j << ", L-F: " << flux_lf[j] << ", result: " << normal_flux[j] << std::endl;
          }
      }
    return;
  }

  // Determine Alfven and middle speeds
  double sdl = spd[0] - ul[1] * hl[0];
  double sdr = spd[4] - ur[1] * hr[0];
  spd[2] = (ur[1] * sdr - ul[1] * sdl - ptr + ptl) / (ur[0] * sdr - ul[0] * sdl);
  double sdml = spd[0] - spd[2];
  double sdmr = spd[4] - spd[2];

  Ulst[0] = ul[0] * sdl / sdml;
  Urst[0] = ur[0] * sdr / sdmr;

  Ulst[5] = Uldst[5] = ul[5];
  Urst[5] = Urdst[5] = ur[5];

  double sqrtdl = sqrt(Ulst[0]);
  double sqrtdr = sqrt(Urst[0]);

  // Sl*, Sr*
  spd[1] = spd[2] - fabs(Ulst[5]) / sqrtdl;
  spd[3] = spd[2] + fabs(Urst[5]) / sqrtdr;

  double ptst = ptl + ul[0] * sdl*(sdl - sdml);
  double ptstr = ptr + ur[0] * sdr*(sdr - sdmr);

  // F*_L
  Ulst[1] = Ulst[0] * spd[2];
  cl = ul[0] * sdl * sdml - Ulst[5] * Ulst[5];
  if (fabs(cl) < NEGLIGIBLE * ptst) {
    Ulst[2] = Ulst[0] * ul[2] * hl[0];
    Ulst[3] = Ulst[0] * ul[3] * hl[0];

    Ulst[6] = ul[6];
    Ulst[7] = ul[7];
  }
  else {
    cl = 1.0 / cl;
    cm = Ulst[5] * (sdl - sdml) * cl;
    Ulst[2] = Ulst[0] * (ul[2] * hl[0] - ul[6] * cm);
    Ulst[3] = Ulst[0] * (ul[3] * hl[0] - ul[7] * cm);
    cm = (ul[0] * sdl * sdl - Ulst[5] * Ulst[5]) * cl;
    Ulst[6] = ul[6] * cm;
    Ulst[7] = ul[7] * cm;
  }
  vbstl = (Ulst[1] * Ulst[5] + Ulst[2] * Ulst[6] + Ulst[3] * Ulst[7]) / Ulst[0];
  Ulst[4] = (sdl * ul[4] - ptl * ul[1] * hl[0] + ptst * spd[2] + Ulst[5] *
    ((ul[1] * ul[5] + ul[2] * ul[6] + ul[3] * ul[7]) * hl[0] - vbstl)) / sdml;

  // F*_R
  Urst[1] = Urst[0] * spd[2];
  cl = ur[0] * sdr * sdmr - Urst[5] * Urst[5];
  if (fabs(cl) < NEGLIGIBLE * ptstr) {
    Urst[2] = Urst[0] * ur[2] * hr[0];
    Urst[3] = Urst[0] * ur[3] * hr[0];

    Urst[6] = ur[6];
    Urst[7] = ur[7];
  }
  else {
    cl = 1.0 / cl;
    cm = Urst[5] * (sdr - sdmr) * cl;
    Urst[2] = Urst[0] * (ur[2] * hr[0] - ur[6] * cm);
    Urst[3] = Urst[0] * (ur[3] * hr[0] - ur[7] * cm);
    cm = (ur[0] * sdr * sdr - Urst[5] * Urst[5])*cl;
    Urst[6] = ur[6] * cm;
    Urst[7] = ur[7] * cm;
  }
  vbstr = (Urst[1] * Urst[5] + Urst[2] * Urst[6] + Urst[3] * Urst[7]) / Urst[0];
  Urst[4] = (sdr * ur[4] - ptr * ur[1] * hr[0] + ptstr * spd[2] + Urst[5] *
    ((ur[1] * Urst[5] + ur[2] * ur[6] + ur[3] * ur[7]) * hr[0] - vbstr)) / sdmr;

  if (spd[1] >= 0.0) {
    for (int j = 0; j < n_components; j++)
      normal_flux[j] = Fl[j] + spd[0] * (Ulst[j] - ul[j]);
    Q_inv(normal_flux, normal_flux, normal);
    if (parameters.debug)
      for (int j = 0; j < n_components; j++)
      {
        if ((std::abs(flux_lf[j]) > 1e-10) && (std::abs(normal_flux[j]) > 1e-10))
          if (std::abs(flux_lf[j] - normal_flux[j]) > std::min(std::abs(flux_lf[j]), std::abs(normal_flux[j])))
          {
            std::cout << "n: " << normal << ", component: " << j << ", L-F: " << flux_lf[j] << ", result: " << normal_flux[j] << std::endl;
          }
      }
    return;
  }
  if (spd[3] <= 0.0 && spd[2] < 0.0) {
    for (int j = 0; j < n_components; j++)
      normal_flux[j] = Fr[j] + spd[4] * (Urst[j] - ur[j]);
    Q_inv(normal_flux, normal_flux, normal);
    if (parameters.debug)
      for (int j = 0; j < n_components; j++)
      {
        if ((std::abs(flux_lf[j]) > 1e-10) && (std::abs(normal_flux[j]) > 1e-10))
          if (std::abs(flux_lf[j] - normal_flux[j]) > std::min(std::abs(flux_lf[j]), std::abs(normal_flux[j])))
          {
            std::cout << "n: " << normal << ", component: " << j << ", L-F: " << flux_lf[j] << ", result: " << normal_flux[j] << std::endl;
          }
      }
    return;
  }

  // F**_L and F**_R
  if (B2 < NEGLIGIBLE * (ptst + ptstr)) {
    for (int j = 0; j < n_components; j++) {
      Uldst[j] = Ulst[j];
      Urdst[j] = Urst[j];
    }
  }
  else {
    invsumd = 1.0 / (sqrtdl + sqrtdr);
    Bsgnl = (Ulst[5] > 0.0) ? 1.0 : -1.0;
    Bsgnr = (Urst[5] > 0.0) ? 1.0 : -1.0;

    Uldst[0] = Ulst[0];
    Urdst[0] = Urst[0];

    Uldst[1] = Ulst[1];
    Urdst[1] = Urst[1];

    cm = invsumd * (sqrtdl * Ulst[2] / Ulst[0] + sqrtdr * Urst[2] / Urst[0]);
    cl = invsumd * (Urst[6] - Ulst[6]);
    Uldst[2] = Ulst[0] * (cm + Bsgnl*cl);
    Urdst[2] = Urst[0] * (cm + Bsgnr*cl);

    cm = invsumd * (sqrtdl * Ulst[3] / Ulst[0] + sqrtdr * Urst[3] / Urst[0]);
    cl = invsumd * (Urst[7] - Ulst[7]);
    Uldst[3] = Ulst[0] * (cm + Bsgnl*cl);
    Urdst[3] = Urst[0] * (cm + Bsgnr*cl);

    cm = invsumd * (sqrtdl *Urst[6] + sqrtdr * Ulst[6]);
    cl = invsumd * sqrtdl * sqrtdr * (Urst[2] / Urst[0] - Ulst[2] / Ulst[0]);
    Uldst[6] = cm + Bsgnl * cl;
    Urdst[6] = cm + Bsgnr * cl;

    cm = invsumd * (sqrtdl * Urst[7] + sqrtdr * Ulst[7]);
    cl = invsumd * sqrtdl * sqrtdr * (Urst[3] / Urst[0] - Ulst[3] / Ulst[0]);
    Uldst[7] = cm + Bsgnl * cl;
    Urdst[7] = cm + Bsgnr * cl;

    Uldst[4] = Ulst[4] - sqrtdl * Bsgnl * (vbstl - spd[2] * Ulst[5] - (Uldst[2] * Uldst[6] + Uldst[3] * Uldst[7]) / Uldst[0]);
    Urdst[4] = Urst[4] + sqrtdr * Bsgnr * (vbstr - spd[2] * Urst[5] - (Urdst[2] * Urdst[6] + Urdst[3] * Urdst[7]) / Urdst[0]);
  }

  if (spd[2] >= 0.0) {
    cm = spd[1] - spd[0];
    for (int j = 0; j < n_components; j++)
      normal_flux[j] = Fl[j] + spd[1] * Uldst[j] - spd[0] * ul[j] - cm*Ulst[j];
    Q_inv(normal_flux, normal_flux, normal);
    if (parameters.debug)
      for (int j = 0; j < n_components; j++)
      {
        if ((std::abs(flux_lf[j]) > 1e-10) && (std::abs(normal_flux[j]) > 1e-10))
          if (std::abs(flux_lf[j] - normal_flux[j]) > std::min(std::abs(flux_lf[j]), std::abs(normal_flux[j])))
          {
            std::cout << "n: " << normal << ", component: " << j << ", L-F: " << flux_lf[j] << ", result: " << normal_flux[j] << std::endl;
          }
      }
    return;
  }
  else
  {
    cm = spd[3] - spd[4];
    for (int j = 0; j < n_components; j++)
      normal_flux[j] = Fr[j] + spd[3] * Urdst[j] - spd[4] * ur[j] - cm*Urst[j];
    Q_inv(normal_flux, normal_flux, normal);
    if (parameters.debug)
      for (int j = 0; j < n_components; j++)
      {
        if ((std::abs(flux_lf[j]) > 1e-10) && (std::abs(normal_flux[j]) > 1e-10))
          if (std::abs(flux_lf[j] - normal_flux[j]) > std::min(std::abs(flux_lf[j]), std::abs(normal_flux[j])))
          {
            std::cout << "n: " << normal << ", component: " << j << ", L-F: " << flux_lf[j] << ", result: " << normal_flux[j] << std::endl;
          }
      }
    return;
  }
  exit(1);
}

template class NumFlux<EquationsTypeMhd, 3>;
template class NumFluxLaxFriedrich<EquationsTypeMhd, 3>;
template class NumFluxHLLD<EquationsTypeMhd, 3>;
