#include "util.h"
#include "equationsMhd.h"

template <int dim>
Equations<EquationsTypeMhd, dim>::Equations(Parameters<dim>& parameters) : parameters(parameters)
{
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::component_names()
{
  return{ "density", "momentum", "momentum", "momentum", "energy", "magnetic_field", "magnetic_field", "magnetic_field" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  return data_component_interpretation;
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_kinetic_energy(const InputVector &W)
{
  return 0.5 * (W[1] * W[1] + W[2] * W[2] + W[3] * W[3]) / W[0];
}

template <>
double Equations<EquationsTypeMhd, 3>::compute_magnetic_energy(const InputVector &W)
{
  return 0.5 * (W[5] * W[5] + W[6] * W[6] + W[7] * W[7]);
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W) const
{
  return std::max(0., (this->parameters.gas_gamma - 1.0) * (W[4] - compute_kinetic_energy(W) - compute_magnetic_energy(W)));
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_kinetic_energy(const std_cxx11::array<double, n_components> &W)
{
  return 0.5 * (W[1] * W[1] + W[2] * W[2] + W[3] * W[3]) / W[0];
}

template <>
double Equations<EquationsTypeMhd, 3>::compute_magnetic_energy(const std_cxx11::array<double, n_components> &W)
{
  return 0.5 * (W[5] * W[5] + W[6] * W[6] + W[7] * W[7]);
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_pressure(const std_cxx11::array<double, n_components> &W) const
{
  return std::max(0., (this->parameters.gas_gamma - 1.0) * (W[4] - compute_kinetic_energy(W) - compute_magnetic_energy(W)));
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_total_pressure(const InputVector &W) const
{
  const double magnetic_energy = compute_magnetic_energy(W);
  return compute_pressure(W, compute_kinetic_energy(W), magnetic_energy) + magnetic_energy;
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W, const double& Uk, const double& Um) const
{
  return (this->parameters.gas_gamma - 1.0) * (W[4] - Uk - Um);
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_magnetic_field_divergence(const std::vector<Tensor<1, dim> > &W) const
{
  double divergence = 0.;

  for (unsigned int d = 0; d < dim; ++d)
    divergence += W[d + 5][d];

  return divergence;
}

template <int dim>
void Equations<EquationsTypeMhd, dim>::store_max_signal_speed(typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1>::value_type val)
{
  this->max_signal_speed = val.val();
}

template <int dim>
void Equations<EquationsTypeMhd, dim>::store_max_signal_speed(double val)
{
  this->max_signal_speed = val;
}

template <>
void Equations<EquationsTypeMhd, 3>::Q(std_cxx11::array<double, n_components> &result, const InputVector &W, const Tensor<1, 3> &normal) const
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

template <int dim>
void Equations<EquationsTypeMhd, dim>::Q_inv(std_cxx11::array<double, n_components> &result, std_cxx11::array<double, n_components> &W, const Tensor<1, dim> &normal) const
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

template <int dim>
void Equations<EquationsTypeMhd, dim>::compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <double, dim>, n_components > &flux) const
{
  const double mag_energy = compute_magnetic_energy(W);
  const double pressure = compute_pressure(W);
  const double E = W[4] + mag_energy;
  const double total_pressure = pressure + mag_energy;
  const double UB = W[1] * W[5] + W[2] * W[6] + W[3] * W[7];

  flux[0][0] = W[1];
  flux[1][0] = (W[1] * W[1] / W[0]) - W[5] * W[5] + total_pressure;
  flux[2][0] = (W[1] * W[2] / W[0]) - W[5] * W[6];
  flux[3][0] = (W[1] * W[3] / W[0]) - W[5] * W[7];
  flux[4][0] = (E + total_pressure) * (W[1] / W[0]) - (W[5] * UB);
  flux[5][0] = 0.0;
  flux[6][0] = ((W[1] / W[0]) * W[6]) - ((W[2] / W[0]) * W[5]);
  flux[7][0] = ((W[1] / W[0]) * W[7]) - ((W[3] / W[0]) * W[5]);

  flux[0][1] = W[2];
  flux[1][1] = (W[2] * W[1] / W[0]) - W[6] * W[5];
  flux[2][1] = (W[2] * W[2] / W[0]) - W[6] * W[6] + total_pressure;
  flux[3][1] = (W[2] * W[3] / W[0]) - W[6] * W[7];
  flux[4][1] = (E + total_pressure) * (W[2] / W[0]) - (W[6] * UB);
  flux[5][1] = ((W[2] / W[0]) * W[5]) - ((W[1] / W[0]) * W[6]);
  flux[6][1] = 0.0;
  flux[7][1] = ((W[2] / W[0]) * W[7]) - ((W[3] / W[0]) * W[6]);

  flux[0][2] = W[3];
  flux[1][2] = (W[3] * W[1] / W[0]) - W[7] * W[5];
  flux[2][2] = (W[3] * W[2] / W[0]) - W[7] * W[6];
  flux[3][2] = (W[3] * W[3] / W[0]) - W[7] * W[7] + total_pressure;
  flux[4][2] = (E + total_pressure) * (W[3] / W[0]) - (W[7] * UB);
  flux[5][2] = ((W[3] / W[0]) * W[5]) - ((W[1] / W[0]) * W[7]);
  flux[6][2] = ((W[3] / W[0]) * W[6]) - ((W[2] / W[0]) * W[7]);
  flux[7][2] = 0.0;
}

template <int dim>
void Equations<EquationsTypeMhd, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus_, const InputVector &Wminus_,
  std_cxx11::array<double, n_components> &normal_flux)
{
  double srdl, srdr, Fl[n_components], Fr[n_components], hl[2], hr[2], Uk, Um, E2, E3, Sl, Sr, pml, pmr, B, B2;
  double Udl[n_components], Udr[n_components], Ul[n_components], Ur[n_components], cl, cm, cr, ptl, ptr;
  double sp[5], sml, smr, ptst, ptstr, vbstl, vbstr, Bsgnl, Bsgnr, invsumd;

  std_cxx11::array<double, n_components> ul, ur;
  Q(ul, Wplus_, normal);
  Q(ur, Wminus_, normal);

  B = 0.5*(ul[5] + ur[5]);  // Simple average of mag. field in direction of normal vector
  B2 = B*B;

  // Calculate left flux
  hl[0] = 1.0 / ul[0];
  Uk = 0.5*hl[0] * (ul[1] * ul[1] + ul[2] * ul[2] + ul[3] * ul[3]);
  Um = 0.5*(ul[5] * ul[5] + ul[6] * ul[6] + ul[7] * ul[7]);       // ifdef _Pb_ then use B2 instead of ul[5]^2
  hl[1] = (parameters.gas_gamma - 1)*(ul[4] - Uk - Um);
  E2 = hl[0] * (ul[1] * ul[7] - ul[3] * ul[5]);
  E3 = hl[0] * (ul[2] * ul[5] - ul[1] * ul[6]);

  Fl[0] = ul[1];
  Fl[1] = hl[0] * ul[1] * ul[1] - ul[5] * ul[5] + Um + hl[1];
  Fl[2] = hl[0] * ul[1] * ul[2] - ul[5] * ul[6];
  Fl[3] = hl[0] * ul[1] * ul[3] - ul[5] * ul[7];
  Fl[5] = 0.0;
  Fl[6] = -E3;
  Fl[7] = E2;
  Fl[4] = hl[0] * ul[1] * (hl[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Uk) + (E2*ul[7] - E3*ul[6]);

  // Calculate right flux
  hr[0] = 1.0 / ur[0];
  Uk = 0.5*hr[0] * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]);
  Um = 0.5*(ur[5] * ur[5] + ur[6] * ur[6] + ur[7] * ur[7]);       // ifdef _Pb_ then use B2 instead of ur[5]^2
  hr[1] = (parameters.gas_gamma - 1)*(ur[4] - Uk - Um);
  E2 = hr[0] * (ur[1] * ur[7] - ur[3] * ur[5]);
  E3 = hr[0] * (ur[2] * ur[5] - ur[1] * ur[6]);

  Fr[0] = ur[1];
  Fr[1] = hr[0] * ur[1] * ur[1] - ur[5] * ur[5] + Um + hr[1];
  Fr[2] = hr[0] * ur[1] * ur[2] - ur[5] * ur[6];
  Fr[3] = hr[0] * ur[1] * ur[3] - ur[5] * ur[7];
  Fr[5] = 0.0;
  Fr[6] = -E3;
  Fr[7] = E2;
  Fr[4] = hr[0] * ur[1] * (hr[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Uk) + (E2*ur[7] - E3*ur[6]);

  pml = 0.5*(ul[5] * ul[5] + ul[6] * ul[6] + ul[7] * ul[7]);
  pmr = 0.5*(ur[5] * ur[5] + ur[6] * ur[6] + ur[7] * ur[7]);
  // fast magnetoacoustic speed
  cl = parameters.gas_gamma*hl[1] + 2.0*pml;
  cl = sqrt(0.5*hl[0] * (cl + sqrt(cl*cl - 4.0*parameters.gas_gamma*hl[1] * ul[5] * ul[5])));
  cr = parameters.gas_gamma*hr[1] + 2.0*pmr;
  cr = sqrt(0.5*hr[0] * (cr + sqrt(cr*cr - 4.0*parameters.gas_gamma*hr[1] * ur[5] * ur[5])));

  ptl = hl[1] + pml;  // total pressure
  ptr = hr[1] + pmr;

  // maximum of fast magnetoacoustic speeds L/R
  cm = (cl>cr) ? cl : cr;
  if (ul[1] * hl[0] <= ur[1] * hr[0]) {
    sp[0] = ul[1] * hl[0] - cm;
    sp[4] = ur[1] * hr[0] + cm;
  }
  else {
    sp[0] = ur[1] * hr[0] - cm;
    sp[4] = ul[1] * hl[0] + cm;
  }

  this->store_max_signal_speed(std::max(std::abs(sp[0]), std::abs(sp[4])));

  // Upwind flux in the case of supersonic flow
  if (sp[0] >= 0.0) {  // use F_L
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
  Sl = sp[0] - ul[1] * hl[0];
  Sr = sp[4] - ur[1] * hr[0];
  sp[2] = (ul[1] * Sl - ur[1] * Sr - ptl + ptr) / (ul[0] * Sl - ur[0] * Sr);
  sml = sp[0] - sp[2];
  smr = sp[4] - sp[2];

  Ul[0] = ul[0] * Sl / sml;  // Density
  Ur[0] = ur[0] * Sr / smr;

  Ul[5] = Udl[5] = ul[5];
  Ur[5] = Udr[5] = ur[5];
  Ul[5] = Ur[5] = Udl[5] = Udr[5] = B; // Magnetic field Bx (normal direction)
  srdl = sqrt(Ul[0]);
  srdr = sqrt(Ur[0]);

  sp[1] = sp[2] - fabs(Ul[5]) / srdl;  // Sl*
  sp[3] = sp[2] + fabs(Ur[5]) / srdr;  // Sr*

  ptst = ptl + ul[0] * Sl*(Sl - sml);
  ptstr = ptr + ur[0] * Sr*(Sr - smr);

  // F*_L
  Ul[1] = Ul[0] * sp[2];

  cl = ul[0] * Sl*sml - Ul[5] * Ul[5];
  if (fabs(cl)<NEGLIGIBLE*ptst) {
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
  if (fabs(cl)<NEGLIGIBLE*ptstr) {
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
  if (sp[3] <= 0.0 && sp[2]<0.0) {
    for (int j = 0; j < 8; j++)
      normal_flux[j] = Fr[j] + sp[4] * (Ur[j] - ur[j]);
    Q_inv(normal_flux, normal_flux, normal);
    return;
  }

  // F**_L and F**_R
  if (B2<NEGLIGIBLE*(ptst + ptstr)) {
    for (int j = 0; j < 8; j++) {
      Udl[j] = Ul[j];
      Udr[j] = Ur[j];
    }
  }
  else {
    invsumd = 1.0 / (srdl + srdr);
    Bsgnl = (Ul[5]>0.0) ? 1.0 : -1.0;
    Bsgnr = (Ur[5]>0.0) ? 1.0 : -1.0;

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

template <int dim>
Equations<EquationsTypeMhd, dim>::Postprocessor::Postprocessor(Equations<EquationsTypeMhd, dim>& equations) : equations(equations)
{}

#if DEAL_II_VERSION_MAJOR > 8 || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR > 5) || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR == 5 && DEAL_II_VERSION_SUBMINOR > 0)
template <int dim>
void
Equations<EquationsTypeMhd, dim>::Postprocessor::evaluate_vector_field(
  const DataPostprocessorInputs::Vector<dim> &inputs,
  std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = inputs.solution_values.size();

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = inputs.solution_values[q](0);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = inputs.solution_values[q](1 + d) / density;

    computed_quantities[q](dim) = equations.compute_pressure(inputs.solution_values[q]);
    computed_quantities[q](dim + 1) = equations.compute_magnetic_field_divergence(inputs.solution_gradients[q]);
  }
}
#else
template <int dim>
void
Equations<EquationsTypeMhd, dim>::Postprocessor::compute_derived_quantities_vector(
  const std::vector<Vector<double> > &uh,
  const std::vector<std::vector<Tensor<1, dim> > > &duh,
  const std::vector<std::vector<Tensor<2, dim> > > &,
  const std::vector<Point<dim> > &,
  const std::vector<Point<dim> > &,
  std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size();

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = uh[q](0);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = uh[q](1 + d) / density;

    std_cxx11::array<double, n_components> uh_q;
    for (unsigned int d = 0; d < n_components; ++d)
      uh_q[d] = uh[q][d];

    computed_quantities[q](dim) = equations.compute_pressure(uh_q);
    computed_quantities[q](dim + 1) = equations.compute_magnetic_field_divergence(duh[q]);
  }
}
#endif

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::Postprocessor::get_names() const
{
  return{ "velocity", "velocity", "velocity", "pressure", "mag_field_divergence" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  return interpretation;
}

template <int dim>
UpdateFlags Equations<EquationsTypeMhd, dim>::Postprocessor::get_needed_update_flags() const
{
  return update_values | update_gradients;
}


template class Equations<EquationsTypeMhd, 3>;
