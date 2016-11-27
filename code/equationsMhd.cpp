#include "util.h"
#include "equationsMhd.h"

template <int dim>
Equations<EquationsTypeMhd, dim>::Equations(Parameters<dim>& parameters) : parameters(parameters)
{
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::component_names()
{
  return{ "momentum", "momentum", "momentum", "density", "energy", "magnetic_field", "magnetic_field", "magnetic_field" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);

  return data_component_interpretation;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_kinetic_energy(const InputVector &W) const
{
  typename InputVector::value_type kinetic_energy = 0;

  for (unsigned int d = 0; d < dim; ++d)
    kinetic_energy += W[first_momentum_component + d] * W[first_momentum_component + d];

  kinetic_energy *= 1. / (2 * W[density_component]);

  return kinetic_energy;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W) const
{
  return ((this->parameters.gas_gamma - 1.0) * (W[energy_component] - compute_kinetic_energy(W)));
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > &flux) const
{
  const typename InputVector::value_type pressure = compute_pressure(W);

  for (unsigned int d = 0; d < dim; ++d)
  {
    for (unsigned int e = 0; e < dim; ++e)
      flux[first_momentum_component + d][e] = W[first_momentum_component + d] * W[first_momentum_component + e] / W[density_component];

    flux[first_momentum_component + d][d] += pressure;
  }

  for (unsigned int d = 0; d < dim; ++d)
    flux[density_component][d] = W[first_momentum_component + d];

  for (unsigned int d = 0; d < dim; ++d)
    flux[energy_component][d] = W[first_momentum_component + d] / W[density_component] * (W[energy_component] + pressure);
}

template <int dim>
template <typename InputVector, typename ValueType>
void Equations<EquationsTypeMhd, dim>::compute_jacobian_addition(double cell_diameter, const InputVector& grad_W, std_cxx11::array <std_cxx11::array <ValueType, dim>, n_components > &jacobian_addition) const
{
  for (unsigned int i = 0; i < n_components; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      jacobian_addition[i][d] = 0;
}


template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::Q(std_cxx11::array<typename InputVector::value_type, n_components> &result, const InputVector &W, const Tensor<1, dim> &normal) const
{
  std_cxx11::array<typename InputVector::value_type, n_components> forResult;
  for (unsigned int d = 0; d < n_components; d++)
    forResult[d] = 0;
  for (unsigned int d = 0; d < dim; d++)
  {
    forResult[first_momentum_component] += normal[d] * W[first_momentum_component + d];
    if (d > 0)
    {
      forResult[first_momentum_component + d] -= normal[d] * W[first_momentum_component];
      forResult[first_momentum_component + d] += normal[0] * W[first_momentum_component + d];
    }

    forResult[first_magnetic_flux_component] += normal[d] * W[first_magnetic_flux_component + d];
    if (d > 0)
    {
      forResult[first_magnetic_flux_component + d] -= normal[d] * W[first_magnetic_flux_component];
      forResult[first_magnetic_flux_component + d] += normal[0] * W[first_magnetic_flux_component + d];
    }
  }
  forResult[density_component] = W[density_component];
  forResult[energy_component] = W[energy_component];

  result[0] = W[density_component];
  for (unsigned int d = 0; d < dim; d++)
    result[d + 1] = forResult[d];
  for (unsigned int d = 0; d < dim; d++)
    result[d + first_momentum_component + dim + 1] = forResult[d];
  result[n_components - 1] = W[energy_component];
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::Q_inv(std_cxx11::array<typename InputVector::value_type, n_components> &result, std_cxx11::array<typename InputVector::value_type, n_components> &F, const Tensor<1, dim> &normal) const
{
  std_cxx11::array<typename InputVector::value_type, n_components> forResult;
  for (unsigned int d = 0; d < n_components; d++)
    forResult[d] = 0;
  for (unsigned int d = 0; d < dim; d++)
  {
    if (d == 0)
      forResult[first_momentum_component] += normal[d] * F[d + 1];
    else
    {
      forResult[first_momentum_component] -= normal[d] * F[d + 1];
      forResult[first_momentum_component + d] -= normal[d] * F[1];
      forResult[first_momentum_component + d] += normal[0] * F[d + 1];
    }

    if (d == 0)
      forResult[first_magnetic_flux_component] += normal[d] * F[dim + d + 1];
    else
    {
      forResult[first_magnetic_flux_component] -= normal[d] * F[dim + d + 1];
      forResult[first_magnetic_flux_component + d] -= normal[d] * F[dim + 1];
      forResult[first_magnetic_flux_component + d] += normal[0] * F[dim + d + 1];
    }
  }
  forResult[density_component] = F[0];
  forResult[energy_component] = F[n_components - 1];

  for (unsigned int d = 0; d < n_components; d++)
    result[d] = forResult[d];
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
  std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux) const
{
  typename InputVector::value_type srdl, srdr, hl[2], hr[2], Uk, Um, E2, E3, Sl, Sr, pml, pmr, B, B2, cl, cm, cr, ptl, ptr;
  typename InputVector::value_type sp[5], sml, smr, ptst, ptstr, vbstl, vbstr, Bsgnl, Bsgnr, invsumd;

  std_cxx11::array<typename InputVector::value_type, n_components> Fl;
  std_cxx11::array<typename InputVector::value_type, n_components> Fr;
  std_cxx11::array<typename InputVector::value_type, n_components> Udl;
  std_cxx11::array<typename InputVector::value_type, n_components> Udr;
  std_cxx11::array<typename InputVector::value_type, n_components> Ul;
  std_cxx11::array<typename InputVector::value_type, n_components> Ur;

  std_cxx11::array<typename InputVector::value_type, n_components> QedWplus;
  std_cxx11::array<typename InputVector::value_type, n_components> QedWminus;

  Q(QedWplus, Wplus, normal);
  Q(QedWminus, Wminus, normal);

  // Simple average of mag. field in direction of normal vector
  B = 0.5*(QedWplus[4] + QedWminus[4]);
  B2 = B*B;

  // Calculate left flux
  hl[0] = 1.0 / QedWminus[0];
  Uk = 0.5*hl[0] * (QedWminus[1] * QedWminus[1] + QedWminus[2] * QedWminus[2] + QedWminus[3] * QedWminus[3]);
  Um = 0.5*(B2 + QedWminus[5] * QedWminus[5] + QedWminus[6] * QedWminus[6]);
  hl[1] = (parameters.gas_gamma - 1)*(QedWminus[7] - Uk - Um);
  E2 = hl[0] * (QedWminus[1] * QedWminus[6] - QedWminus[3] * QedWminus[4]);
  E3 = hl[0] * (QedWminus[2] * QedWminus[4] - QedWminus[1] * QedWminus[5]);

  Fl[0] = QedWminus[1];
  Fl[1] = hl[0] * QedWminus[1] * QedWminus[1] - QedWminus[4] * QedWminus[4] + Um + hl[1];
  Fl[2] = hl[0] * QedWminus[1] * QedWminus[2] - QedWminus[4] * QedWminus[5];
  Fl[3] = hl[0] * QedWminus[1] * QedWminus[3] - QedWminus[4] * QedWminus[6];
  Fl[4] = 0.0;
  Fl[5] = -E3;
  Fl[6] = E2;
  Fl[7] = hl[0] * QedWminus[1] * (hl[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Uk) + (E2*QedWminus[6] - E3*QedWminus[5]);

  // Calculate right flux
  hr[0] = 1.0 / QedWplus[0];
  Uk = 0.5*hr[0] * (QedWplus[1] * QedWplus[1] + QedWplus[2] * QedWplus[2] + QedWplus[3] * QedWplus[3]);
  Um = 0.5*(B2 + QedWplus[5] * QedWplus[5] + QedWplus[6] * QedWplus[6]);
  hr[1] = (parameters.gas_gamma - 1)*(QedWplus[7] - Uk - Um);
  E2 = hr[0] * (QedWplus[1] * QedWplus[6] - QedWplus[3] * QedWplus[4]);
  E3 = hr[0] * (QedWplus[2] * QedWplus[4] - QedWplus[1] * QedWplus[5]);

  Fr[0] = QedWplus[1];
  Fr[1] = hr[0] * QedWplus[1] * QedWplus[1] - QedWplus[4] * QedWplus[4] + Um + hr[1];
  Fr[2] = hr[0] * QedWplus[1] * QedWplus[2] - QedWplus[4] * QedWplus[5];
  Fr[3] = hr[0] * QedWplus[1] * QedWplus[3] - QedWplus[4] * QedWplus[6];
  Fr[4] = 0.0;
  Fr[5] = -E3;
  Fr[6] = E2;
  Fr[7] = hr[0] * QedWplus[1] * (hr[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Uk) + (E2*QedWplus[6] - E3*QedWplus[5]);

  pml = 0.5*(QedWminus[4] * QedWminus[4] + QedWminus[5] * QedWminus[5] + QedWminus[6] * QedWminus[6]);
  pmr = 0.5*(QedWplus[4] * QedWplus[4] + QedWplus[5] * QedWplus[5] + QedWplus[6] * QedWplus[6]);
  // fast magnetoacoustic speed
  cl = parameters.gas_gamma*hl[1] + 2.0*pml;
  cl = sqrt(0.5*hl[0] * (cl + sqrt(cl*cl - 4.0*parameters.gas_gamma*hl[1] * QedWminus[4] * QedWminus[4])));
  cr = parameters.gas_gamma*hr[1] + 2.0*pmr;
  cr = sqrt(0.5*hr[0] * (cr + sqrt(cr*cr - 4.0*parameters.gas_gamma*hr[1] * QedWplus[4] * QedWplus[4])));

  // total pressure
  ptl = hl[1] + pml;
  ptr = hr[1] + pmr;

  // maximum of fast magnetoacoustic speeds L/R
  cm = (cl > cr) ? cl : cr;
  if (QedWminus[1] * hl[0] <= QedWplus[1] * hr[0])
  {
    sp[0] = QedWminus[1] * hl[0] - cm;
    sp[4] = QedWplus[1] * hr[0] + cm;
  }
  else
  {
    sp[0] = QedWplus[1] * hr[0] - cm;
    sp[4] = QedWminus[1] * hl[0] + cm;
  }

  // Upwind flux in the case of supersonic flow
  if (sp[0] >= 0.0)
  {
    // use F_L
    Q_inv<InputVector>(normal_flux, Fl, normal);
    return;
  }
  if (sp[4] <= 0.0)
  {
    // use F_R
    Q_inv<InputVector>(normal_flux, Fr, normal);
    return;
  }

  // Determine Alfven and middle speeds
  Sl = sp[0] - QedWminus[1] * hl[0];
  Sr = sp[4] - QedWplus[1] * hr[0];
  sp[2] = (QedWminus[1] * Sl - QedWplus[1] * Sr - ptl + ptr) / (QedWminus[0] * Sl - QedWplus[0] * Sr);
  sml = sp[0] - sp[2];
  smr = sp[4] - sp[2];

  // Density
  Ul[0] = QedWminus[0] * Sl / sml;
  Ur[0] = QedWplus[0] * Sr / smr;

  Ul[4] = Udl[4] = QedWminus[4];
  Ur[4] = Udr[4] = QedWplus[4];

  srdl = sqrt(Ul[0]);
  srdr = sqrt(Ur[0]);

  // Sl*
  sp[1] = sp[2] - fabs(Ul[4]) / srdl;
  // Sr*
  sp[3] = sp[2] + fabs(Ur[4]) / srdr;

  ptst = ptl + QedWminus[0] * Sl*(Sl - sml);
  ptstr = ptr + QedWplus[0] * Sr*(Sr - smr);

  // F*_L
  Ul[1] = Ul[0] * sp[2];

  cl = QedWminus[0] * Sl*sml - Ul[4] * Ul[4];
  if (fabs(cl) < 1e-8 * ptst)
  {
    Ul[2] = Ul[0] * QedWminus[2] * hl[0];
    Ul[3] = Ul[0] * QedWminus[3] * hl[0];
    Ul[5] = QedWminus[5];
    Ul[6] = QedWminus[6];
  }
  else
  {
    cl = 1.0 / cl;
    cm = Ul[4] * (Sl - sml)*cl;
    Ul[2] = Ul[0] * (QedWminus[2] * hl[0] - QedWminus[5] * cm);
    Ul[3] = Ul[0] * (QedWminus[3] * hl[0] - QedWminus[6] * cm);
    cm = (QedWminus[0] * Sl*Sl - Ul[4] * Ul[4])*cl;
    Ul[5] = QedWminus[5] * cm;
    Ul[6] = QedWminus[6] * cm;
  }
  vbstl = (Ul[1] * Ul[4] + Ul[2] * Ul[5] + Ul[3] * Ul[6]) / Ul[0];

  Ul[7] = (Sl*QedWminus[7] - ptl*QedWminus[1] * hl[0] + ptst*sp[2] + Ul[4] *
    ((QedWminus[1] * QedWminus[4] + QedWminus[2] * QedWminus[5] + QedWminus[3] * QedWminus[6])*hl[0] - vbstl)) / sml;

  // F*_R
  Ur[1] = Ur[0] * sp[2];
  cl = QedWplus[0] * Sr*smr - Ur[4] * Ur[4];
  if (fabs(cl) < 1e-8*ptstr)
  {
    Ur[2] = Ur[0] * QedWplus[2] * hr[0];
    Ur[3] = Ur[0] * QedWplus[3] * hr[0];
    Ur[5] = QedWplus[5];
    Ur[6] = QedWplus[6];
  }
  else
  {
    cl = 1.0 / cl;
    cm = Ur[4] * (Sr - smr)*cl;
    Ur[2] = Ur[0] * (QedWplus[2] * hr[0] - QedWplus[5] * cm);
    Ur[3] = Ur[0] * (QedWplus[3] * hr[0] - QedWplus[6] * cm);
    cm = (QedWplus[0] * Sr*Sr - Ur[4] * Ur[4])*cl;
    Ur[5] = QedWplus[5] * cm;
    Ur[6] = QedWplus[6] * cm;
  }
  vbstr = (Ur[1] * Ur[4] + Ur[2] * Ur[5] + Ur[3] * Ur[6]) / Ur[0];

  Ur[7] = (Sr*QedWplus[7] - ptr*QedWplus[1] * hr[0] + ptstr*sp[2] + Ur[4] *
    ((QedWplus[1] * Ur[4] + QedWplus[2] * QedWplus[5] + QedWplus[3] * QedWplus[6])*hr[0] - vbstr)) / smr;

  if (sp[1] >= 0.0)
  {
    for (register int j = 0; j < n_components; j++)
      normal_flux[j] = Fl[j] + sp[0] * (Ul[j] - QedWminus[j]);
    Q_inv<InputVector>(normal_flux, normal_flux, normal);
    return;
  }
  if (sp[3] <= 0.0 && sp[2] < 0.0)
  {
    for (register int j = 0; j < n_components; j++)
      normal_flux[j] = Fr[j] + sp[4] * (Ur[j] - QedWplus[j]);
    Q_inv<InputVector>(normal_flux, normal_flux, normal);
    return;
  }

  // F**_L and F**_R
  if (B2 < 1e-8*(ptst + ptstr))
  {
    for (register int j = 0; j < n_components; j++)
    {
      Udl[j] = Ul[j];
      Udr[j] = Ur[j];
    }
  }
  else
  {
    invsumd = 1.0 / (srdl + srdr);
    Bsgnl = (Ul[4] > 0.0) ? 1.0 : -1.0;
    Bsgnr = (Ur[4] > 0.0) ? 1.0 : -1.0;

    Udl[0] = Ul[0];
    Udr[0] = Ur[0];

    Udl[1] = Ul[1];
    Udr[1] = Ur[1];

    cm = invsumd*(srdl*Ul[2] / Ul[0] + srdr*Ur[2] / Ur[0]);
    cl = invsumd*(Ur[5] - Ul[5]);
    Udl[2] = Ul[0] * (cm + Bsgnl*cl);
    Udr[2] = Ur[0] * (cm + Bsgnr*cl);

    cm = invsumd*(srdl*Ul[3] / Ul[0] + srdr*Ur[3] / Ur[0]);
    cl = invsumd*(Ur[6] - Ul[6]);
    Udl[3] = Ul[0] * (cm + Bsgnl*cl);
    Udr[3] = Ur[0] * (cm + Bsgnr*cl);

    cm = invsumd*(srdl*Ur[5] + srdr*Ul[5]);
    cl = invsumd*srdl*srdr*(Ur[2] / Ur[0] - Ul[2] / Ul[0]);
    Udl[5] = cm + Bsgnl*cl;
    Udr[5] = cm + Bsgnr*cl;

    cm = invsumd*(srdl*Ur[6] + srdr*Ul[6]);
    cl = invsumd*srdl*srdr*(Ur[3] / Ur[0] - Ul[3] / Ul[0]);
    Udl[6] = cm + Bsgnl*cl;
    Udr[6] = cm + Bsgnr*cl;

    Udl[7] = Ul[7] - srdl*Bsgnl*(vbstl - sp[2] * Ul[4] - (Udl[2] * Udl[5] + Udl[3] * Udl[6]) / Udl[0]);
    Udr[7] = Ur[7] + srdr*Bsgnr*(vbstr - sp[2] * Ur[4] - (Udr[2] * Udr[5] + Udr[3] * Udr[6]) / Udr[0]);
  }

  if (sp[2] >= 0.0)
  {
    cm = sp[1] - sp[0];
    for (register int j = 0; j < n_components; j++)
      normal_flux[j] = Fl[j] + sp[1] * Udl[j] - sp[0] * QedWminus[j] - cm*Ul[j];
    Q_inv<InputVector>(normal_flux, normal_flux, normal);
    return;
  }
  cm = sp[3] - sp[4];
  for (register int j = 0; j < n_components; j++)
    normal_flux[j] = Fr[j] + sp[3] * Udr[j] - sp[4] * QedWplus[j] - cm*Ur[j];
  Q_inv<InputVector>(normal_flux, normal_flux, normal);
  return;
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing) const
{
  // LK: Tohle u me budou samy nuly
  const double gravity = -1.0;

  for (unsigned int c = 0; c < n_components; ++c)
  {
    switch (c)
    {
    case first_momentum_component + 1:
      forcing[c] = gravity * W[density_component];
      break;
    case energy_component:
      forcing[c] = gravity * W[first_momentum_component + 1];
      break;
    default:
      forcing[c] = 0;
    }
  }
}

template <int dim>
template <typename DataVector>
void Equations<EquationsTypeMhd, dim>::compute_Wminus(const BoundaryKind(&boundary_kind)[n_components],
  const Tensor<1, dim> &normal_vector,
  const DataVector     &Wplus,
  const Vector<double> &boundary_values,
  const DataVector     &Wminus) const
{
  for (unsigned int c = 0; c < n_components; c++)
  {
    switch (boundary_kind[c])
    {
    case inflow_boundary:
    {
      Wminus[c] = boundary_values(c);
      break;
    }

    case outflow_boundary:
    default:
    {
      Wminus[c] = Wplus[c];
      break;
    }
    }
  }
}

template <int dim>
Equations<EquationsTypeMhd, dim>::Postprocessor::Postprocessor(Equations<EquationsTypeMhd, dim>& equations) : equations(equations)
{}

template <int dim>
void
Equations<EquationsTypeMhd, dim>::Postprocessor::compute_derived_quantities_vector(
  const std::vector<Vector<double> > &uh,
  const std::vector<std::vector<Tensor<1, dim> > > &duh,
  const std::vector<std::vector<Tensor<2, dim> > > &/*dduh*/,
  const std::vector<Point<dim> > &/*normals*/,
  const std::vector<Point<dim> > &/*evaluation_points*/,
  std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size();

  Assert(computed_quantities.size() == n_quadrature_points,
    ExcInternalError());

  Assert(uh[0].size() == n_components,
    ExcInternalError());

  Assert(computed_quantities[0].size() == dim + 1, ExcInternalError());

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = uh[q](density_component);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = uh[q](first_momentum_component + d) / density;

    computed_quantities[q](dim) = equations.compute_pressure(uh[q]);
  }
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::Postprocessor::get_names() const
{
  return{ "velocity", "velocity", "velocity", "pressure" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);

  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  return interpretation;
}

template <int dim>
UpdateFlags Equations<EquationsTypeMhd, dim>::Postprocessor::get_needed_update_flags() const
{
  return update_values;
}

template class Equations<EquationsTypeMhd, 2>;
template class Equations<EquationsTypeMhd, 3>;

#include "equationsMhd_inst.cpp"