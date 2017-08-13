#include "util.h"
#include "equationsMhd.h"

template <int dim>
Equations<EquationsTypeMhd, dim>::Equations(Parameters<dim>& parameters) : parameters(parameters)
{
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::component_names()
{
  return{ "density", "momentum", "momentum", "momentum", "magnetic_field", "magnetic_field", "magnetic_field", "energy" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  return data_component_interpretation;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_kinetic_energy(const InputVector &W)
{
  typename InputVector::value_type kinetic_energy = 0.;

  for (unsigned int d = 0; d < dim; ++d)
    kinetic_energy += W[first_momentum_component + d] * W[first_momentum_component + d];

  kinetic_energy = (0.5 * kinetic_energy) / W[density_component];

  return kinetic_energy;
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_magnetic_field_divergence(const std::vector<Tensor<1, dim> > &W) const
{
  double divergence = 0.;

  for (unsigned int d = 0; d < dim; ++d)
    divergence += W[first_magnetic_flux_component + d][d];

  return divergence;
}

template <>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, 3>::compute_magnetic_energy(const InputVector &W)
{
  return 0.5 *
    (
      W[first_magnetic_flux_component] * W[first_magnetic_flux_component]
      + W[first_magnetic_flux_component + 1] * W[first_magnetic_flux_component + 1]
      + W[first_magnetic_flux_component + 2] * W[first_magnetic_flux_component + 2]
      );
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W) const
{
  return (this->parameters.gas_gamma - 1.0) * (W[energy_component] - compute_kinetic_energy(W) - compute_magnetic_energy(W));
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W, const typename InputVector::value_type& Uk, const typename InputVector::value_type& Um) const
{
  return (this->parameters.gas_gamma - 1.0) * (W[energy_component] - Uk - Um);
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_flux_vector(const int derivative, const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components > &flux) const
{
  std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > flux_matrix;
  compute_flux_matrix(W, flux_matrix);
  for (unsigned int d = 0; d < n_components; ++d)
    flux[d] = flux_matrix[d][derivative];
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > &flux) const
{
  const typename InputVector::value_type pressure = compute_pressure(W);
  typename InputVector::value_type kinetic_energy = compute_kinetic_energy(W);
  typename InputVector::value_type magnetic_energy = compute_magnetic_energy(W);
  typename InputVector::value_type energy_item[dim];

  // The result variable "flux" has two dimensions [][] - the first is the component (density, momentum x, ...), the second is the spatial derivative.
  // That is why the loop through the second dimension has always as many steps as there are spatial dimensions in the physical space.

  for (unsigned int d = 0; d < dim; ++d)
    flux[density_component][d] = W[first_momentum_component + d];

  for (unsigned int d = 0; d < dim; ++d)
  {
    for (unsigned int e = 0; e < dim; ++e)
    {
      flux[first_momentum_component + d][e] = W[first_momentum_component + d] * W[first_momentum_component + e] / W[density_component];
      flux[first_momentum_component + d][e] -= W[first_magnetic_flux_component + d] * W[first_magnetic_flux_component + e];
    }

    flux[first_momentum_component + d][d] += pressure;
    flux[first_momentum_component + d][d] += magnetic_energy;
  }

  for (unsigned int d = 0; d < dim; ++d)
  {
    for (unsigned int e = 0; e < dim; ++e)
    {
      if (d == e)
        flux[first_magnetic_flux_component + d][e] = 0.;
      else
      {
        flux[first_magnetic_flux_component + d][e] = W[first_momentum_component + e] * W[first_magnetic_flux_component + d] / W[density_component];
        flux[first_magnetic_flux_component + d][e] -= W[first_magnetic_flux_component + e] * W[first_momentum_component + d] / W[density_component];
      }
    }
  }

  for (unsigned int d = 0; d < dim; ++d)
  {
    flux[energy_component][d] = (W[energy_component] + pressure + magnetic_energy) * W[first_momentum_component + d] / W[density_component];
    for (unsigned int e = 0; e < dim; ++e)
      flux[energy_component][d] -= (W[first_magnetic_flux_component + d] * W[first_momentum_component + e] * W[first_magnetic_flux_component + e]) / W[density_component];
  }
}

template <int dim>
template <typename InputVector, typename ValueType>
void Equations<EquationsTypeMhd, dim>::compute_jacobian_addition(double cell_diameter, const InputVector& grad_W, std_cxx11::array <std_cxx11::array <ValueType, dim>, n_components > &jacobian_addition) const
{
  for (unsigned int i = 0; i < n_components; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      jacobian_addition[i][d] = 0.;
}

template <>
template <typename InputVector>
void Equations<EquationsTypeMhd, 3>::Q(std_cxx11::array<typename InputVector::value_type, n_components> &result, const InputVector &W, const Tensor<1, 3> &normal) const
{
  std_cxx11::array<typename InputVector::value_type, n_components> forResult;
  typename InputVector::value_type b = asin(normal[2]);
  typename InputVector::value_type cb = cos(b);
  typename InputVector::value_type a;
  if (std::abs(normal[0]) < 1e-8)
  {
    if (std::abs(cb) < 1e-8)
      a = 0.;
    else
      a = asin(normal[1] / cb);
  }
  else
    a = acos(normal[0] / cb);

  typename InputVector::value_type sa = sin(a);
  typename InputVector::value_type sb = normal[2];
  typename InputVector::value_type ca = cos(a);

  forResult[density_component] = W[density_component];

  forResult[first_momentum_component]
    = (ca * cb * W[first_momentum_component])
    + (sa * cb * W[first_momentum_component + 1])
    + (sb * W[first_momentum_component + 2]);

  forResult[first_momentum_component + 1]
    = (-sa * W[first_momentum_component])
    + (ca * W[first_momentum_component + 1]);

  forResult[first_momentum_component + 2]
    = (-ca * sb * W[first_momentum_component])
    - (sa * sb * W[first_momentum_component + 1])
    + (cb * W[first_momentum_component + 2]);

  forResult[first_magnetic_flux_component]
    = (ca * cb * W[first_magnetic_flux_component])
    + (sa * cb * W[first_magnetic_flux_component + 1])
    + (sb * W[first_magnetic_flux_component + 2]);

  forResult[first_magnetic_flux_component + 1]
    = (-sa * W[first_magnetic_flux_component])
    + (ca * W[first_magnetic_flux_component + 1]);

  forResult[first_magnetic_flux_component + 2]
    = (-ca * sb * W[first_magnetic_flux_component])
    - (sa * sb * W[first_magnetic_flux_component + 1])
    + (cb * W[first_magnetic_flux_component + 2]);

  forResult[energy_component] = W[energy_component];

  for (unsigned int d = 0; d < n_components; d++)
    result[d] = forResult[d];
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::Q_inv(std_cxx11::array<typename InputVector::value_type, n_components> &result, std_cxx11::array<typename InputVector::value_type, n_components> &F, const Tensor<1, dim> &normal) const
{
  std_cxx11::array<typename InputVector::value_type, n_components> forResult;
  typename InputVector::value_type b = asin(normal[2]);
  typename InputVector::value_type cb = cos(b);
  typename InputVector::value_type a;
  if (std::abs(normal[0]) < 1e-8)
  {
    if (std::abs(cb) < 1e-8)
      a = 0.;
    else
      a = asin(normal[1] / cb);
  }
  else
    a = acos(normal[0] / cb);

  typename InputVector::value_type sa = sin(a);
  typename InputVector::value_type sb = normal[2];
  typename InputVector::value_type ca = cos(a);

  forResult[density_component] = F[density_component];

  forResult[first_momentum_component]
    = (ca * cb * F[first_momentum_component])
    - (sa * F[first_momentum_component + 1])
    - (ca * sb * F[first_momentum_component + 2]);

  forResult[first_momentum_component + 1]
    = (sa * cb * F[first_momentum_component])
    + (ca * F[first_momentum_component + 1])
    - (sa * sb * F[first_momentum_component + 2]);

  forResult[first_momentum_component + 2]
    = sb * F[first_momentum_component]
    + cb * F[first_momentum_component + 2];

  forResult[first_magnetic_flux_component]
    = (ca * cb * F[first_magnetic_flux_component])
    - (sa * F[first_magnetic_flux_component + 1])
    - (ca * sb * F[first_magnetic_flux_component + 2]);

  forResult[first_magnetic_flux_component + 1]
    = (sa * cb * F[first_magnetic_flux_component])
    + (ca * F[first_magnetic_flux_component + 1])
    - (sa * sb * F[first_magnetic_flux_component + 2]);

  forResult[first_magnetic_flux_component + 2]
    = sb * F[first_magnetic_flux_component]
    + cb * F[first_magnetic_flux_component + 2];

  forResult[energy_component] = F[energy_component];

  for (unsigned int d = 0; d < n_components; d++)
    result[d] = forResult[d];
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::smallest_eigenvalue(const InputVector &W) const
{
  return std::sqrt(
    (0.5 / W[0]) * (
    (this->parameters.gas_gamma * this->compute_pressure(W)) + (2. * this->compute_magnetic_energy(W))
      - std::sqrt(
      ((this->parameters.gas_gamma * this->compute_pressure(W)) + (2. * this->compute_magnetic_energy(W))) * ((this->parameters.gas_gamma * this->compute_pressure(W)) + (2. * this->compute_magnetic_energy(W)))
        - (4. * this->parameters.gas_gamma * this->compute_pressure(W) * W[4] * W[4])
      )
      )
  );
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::largest_eigenvalue(const InputVector &W) const
{
  return std::sqrt(
    (0.5 / W[0]) * (
    (this->parameters.gas_gamma * this->compute_pressure(W)) + (2. * this->compute_magnetic_energy(W))
      + std::sqrt(
      ((this->parameters.gas_gamma * this->compute_pressure(W)) + (2. * this->compute_magnetic_energy(W))) * ((this->parameters.gas_gamma * this->compute_pressure(W)) + (2. * this->compute_magnetic_energy(W)))
        - (4. * this->parameters.gas_gamma * this->compute_pressure(W) * W[4] * W[4])
      )
      )
  );
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::left_signal_speed(const InputVector &WL, const InputVector &WR) const
{
  return std::min(smallest_eigenvalue(WL), smallest_eigenvalue(WR));
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::right_signal_speed(const InputVector &WL, const InputVector &WR) const
{
  return std::max(largest_eigenvalue(WL), largest_eigenvalue(WR));
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
  std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux) const
{
  // If we use Lax Friedrich's flux.
  if (this->parameters.num_flux_type == Parameters<dim>::lax_friedrich)
  {
    std_cxx11::array<std_cxx11::array <typename InputVector::value_type, dim>, n_components > iflux, oflux;

    compute_flux_matrix(Wplus, iflux);
    compute_flux_matrix(Wminus, oflux);

    for (unsigned int di = 0; di < n_components; ++di)
    {
      normal_flux[di] = 0;
      for (unsigned int d = 0; d < dim; ++d)
        normal_flux[di] += 0.5*(iflux[di][d] + oflux[di][d]) * normal[d];

      normal_flux[di] += 0.5*this->parameters.lax_friedrich_stabilization_value*(Wplus[di] - Wminus[di]);
    }
  }

  // If we use HLLD
  if (this->parameters.num_flux_type == Parameters<dim>::hlld)
  {
    typename InputVector::value_type S_L, S_R;
    std_cxx11::array<typename InputVector::value_type, n_components> QedWL, QedWR;

    Q(QedWL, Wplus, normal);
    Q(QedWR, Wminus, normal);

    // Edge states.
    S_L = this->left_signal_speed(QedWL, QedWR);
    if (S_L >= 0.)
    {
      this->compute_flux_vector(0, QedWL, normal_flux);
      Q_inv<InputVector>(normal_flux, normal_flux, normal);
      return;
    }

    S_R = this->right_signal_speed(QedWL, QedWR);
    if (S_R <= 0.)
    {
      this->compute_flux_vector(0, QedWR, normal_flux);
      Q_inv<InputVector>(normal_flux, normal_flux, normal);
      return;
    }

    // Now we need to calculate S_M & all other stuff.
    typename InputVector::value_type u_L, u_R, rho_L, rho_R, p_L, p_R, p_T_R, p_T_L, S_M, B_x, psT;
    rho_L = QedWL[0];
    u_L = QedWL[1] / QedWL[0];
    p_L = this->compute_pressure(QedWL);
    p_T_L = p_L + this->compute_magnetic_energy(QedWL);
    rho_R = QedWR[0];
    u_R = QedWR[1] / QedWR[0];
    p_R = this->compute_pressure(QedWR);
    p_T_R = p_R + this->compute_magnetic_energy(QedWR);
    S_M = (((S_R - u_R) * (rho_R * u_R)) - ((S_L - u_L) * (rho_L * u_L)) - (p_T_R) + (p_T_L)) / (((S_R - u_R) * rho_R) - ((S_L - u_L) * rho_L));
    B_x = 0.5 * (QedWL[4] + QedWR[4]);
    psT = (((S_R - u_R) * (rho_R * p_T_L)) - ((S_L - u_L) * (rho_L * p_T_R)) + ((rho_L * rho_R) * (S_R - u_R) * (S_L - u_L) * (u_R - u_L))) / (((S_R - u_R) * rho_R) - ((S_L - u_L) * rho_L));

    std_cxx11::array<typename InputVector::value_type, n_components> Us_L;
    {
      Us_L[0] = rho_L * (S_L - u_L) / (S_L - S_M);
      Us_L[1] = Us_L[0] * S_M;
      Us_L[2] = Us_L[0] * (QedWL[2] / rho_L) - (B_x * QedWL[5] * ((S_M - u_L) / ((rho_L * (S_L - u_L) * (S_L - S_M)) - (B_x * B_x))));
      Us_L[3] = Us_L[0] * (QedWL[3] / rho_L) - (B_x * QedWL[6] * ((S_M - u_L) / ((rho_L * (S_L - u_L) * (S_L - S_M)) - (B_x * B_x))));
      // ???
      Us_L[4] = B_x;
      Us_L[5] = QedWL[5] * ((rho_L * (S_L - u_L) * (S_L - u_L)) - (B_x * B_x)) / ((rho_L * (S_L - u_L) * (S_L - S_M)) - (B_x * B_x));
      Us_L[6] = QedWL[6] * ((rho_L * (S_L - u_L) * (S_L - u_L)) - (B_x * B_x)) / ((rho_L * (S_L - u_L) * (S_L - S_M)) - (B_x * B_x));
      Us_L[7] = (((S_L - u_L) * QedWL[7]) - (p_T_L * u_L) + (psT * S_M) + (B_x * (((QedWL[1] * QedWL[4] + QedWL[2] * QedWL[5] + QedWL[3] * QedWL[6]) / rho_L) - ((Us_L[1] * Us_L[4] + Us_L[2] * Us_L[5] + Us_L[3] * Us_L[6]) / Us_L[0])))) / (S_L - S_M);
    }

    std_cxx11::array<typename InputVector::value_type, n_components> Us_R;
    {
      Us_R[0] = rho_R * (S_R - u_R) / (S_R - S_M);
      Us_R[1] = Us_R[0] * S_M;
      Us_R[2] = Us_R[0] * (QedWR[2] / rho_R) - (B_x * QedWR[5] * ((S_M - u_R) / ((rho_R * (S_R - u_R) * (S_R - S_M)) - (B_x * B_x))));
      Us_R[3] = Us_R[0] * (QedWR[3] / rho_R) - (B_x * QedWR[6] * ((S_M - u_R) / ((rho_R * (S_R - u_R) * (S_R - S_M)) - (B_x * B_x))));
      // ???
      Us_R[4] = B_x;
      Us_R[5] = QedWR[5] * ((rho_R * (S_R - u_R) * (S_R - u_R)) - (B_x * B_x)) / ((rho_R * (S_R - u_R) * (S_R - S_M)) - (B_x * B_x));
      Us_R[6] = QedWR[6] * ((rho_R * (S_R - u_R) * (S_R - u_R)) - (B_x * B_x)) / ((rho_R * (S_R - u_R) * (S_R - S_M)) - (B_x * B_x));
      Us_R[7] = (((S_R - u_R) * QedWR[7]) - (p_T_R * u_R) + (psT * S_M) + (B_x * (((QedWR[1] * QedWR[4] + QedWR[2] * QedWR[5] + QedWR[3] * QedWR[6]) / rho_R) - ((Us_R[1] * Us_R[4] + Us_R[2] * Us_R[5] + Us_R[3] * Us_R[6]) / Us_R[0])))) / (S_R - S_M);
    }

    typename InputVector::value_type Ss_L = S_M - (std::abs(B_x) / std::sqrt(Us_L[0]));
    typename InputVector::value_type Ss_R = S_M + (std::abs(B_x) / std::sqrt(Us_R[0]));

    // Fs_L, Fss_L
    if (S_M >= 0.)
    {
      std_cxx11::array<typename InputVector::value_type, n_components> F_L;
      this->compute_flux_vector(0, QedWL, F_L);
      // THIS is different to RIGHT state.
      if (Ss_L >= 0.)
      {
        std_cxx11::array<typename InputVector::value_type, n_components> Fs_L;
        for (int k = 0; k < n_components; k++)
          Fs_L[k] = F_L[k] + S_L * (Us_L[k] - QedWL[k]);
        Q_inv<InputVector>(normal_flux, Fs_L, normal);
        return;
      }
      else
      {
        std_cxx11::array<typename InputVector::value_type, n_components> Uss_L;
        Uss_L[0] = Us_L[0];
        Uss_L[1] = Us_L[1];
        Uss_L[2] = Uss_L[0] * (std::sqrt(Us_L[0]) * (Us_L[2] / Us_L[0]) + std::sqrt(Us_R[0]) * (Us_R[2] / Us_R[0]) + ((Us_R[5] - Us_L[5]) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_L[0]) + std::sqrt(Us_R[0]));
        Uss_L[3] = Uss_L[0] * (std::sqrt(Us_L[0]) * (Us_L[3] / Us_L[0]) + std::sqrt(Us_R[0]) * (Us_R[3] / Us_R[0]) + ((Us_R[6] - Us_L[6]) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_L[0]) + std::sqrt(Us_R[0]));
        // ???
        Uss_L[4] = B_x;
        Uss_L[5] = ((std::sqrt(Us_L[0]) * Us_L[5]) + (std::sqrt(Us_R[0]) * Us_R[5]) + std::sqrt(Us_L[0] * Us_R[0]) * (((Us_R[2] / Us_R[0]) - (Us_L[2] / Us_L[0])) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_L[0]) + std::sqrt(Us_R[0]));
        Uss_L[6] = ((std::sqrt(Us_L[0]) * Us_L[6]) + (std::sqrt(Us_R[0]) * Us_R[6]) + std::sqrt(Us_L[0] * Us_R[0]) * (((Us_R[3] / Us_R[0]) - (Us_L[3] / Us_L[0])) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_L[0]) + std::sqrt(Us_R[0]));
        // THIS is different to RIGHT state.
        Uss_L[7] = Us_L[7] - (std::sqrt(Us_L[0]) * (((Us_L[1] * Us_L[4] + Us_L[2] * Us_L[5] + Us_L[2] * Us_L[6]) / Us_L[0]) - ((Uss_L[1] * Uss_L[4] + Uss_L[2] * Uss_L[5] + Uss_L[2] * Uss_L[6]) / Uss_L[0])) * (B_x > 0. ? 1. : -1.));

        std_cxx11::array<typename InputVector::value_type, n_components> Fss_L;
        for (int k = 0; k < n_components; k++)
          Fss_L[k] = F_L[k] + Ss_L * Uss_L[k] - ((Ss_L - S_L) * Us_L[k]) - (S_L * QedWL[k]);
        Q_inv<InputVector>(normal_flux, Fss_L, normal);
        return;
      }
    }
    // Fs_R, Fss_R
    else
    {
      std_cxx11::array<typename InputVector::value_type, n_components> F_R;
      this->compute_flux_vector(0, QedWR, F_R);
      // THIS is different to LEFT state.
      if (Ss_R <= 0.)
      {
        std_cxx11::array<typename InputVector::value_type, n_components> Fs_R;
        for (int k = 0; k < n_components; k++)
          Fs_R[k] = F_R[k] + S_R * (Us_R[k] - QedWR[k]);
        Q_inv<InputVector>(normal_flux, Fs_R, normal);
        return;
      }
      else
      {
        std_cxx11::array<typename InputVector::value_type, n_components> Uss_R;
        Uss_R[0] = Us_R[0];
        Uss_R[1] = Us_R[1];
        Uss_R[2] = Uss_R[0] * (std::sqrt(Us_R[0]) * (Us_R[2] / Us_R[0]) + std::sqrt(Us_R[0]) * (Us_R[2] / Us_R[0]) + ((Us_R[5] - Us_R[5]) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_R[0]) + std::sqrt(Us_R[0]));
        Uss_R[3] = Uss_R[0] * (std::sqrt(Us_R[0]) * (Us_R[3] / Us_R[0]) + std::sqrt(Us_R[0]) * (Us_R[3] / Us_R[0]) + ((Us_R[6] - Us_R[6]) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_R[0]) + std::sqrt(Us_R[0]));
        // ???
        Uss_R[4] = B_x;
        Uss_R[5] = ((std::sqrt(Us_R[0]) * Us_R[5]) + (std::sqrt(Us_R[0]) * Us_R[5]) + std::sqrt(Us_R[0] * Us_R[0]) * (((Us_R[2] / Us_R[0]) - (Us_R[2] / Us_R[0])) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_R[0]) + std::sqrt(Us_R[0]));
        Uss_R[6] = ((std::sqrt(Us_R[0]) * Us_R[6]) + (std::sqrt(Us_R[0]) * Us_R[6]) + std::sqrt(Us_R[0] * Us_R[0]) * (((Us_R[3] / Us_R[0]) - (Us_R[3] / Us_R[0])) * (B_x > 0. ? 1. : -1.))) / (std::sqrt(Us_R[0]) + std::sqrt(Us_R[0]));
        // THIS is different to RIGHT state.
        Uss_R[7] = Us_R[7] + (std::sqrt(Us_R[0]) * (((Us_R[1] * Us_R[4] + Us_R[2] * Us_R[5] + Us_R[2] * Us_R[6]) / Us_R[0]) - ((Uss_R[1] * Uss_R[4] + Uss_R[2] * Uss_R[5] + Uss_R[2] * Uss_R[6]) / Uss_R[0])) * (B_x > 0. ? 1. : -1.));

        std_cxx11::array<typename InputVector::value_type, n_components> Fss_R;
        for (int k = 0; k < n_components; k++)
          Fss_R[k] = F_R[k] + Ss_R * Uss_R[k] - ((Ss_R - S_R) * Us_R[k]) - (S_R * QedWR[k]);
        Q_inv<InputVector>(normal_flux, Fss_R, normal);
        return;
      }
    }
  }
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing) const
{
  for (unsigned int c = 0; c < n_components; ++c)
    forcing[c] = 0.;
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
    const double density = inputs.solution_values[q](density_component);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = inputs.solution_values[q](first_momentum_component + d) / density;

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

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = uh[q](density_component);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = uh[q](first_momentum_component + d) / density;

    computed_quantities[q](dim) = equations.compute_pressure(uh[q]);
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
  return update_values;
}

template class Equations<EquationsTypeMhd, 3>;

#include "equationsMhd_inst.cpp"
