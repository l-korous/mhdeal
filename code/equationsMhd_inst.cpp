

template typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type Equations<2>::compute_kinetic_energy(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W);

template typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type Equations<2>::compute_pressure(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W);

template void Equations<2>::compute_flux_matrix(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W, std_cxx11::array <std_cxx11::array <typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type, 2>, Equations<2>::n_components > &flux);

template void Equations<2>::numerical_normal_flux(const Tensor<1, 2> &normal, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wplus, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wminus, double alpha,
  std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type, n_components> &normal_flux);

template void Equations<2>::compute_forcing_vector(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W, std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type, n_components> &forcing);

template void Equations<2>::compute_Wminus(const BoundaryKind(&boundary_kind)[n_components], const Tensor<1, 2> &normal_vector, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wplus, const Vector<double> &boundary_values, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wminus);




template typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type Equations<3>::compute_kinetic_energy(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W);

template typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type Equations<3>::compute_pressure(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W);

template void Equations<3>::compute_flux_matrix(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W, std_cxx11::array <std_cxx11::array <typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type, 3>, Equations<3>::n_components > &flux);

template void Equations<3>::numerical_normal_flux(const Tensor<1, 3> &normal, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wplus, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wminus, double alpha,
  std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type, n_components> &normal_flux);

template void Equations<3>::compute_forcing_vector(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W, std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type, n_components> &forcing);

template void Equations<3>::compute_Wminus(const BoundaryKind(&boundary_kind)[n_components], const Tensor<1, 3> &normal_vector, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wplus, const Vector<double> &boundary_values, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &Wminus);



template typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type Equations<2>::compute_kinetic_energy(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W);

template typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type Equations<2>::compute_pressure(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W);

template void Equations<2>::compute_flux_matrix(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W, std_cxx11::array <std_cxx11::array <typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type, 2>, Equations<2>::n_components > &flux);

template void Equations<2>::numerical_normal_flux(const Tensor<1, 2> &normal, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &Wplus, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &Wminus, double alpha,
  std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type, n_components> &normal_flux);

template void Equations<2>::compute_forcing_vector(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W, std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type, n_components> &forcing);

template void Equations<2>::compute_Wminus(const BoundaryKind(&boundary_kind)[n_components], const Tensor<1, 2> &normal_vector, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> &Wplus, const Vector<double> &boundary_values, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> &Wminus);




template typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type Equations<3>::compute_kinetic_energy(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W);

template typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type Equations<3>::compute_pressure(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W);

template void Equations<3>::compute_flux_matrix(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W, std_cxx11::array <std_cxx11::array <typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type, 3>, Equations<3>::n_components > &flux);

template void Equations<3>::numerical_normal_flux(const Tensor<1, 3> &normal, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &Wplus, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &Wminus, double alpha,
  std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type, n_components> &normal_flux);

template void Equations<3>::compute_forcing_vector(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1> &W, std_cxx11::array<typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double> , false, 1>::value_type, n_components> &forcing);

template void Equations<3>::compute_Wminus(const BoundaryKind(&boundary_kind)[n_components], const Tensor<1, 3> &normal_vector, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> &Wplus, const Vector<double> &boundary_values, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> &Wminus);

