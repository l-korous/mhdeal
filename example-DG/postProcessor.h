#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H

#include <deal.II/numerics/data_out.h>

#include "definitions.h"

class Postprocessor : public dealii::DataPostprocessor < DIM >
{
public:
  Postprocessor();

  void compute_derived_quantities_scalar(const std::vector<d>         &/*uh*/,
      const std::vector<dealii::Tensor<1, DIM> > &/*duh*/,
      const std::vector<dealii::Tensor<2, DIM> > &/*dduh*/,
      const std::vector<dealii::Point<DIM> >    &/*normals*/,
      const std::vector<dealii::Point<DIM> >    &/*evaluation_points*/,
      std::vector<dealii::Vector<double> >      &computed_quantities) const;

  void compute_derived_quantities_vector(const std::vector<dealii::Vector<d> > &uh,
    const std::vector<vecDimVec> &duh,
    const std::vector<std::vector<dealii::Tensor<2, DIM> > > &dduh,
    const std::vector<dealii::Point<DIM> > &normals,
    const std::vector<dealii::Point<DIM> > &evaluation_points,
    std::vector<dealii::Vector<d> > &computed_quantities) const;

  virtual std::vector<std::string> get_names() const;
  virtual std::vector < dealii::DataComponentInterpretation::DataComponentInterpretation > get_data_component_interpretation() const;
  virtual dealii::UpdateFlags get_needed_update_flags() const;
};

#endif
