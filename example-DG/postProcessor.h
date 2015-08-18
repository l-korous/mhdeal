#ifndef POST_PROCESSOR_H
#define POST_PROCESSOR_H

#include "definitions.h"

class Postprocessor : public DataPostprocessor < DIM >
{
public:
  Postprocessor();
  virtual void compute_derived_quantities_vector(const std::vector<Vector<d> > &uh,
    const std::vector<vecDimVec> &duh,
    const std::vector<std::vector<Tensor<2, DIM> > > &dduh,
    const std::vector<Point<DIM> > &normals,
    const std::vector<Point<DIM> > &evaluation_points,
    std::vector<Vector<d> > &computed_quantities) const;

  virtual std::vector<std::string> get_names() const;
  virtual std::vector < DataComponentInterpretation::DataComponentInterpretation > get_data_component_interpretation() const;
  virtual UpdateFlags get_needed_update_flags() const;
};

#endif // POST_PROCESSOR_H
