#include "util.h"
#include "fe_is_constant_interface.h"

template <int dim, int spacedim = dim>
class FE_DG_DivFree : public FiniteElement<dim, spacedim>, public FiniteElementIsConstantInterface<dim>
{
public:
  FE_DG_DivFree();

  virtual std::string get_name() const;

  virtual UpdateFlags requires_update_flags(const UpdateFlags update_flags) const;

  bool is_constant(const unsigned int i) const;

  /**
  * Return the value of the
  * @p ith shape function at the
  * point @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function. This function is not
  * implemented.
  */
  virtual double shape_value(const unsigned int i,
    const Point<dim> &p) const;

  /**
  * Return the value of the
  * @p ith shape function at the
  * point @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  */
  double shape_value
  (const typename Triangulation<dim, spacedim>::cell_iterator & cell,
    const unsigned int i,
    const Point<dim> &p) const;


  /**
  * Compute all shape function values at the specified points @p p.
  * @p values[i][q] is the @p i th shape function value at point @p p[q].
  */
  void
    shape_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell,
      const std::vector< Point<dim> > &p,
      std::vector< std::vector<double> >& values) const;

  /**
  * Return the value of the
  * @p componentth vector
  * component of the @p ith shape
  * function at the point
  * @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  *
  * Since this element is scalar,
  * the returned value is the same
  * as if the function without the
  * @p _component suffix were
  * called, provided that the
  * specified component is zero.
  * This function is not implemented.
  */
  virtual double shape_value_component(const unsigned int i,
    const Point<dim> &p,
    const unsigned int component) const;

  /**
  * Return the value of the
  * @p componentth vector
  * component of the @p ith shape
  * function at the point
  * @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  *
  * Since this element is scalar,
  * the returned value is the same
  * as if the function without the
  * @p _component suffix were
  * called, provided that the
  * specified component is zero.
  */
  double shape_value_component
  (const typename Triangulation<dim, spacedim>::cell_iterator & cell,
    const unsigned int i,
    const Point<dim> &p,
    const unsigned int component) const;

  double shape_value_component
  (const typename Triangulation<dim, spacedim>::cell_iterator & cell,
    const unsigned int i,
    const Point<dim> &p,
    const unsigned int component, double cell_diameter) const;

  /**
  * Return the gradient of the
  * @p ith shape function at the
  * point @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function. This function is not
  * implemented.
  */
  virtual Tensor<1, dim> shape_grad(const unsigned int  i,
    const Point<dim>   &p) const;

  /**
  * Return the gradient of the
  * @p ith shape function at the
  * point @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  */
  Tensor<1, dim> shape_grad
  (const typename Triangulation<dim, spacedim>::cell_iterator & cell,
    const unsigned int  i,
    const Point<dim>   &p) const;

  /**
  * Return the gradient of the
  * @p componentth vector
  * component of the @p ith shape
  * function at the point
  * @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  *
  * Since this element is scalar,
  * the returned value is the same
  * as if the function without the
  * @p _component suffix were
  * called, provided that the
  * specified component is zero.
  * This function is not implemented.
  */
  virtual Tensor<1, dim> shape_grad_component(const unsigned int i,
    const Point<dim> &p,
    const unsigned int component) const;

  /**
  * Return the gradient of the
  * @p componentth vector
  * component of the @p ith shape
  * function at the point
  * @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  *
  * Since this element is scalar,
  * the returned value is the same
  * as if the function without the
  * @p _component suffix were
  * called, provided that the
  * specified component is zero.
  */
  Tensor<1, dim> shape_grad_component
  (const typename Triangulation<dim, spacedim>::cell_iterator & cell,
    const unsigned int i,
    const Point<dim> &p,
    const unsigned int component) const;

  /**
  * X Return the tensor of second
  * X derivatives of the @p ith
  * X shape function at point @p p
  * X on the unit cell.  See the
  * X FiniteElement base
  * X class for more information
  * X about the semantics of this
  * X function. This function is not
  * X implemented.
  */
  virtual Tensor<2, dim> shape_grad_grad(const unsigned int  i,
    const Point<dim> &p) const;

  /**
  * Return the tensor of second
  * derivatives of the @p ith
  * shape function at point @p p
  * on the unit cell.  See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  */
  Tensor<2, dim> shape_grad_grad
  (const typename Triangulation<dim, spacedim>::cell_iterator & cell,
    const unsigned int  i,
    const Point<dim> &p) const;

  /**
  * Return the second derivative
  * of the @p componentth vector
  * component of the @p ith shape
  * function at the point
  * @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  *
  * Since this element is scalar,
  * the returned value is the same
  * as if the function without the
  * @p _component suffix were
  * called, provided that the
  * specified component is zero.
  * This function is not implemented.
  */
  virtual Tensor<2, dim> shape_grad_grad_component(const unsigned int i,
    const Point<dim> &p,
    const unsigned int component) const;

  /**
  * Return the second derivative
  * of the @p componentth vector
  * component of the @p ith shape
  * function at the point
  * @p p. See the
  * FiniteElement base
  * class for more information
  * about the semantics of this
  * function.
  *
  * Since this element is scalar,
  * the returned value is the same
  * as if the function without the
  * @p _component suffix were
  * called, provided that the
  * specified component is zero.
  */
  Tensor<2, dim> shape_grad_grad_component
  (const typename Triangulation<dim, spacedim>::cell_iterator & cell,
    const unsigned int i,
    const Point<dim> &p,
    const unsigned int component) const;

  /**
  * Return the polynomial degree
  * of this finite element,
  * i.e. the value passed to the
  * constructor.
  */
  unsigned int get_degree() const;

  /**
  * Return the matrix
  * interpolating from a face of
  * of one element to the face of
  * the neighboring element.
  * The size of the matrix is
  * then <tt>source.dofs_per_face</tt> times
  * <tt>this->dofs_per_face</tt>.
  *
  * Derived elements will have to
  * implement this function. They
  * may only provide interpolation
  * matrices for certain source
  * finite elements, for example
  * those from the same family. If
  * they don't implement
  * interpolation from a given
  * element, then they must throw
  * an exception of type
  * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
  */
  virtual void
    get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
      FullMatrix<double>       &matrix) const;

  /**
  * Return the matrix
  * interpolating from a face of
  * of one element to the face of
  * the neighboring element.
  * The size of the matrix is
  * then <tt>source.dofs_per_face</tt> times
  * <tt>this->dofs_per_face</tt>.
  *
  * Derived elements will have to
  * implement this function. They
  * may only provide interpolation
  * matrices for certain source
  * finite elements, for example
  * those from the same family. If
  * they don't implement
  * interpolation from a given
  * element, then they must throw
  * an exception of type
  * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
  */
  virtual void
    get_subface_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
      const unsigned int        subface,
      FullMatrix<double>       &matrix) const;

  /**
  * @name Functions to support hp
  * @{
  */

  /**
  * If, on a vertex, several finite elements are active, the hp code
  * first assigns the degrees of freedom of each of these FEs
  * different global indices. It then calls this function to find out
  * which of them should get identical values, and consequently can
  * receive the same global DoF index. This function therefore
  * returns a list of identities between DoFs of the present finite
  * element object with the DoFs of @p fe_other, which is a reference
  * to a finite element object representing one of the other finite
  * elements active on this particular vertex. The function computes
  * which of the degrees of freedom of the two finite element objects
  * are equivalent, both numbered between zero and the corresponding
  * value of dofs_per_vertex of the two finite elements. The first
  * index of each pair denotes one of the vertex dofs of the present
  * element, whereas the second is the corresponding index of the
  * other finite element.
  *
  * This being a discontinuous element, the set of such constraints
  * is of course empty.
  */
  virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_vertex_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const;

  /**
  * Same as hp_vertex_dof_indices(), except that the function treats
  * degrees of freedom on lines.
  *
  * This being a discontinuous element, the set of such constraints
  * is of course empty.
  */
  virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_line_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const;

  /**
  * Same as hp_vertex_dof_indices(), except that the function treats
  * degrees of freedom on quads.
  *
  * This being a discontinuous element, the set of such constraints
  * is of course empty.
  */
  virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const;

  /**
  * Return whether this element
  * implements its hanging node
  * constraints in the new way,
  * which has to be used to make
  * elements "hp compatible".
  *
  * For the FE_DG_DivFree class the
  * result is always true (independent of
  * the degree of the element), as it has
  * no hanging nodes (being a
  * discontinuous element).
  */
  virtual bool hp_constraints_are_implemented() const;

  /**
  * Return whether this element dominates
  * the one given as argument when they
  * meet at a common face,
  * whether it is the other way around,
  * whether neither dominates, or if
  * either could dominate.
  *
  * For a definition of domination, see
  * FiniteElementBase::Domination and in
  * particular the @ref hp_paper "hp paper".
  */
  virtual
    FiniteElementDomination::Domination
    compare_for_face_domination(const FiniteElement<dim, spacedim> &fe_other) const;

  /**
  * @}
  */

  /**
  * Check for non-zero values on a face.
  *
  * This function returns
  * @p true, if the shape
  * function @p shape_index has
  * non-zero values on the face
  * @p face_index.
  *
  * Implementation of the
  * interface in
  * FiniteElement
  */
  virtual bool has_support_on_face(const unsigned int shape_index,
    const unsigned int face_index) const;

  /**
  * Determine an estimate for the
  * memory consumption (in bytes)
  * of this object.
  *
  * This function is made virtual,
  * since finite element objects
  * are usually accessed through
  * pointers to their base class,
  * rather than the class itself.
  */
  virtual std::size_t memory_consumption() const;

protected:

  /**
  * @p clone function instead of
  * a copy constructor.
  *
  * This function is needed by the
  * constructors of @p FESystem.
  */
  virtual FiniteElement<dim, spacedim> *clone() const;

  /**
  * Prepare internal data
  * structures and fill in values
  * independent of the cell.
  */
  virtual
    typename FiniteElement<dim, spacedim>::InternalDataBase *
    get_data(const UpdateFlags                                                    update_flags,
      const Mapping<dim, spacedim>                                         &mapping,
      const Quadrature<dim>                                               &quadrature,
      dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
    void
    fill_fe_values(const typename Triangulation<dim, spacedim>::cell_iterator           &cell,
      const CellSimilarity::Similarity                                     cell_similarity,
      const Quadrature<dim>                                               &quadrature,
      const Mapping<dim, spacedim>                                         &mapping,
      const typename Mapping<dim, spacedim>::InternalDataBase              &mapping_internal,
      const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
      const typename FiniteElement<dim, spacedim>::InternalDataBase        &fe_internal,
      dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
    void
    fill_fe_face_values(const typename Triangulation<dim, spacedim>::cell_iterator           &cell,
      const unsigned int                                                   face_no,
      const Quadrature<dim - 1>                                             &quadrature,
      const Mapping<dim, spacedim>                                         &mapping,
      const typename Mapping<dim, spacedim>::InternalDataBase              &mapping_internal,
      const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
      const typename FiniteElement<dim, spacedim>::InternalDataBase        &fe_internal,
      dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
    void
    fill_fe_subface_values(const typename Triangulation<dim, spacedim>::cell_iterator           &cell,
      const unsigned int                                                   face_no,
      const unsigned int                                                   sub_no,
      const Quadrature<dim - 1>                                             &quadrature,
      const Mapping<dim, spacedim>                                         &mapping,
      const typename Mapping<dim, spacedim>::InternalDataBase              &mapping_internal,
      const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
      const typename FiniteElement<dim, spacedim>::InternalDataBase        &fe_internal,
      dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const;

private:

  /**
  * Only for internal use. Its
  * full name is
  * @p get_dofs_per_object_vector
  * function and it creates the
  * @p dofs_per_object vector that is
  * needed within the constructor to
  * be passed to the constructor of
  * @p FiniteElementData.
  */
  static
    std::vector<unsigned int>
    get_dpo_vector(const unsigned int degree);

  /**
  * Allow access from other dimensions.
  */
  template <int, int> friend class FE_DG_DivFree;
};