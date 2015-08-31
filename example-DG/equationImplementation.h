class EquationImplementation
{
public:

  static d matrixVolValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, Point<DIM> quadPoint);

  static d matrixBoundaryEdgeValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal);

  static d matrixInternalEdgeValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N,
    Point<DIM> quadPoint, Point<DIM> normal);


  static d rhsVolValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, Point<DIM> quadPoint);

  static d rhsBoundaryEdgeValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal);

  static d rhsInternalEdgeValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad,
    d v_valN, vec Un_valN, dimVec v_gradN,
    vecDimVec Un_gradN, bool v_N, Point<DIM> quadPoint, Point<DIM> normal);
  
  static void Jacobians(FullMatrix<double> *J,
                        std::vector<Vector<double> > lv,
                        const unsigned int qp)
};