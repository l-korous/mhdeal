class NumFlux
{
public:
  virtual void calculate(vec U_L, vec U_R, Point<DIM> normal, vec result) = 0;
};

class NumFluxLaxFriedrichs
{
public:
  void calculate(vec U_L, vec U_R, Point<DIM> normal, vec result);
};