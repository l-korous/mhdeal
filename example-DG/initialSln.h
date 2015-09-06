#ifndef INITIAL_SLN_H
#define INITIAL_SLN_H

class InitialSlnRho
{
public:
  static d value(const Point<DIM> &p);
};

class InitialSlnEnergy
{
public:
  static d value(const Point<DIM> &p);
};

class InitialSlnMomentumX
{
public:
  static d value(const Point<DIM> &p);
};

class InitialSlnMomentumY
{
public:
  static d value(const Point<DIM> &p);
};

class InitialSlnMomentumZ
{
public:
  static d value(const Point<DIM> &p);
};

#endif