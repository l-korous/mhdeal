#ifndef INITIAL_SLN_H
#define INITIAL_SLN_H

#include "definitions.h"

class InitialSlnRho
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnMomentumX
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnMomentumY
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnMomentumZ
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnEnergy
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnB1
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnB2
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnB3
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnJ1
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnJ2
{
public:
  static d value(const dealii::Point<DIM> &p);
};

class InitialSlnJ3
{
public:
  static d value(const dealii::Point<DIM> &p);
};

#endif
