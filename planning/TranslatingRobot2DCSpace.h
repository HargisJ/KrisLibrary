#ifndef TRANSLATING_ROBOT_2D_CSPACE_H
#define TRANSLATING_ROBOT_2D_CSPACE_H

#include "Geometric2DCSpace.h"

/** @brief a translating 2D robot in a 2D workspace.
 */
class TranslatingRobot2DCSpace : public ExplicitCSpace
{
public:
  TranslatingRobot2DCSpace();
  void DrawWorkspaceGL() const;
  void DrawRobotGL(const Config& x) const;
  void DrawGL(const Config& q) const;

  virtual void Sample(Config& x);
  virtual void SampleNeighborhood(const Config& c,Real r,Config& x);
  virtual int NumObstacles();
  virtual std::string ObstacleName(int obstacle);
  virtual bool IsFeasible(const Config& x);
  virtual bool IsFeasible(const Config& x,int obstacle);
  virtual EdgePlanner* LocalPlanner(const Config& a,const Config& b);
  virtual EdgePlanner* LocalPlanner(const Config& a,const Config& b,int obstacle);
  virtual Real Distance(const Config& x, const Config& y);

  Real visibilityEpsilon;
  AABB2D domain;
  Geometric2DCollection obstacles;
  Geometric2DCollection robot;
};

#endif
