#ifndef COLLISION_POINT_CLOUD_H
#define COLLISION_POINT_CLOUD_H

#include <meshing/PointCloud.h>
#include <math3d/geometry3d.h>
#include <limits.h>
#include <meshing/Voxelize.h>
#include <meshing/VolumeGrid.h>
#include <geometry/primitives.h>
#include <cmath>
#include <geometry/CollisionMesh.h>
#include <Timer.h>

namespace Geometry {

  using namespace Math3D;

class CollisionPointCloud : public Meshing::PointCloud3D
{
 public:
  CollisionPointCloud();
  CollisionPointCloud(const Meshing::PointCloud3D& pc);
  ///Needs to be called if this point cloud was loaded or set up any other
  ///way than the constructor
  void InitCollisions();

  //TODO: collision accelerators
  AABB3D bblocal;
  RigidTransform currentTransform;
};

void GetBB(const CollisionPointCloud& pc,Box3D& b);
bool WithinDistance(const CollisionPointCloud& pc,const GeometricPrimitive3D& g,Real tol);
Real Distance(const CollisionPointCloud& pc,const GeometricPrimitive3D& g);
void NearbyPoints(const CollisionPointCloud& pc,const GeometricPrimitive3D& g,Real tol,std::vector<int>& points,size_t maxContacts=INT_MAX);

struct SpatialHashPC
{
	std::vector<std::vector<Vector3> > table;
	AABB3D bb;
	Box3D currentBox;
	Real cellSize;
	RigidTransform currentTransform;

	SpatialHashPC(CollisionPointCloud, Real);
	~SpatialHashPC();

	//void getMinMaxPoints(Vector3&, Vector3&, Geometry::CollisionPointCloud3D);
	Vector3 getDimension();
	int getGridX();
	int getGridY();
	int getGridZ();
	int getTableSize();
	void Transform(RigidTransform&);

	int Hash(Vector3);
	int Hash(Real x, Real y, Real z);

	void addPoint(Vector3);
	bool cellEmpty(Vector3);
	bool cellEmpty(Real, Real, Real);
	bool cellEmpty(int);
};

struct EuclideanDistanceTransform
{	
	AABB3D aabb;
	Box3D bb;
	RigidTransform currentTransform;
	Real cellSize;   			 //Needs to be set
	Geometry::CollisionMesh m;   //Needs to be initialized
	Array3D<Real>* distance; 	 //Initialize the size
	Array3D<Vector3>* gradient;
	std::vector<IntTriple> surfaceCells;

	EuclideanDistanceTransform(CollisionMesh, Real);
	~EuclideanDistanceTransform();

	void Transform(RigidTransform&);
	void FitAABB();
};

//bool Collide(CollisionPointCloud&, CollisionMesh&, std::vector<Meshing::Vector3>&);
bool Collide(SpatialHashPC&, EuclideanDistanceTransform&, std::vector<Meshing::Vector3>&);

bool Collide(CollisionPointCloud&, EuclideanDistanceTransform&);

bool NewCollide(SpatialHashPC&, EuclideanDistanceTransform&);


} //namespace Geometry

#endif
