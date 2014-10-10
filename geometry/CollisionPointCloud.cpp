#include "CollisionPointCloud.h"

namespace Geometry {

CollisionPointCloud::CollisionPointCloud()
{
  currentTransform.setIdentity();
}

CollisionPointCloud::CollisionPointCloud(const Meshing::PointCloud3D& _pc)
  :Meshing::PointCloud3D(_pc)
{
  currentTransform.setIdentity();
  InitCollisions();

void CollisionPointCloud::InitCollisions()
{
  bblocal.minimize();
  for(size_t i=0;i<points.size();i++)
    bblocal.expand(points[i]);
}

void GetBB(const CollisionPointCloud& pc,Box3D& b)
{
  b.setTransformed(pc.bblocal,pc.currentTransform);
}

bool WithinDistance(const CollisionPointCloud& pc,const GeometricPrimitive3D& g,Real tol)
{
  Box3D bb;
  GetBB(pc,bb);
  //quick reject test
  if(g.Distance(bb) > tol) return false;

  //test all points, linearly
  for(size_t i=0;i<pc.points.size();i++)
    if(g.Distance(pc.points[i]) <= tol) return true;
  return false;
}

Real Distance(const CollisionPointCloud& pc,const GeometricPrimitive3D& g)
{
  Real dmax = Inf;
  //test all points, linearly
  for(size_t i=0;i<pc.points.size();i++)
    dmax = Min(dmax,g.Distance(pc.points[i]));
  return dmax;
}

void NearbyPoints(const CollisionPointCloud& pc,const GeometricPrimitive3D& g,Real tol,std::vector<int>& points,size_t maxContacts)
{
  Box3D bb;
  GetBB(pc,bb);
  //quick reject test
  if(g.Distance(bb) > tol) return;

  //test all points, linearly
  for(size_t i=0;i<pc.points.size();i++)
    if(g.Distance(pc.points[i]) <= tol) {
      points.push_back(int(i));
      if(points.size()>=maxContacts) return;
    }
}

SpatialHashPC::SpatialHashPC(CollisionPointCloud pc, Real cSize)
{
  cellSize = cSize;
  currentTransform.setIdentity();
  // Vector3 bmin, bmax;
  // getMinMaxPoints(bmin, bmax, pc);
  bb = pc.bblocal;
  table.resize(getTableSize());
  //std::cout << pc.bblocal << "\n" << getTableSize() << "\n";
  currentBox.set(bb);

  std::vector<Vector3>::iterator it = pc.points.begin();
  while(it != pc.points.end())
  {
    addPoint(*it);
    it++;
  }
}

SpatialHashPC::~SpatialHashPC()
{

}

// void SpatialHashPC::getMinMaxPoints(Vector3& min, Vector3& max, Meshing::PointCloud3D pc)
// {
//   Real minX = 0; Real minY = 0; Real minZ = 0;
//   Real maxX = 0; Real maxY = 0; Real maxZ = 0;
//   std::vector<Vector3>::iterator it = pc.points.begin();
//   while(it != pc.points.end())
//   {
//     if (it->x < minX) minX = it->x;
//     else if (it->x > maxX) maxX = it->x;
//     if (it->y < minY) minY = it->y;
//     else if (it->y > maxY) maxY = it->y;
//     if (it->z < minZ) minZ = it->z;
//     else if (it->z > maxZ) maxZ = it->z;

//     it++;
//   }
//   min.set(minX, minY, minZ);
//   max.set(maxX, maxY, maxZ);
// }

Vector3 SpatialHashPC::getDimension()
{
  Vector3 dim (getGridX(),getGridY(),getGridZ());
  return dim;
}

int SpatialHashPC::getGridX()
{
  return(ceil((bb.bmax.x-bb.bmin.x)/cellSize));
}

int SpatialHashPC::getGridY()
{
  return(ceil((bb.bmax.y-bb.bmin.y)/cellSize));
}

int SpatialHashPC::getGridZ()
{
  return(ceil((bb.bmax.z-bb.bmin.z)/cellSize));
}

int SpatialHashPC::getTableSize()
{
  return(getGridX()*getGridY()*getGridZ());
}

void SpatialHashPC::Transform(RigidTransform& b)
{
  currentTransform *= b;
  currentBox.setTransformed(currentBox, b);
}

int SpatialHashPC::Hash(Vector3 point)
{
  if(!currentTransform.isIdentity())
    currentTransform.mulInverse(point, point);
  return(floor((point.x - bb.bmin.x)/cellSize) + getGridX()*floor((point.y - bb.bmin.y)/cellSize) + getGridX()*getGridY()*floor((point.z - bb.bmin.z)/cellSize));
}

int SpatialHashPC::Hash(Real x, Real y, Real z)
{
  Vector3 point(x,y,z);
  return(Hash(point));
}

void SpatialHashPC::addPoint(Vector3 pt)
{
  int index = Hash(pt);
  table[index].push_back(pt);
}

bool SpatialHashPC::cellEmpty(Vector3 pt)
{
  int index = Hash(pt);
  return table[index].empty();
}

bool SpatialHashPC::cellEmpty(Real x, Real y, Real z)
{
  Vector3 pt (x,y,z);
  return cellEmpty(pt);
}

bool SpatialHashPC::cellEmpty(int index)
{
  return table[index].empty();
}

EuclideanDistanceTransform::EuclideanDistanceTransform(CollisionMesh tm, Real cSize)
{
  cellSize = cSize;
  m = tm; 
  currentTransform.setIdentity();
  m.GetAABB(aabb.bmin, aabb.bmax);
  // std::cout << "(" << aabb.bmin.x << ", " << aabb.bmin.y << ", " << aabb.bmin.z << ") ("
  //   << aabb.bmax.x << ", " << aabb.bmax.y << ", " << aabb.bmax.z << ")\n";
  // std::cout.flush();
  FitAABB();
  bb.set(aabb);
  // std::cout << "(" << aabb.bmin.x << ", " << aabb.bmin.y << ", " << aabb.bmin.z << ") ("
  //   << aabb.bmax.x << ", " << aabb.bmax.y << ", " << aabb.bmax.z << ")\n";
  Vector3 gridSize((aabb.bmax.x-aabb.bmin.x)/cellSize,(aabb.bmax.y-aabb.bmin.y)/cellSize, (aabb.bmax.z-aabb.bmin.z)/cellSize);
  //std::cout << (bb.bmax.x-bb.bmin.x)/cellSize << " " << (bb.bmax.y-bb.bmin.y)/cellSize << " " << (bb.bmax.z-bb.bmin.z)/cellSize << "\n";
  distance = new Array3D<Real>(gridSize.x,gridSize.y,gridSize.z);
  gradient = new Array3D<Vector3>(gridSize.x,gridSize.y,gridSize.z);
  FastMarchingMethod(m, *distance, *gradient, aabb, surfaceCells);
}

EuclideanDistanceTransform::~EuclideanDistanceTransform()
{

}

void EuclideanDistanceTransform::Transform(RigidTransform& T)
{
  currentTransform *= T;
  bb.setTransformed(aabb, currentTransform);
}

void EuclideanDistanceTransform::FitAABB()
{
  Vector3 size=aabb.bmax-aabb.bmin;
  int iX, iY, iZ;
  iX = (size.x / cellSize);
  iY = (size.y / cellSize);
  iZ = (size.z / cellSize);
  if(iX*cellSize != size.x)
  {
    size.x = (1+iX)*cellSize;
  }
  if(iY*cellSize != size.y)
  {
    size.y = (1+iY)*cellSize;
  }
  if(iZ*cellSize != size.z)
  {
    size.z = (1+iZ)*cellSize;
  }

  aabb.bmax = aabb.bmin+size;
}


// bool Collide(CollisionPointCloud& pc, CollisionMesh& m, std::vector<Meshing::Vector3>& collisionPoints)
// {
//   Timer timer;
//   std::cout << "In collide" << "\n";
//   std::cout.flush();

//   AABB3D pcbb, mbb;
//   pcbb = pc.bblocal;
//   m.GetAABB(mbb.bmin, mbb.bmax);


//   if(!(pcbb.bmin.x <= mbb.bmax.x && pcbb.bmax.x >= mbb.bmin.x &&
// 	pcbb.bmin.y <= mbb.bmax.y && pcbb.bmax.y >= mbb.bmin.y &&
// 	pcbb.bmin.z <= mbb.bmax.z && pcbb.bmax.z <= mbb.bmin.z))
// 	{	
// 		return false;
// 	}

//   //SET CELL SIZES
//   Real pcCellSize = .4;
//   Real edtCellSize = .25;

//   SpatialHashPC sh(pc,pcCellSize);

//   std::cout << "SH built: " << "(" << sh.bb.bmin.x << " ," << sh.bb.bmin.y << " ," << sh.bb.bmin.z;
//   std::cout << ") (" << sh.bb.bmax.x << " ," << sh.bb.bmax.y << " ," << sh.bb.bmax.z << ")\n";
//   std::cout << sh.getGridX() << ", " << sh.getGridY() << ", " << sh.getGridZ() << "\n";
//   std::cout << "SH build time: " << timer.ElapsedTime() << "\n";
//   std::cout.flush();

//   timer.Reset();
  
//   EuclideanDistanceTransform edt(m, edtCellSize);
//   std::cout << edt.distance->m << " " << edt.distance->n << " " << edt.distance->p << "\n";

//   std::cout << "EDT built: (" << edt.aabb.bmin.x << ", " << edt.aabb.bmin.y << ", " << edt.aabb.bmin.z << ") ("
//     << edt.aabb.bmax.x << ", " << edt.aabb.bmax.y << ", " << edt.aabb.bmax.z << ")\n";
//   std::cout << "EDT build time: " << timer.ElapsedTime() << "\n";
//   std::cout.flush();

//   return Collide(sh, edt, collisionPoints);
// }


/*Function to test for collision between a pointcloud and mesh by first iterating
* through the EDT voxels and hashing the vertices to find the hash buckets that may
* be in collision. The points in the bucket are then tested against the EDT.
*
* This implementation is flawed. The list used to store the bucket indices has to be sorted to remove duplicates.
*Also, a set of surface and interior cells of the EDT should be stored in the object instead of having to be generated
* every time in the collision function.
*/
bool Collide(SpatialHashPC& sh, EuclideanDistanceTransform& edt, std::vector<Meshing::Vector3>& collisionPoints){
  Timer timer;
  bool c = false;
  if(sh.cellSize < 1.733*edt.cellSize)
  {
    //std::cout << "Cannot accurately test. Cell sizes are too similar\n";
  }

  Meshing::VolumeGrid vg;
  vg.Resize(edt.distance->m, edt.distance->n, edt.distance->p);
  vg.value = *(edt.distance);
  vg.bb = edt.aabb; 
  Meshing::VolumeGridIterator<Real> it = vg.getIterator();
  std::list<int> cellCollisions;

  // std::cout << "Before while" << "\n";
  // std::cout.flush();
  int s = sh.getTableSize();

  // std::cout << vg.value.m << "\n";
  // std::cout << vg.value.n << "\n";
  // std::cout << vg.value.p << "\n";

  /*Iterate through the EDT and hash each center point of those cells
  that are on the surface or inside the object to see if there may be a collision*/
  while(!it.isDone())
  {
    // std::cout << "Beginning while" << "\n";
    // std::cout.flush();

    //Test the value in the cell to see if it is in the mesh object
    //std::cout << *it << " ";
    if (*it <= 0)
    {
      AABB3D cell;
      int currentIndex;
      if(edt.currentTransform.isIdentity())
      {
        it.getCell(cell);
      }
      else
      {
        Box3D bcell;
        it.getCell(cell);

        bcell.setTransformed(cell, edt.currentTransform);
        bcell.getAABB(cell);
      }
      //Test the bbmin and bbmax of the cell to see if they hash to the same value
      if (sh.Hash(cell.bmin) == sh.Hash(cell.bmax))
      {
        //Test to see if the cell is in an empty hash bucket
        if(sh.currentBox.contains(cell.bmin))
        {
          currentIndex = sh.Hash(cell.bmin);
          if(!sh.cellEmpty(currentIndex))
          { 
            //std::cout << sh.Hash(cell.bmin) << " ";
            cellCollisions.push_back(currentIndex);
              // int x;
              // std::cout << sh.Hash(cell.bmin) << " ";
              // std::cin >> x;
          }
        }
      }else{
        Vector3 currentPt;
        //Test each of the vertices of the cell
        if(sh.currentBox.contains(cell.bmin))
        {
          currentIndex = sh.Hash(cell.bmin);
          if(!sh.cellEmpty(currentIndex))
          {
            // std::cout << "(" << cell.bmin.x << ", " << cell.bmin.y << ", " << cell.bmin.z << ") ";
            // std::cout << sh.Hash(cell.bmin) << " "
              cellCollisions.push_back(currentIndex);
              // int x;
              // std::cout << sh.Hash(cell.bmin) << " ";
              // std::cout << "Inside the if.....................";
              // std::cin >> x;            
          }
        }
        currentPt.set(cell.bmin.x, cell.bmin.y, cell.bmax.z);
        if(sh.currentBox.contains(currentPt))
        {
          currentIndex = sh.Hash(cell.bmin.x, cell.bmin.y, cell.bmax.z);
          if(!sh.cellEmpty(currentIndex))
          {            
              cellCollisions.push_back(currentIndex);
              // int x;
              // std::cout << sh.Hash(cell.bmin.x, cell.bmin.y, cell.bmax.z) << " ";
              // std::cin >> x;            
          }
        }
        currentPt.set(cell.bmin.x, cell.bmax.y, cell.bmin.z);
        if(sh.currentBox.contains(currentPt))
        {
          currentIndex = sh.Hash(cell.bmin.x, cell.bmax.y, cell.bmin.z);
          if(!sh.cellEmpty(currentIndex))
          {
            cellCollisions.push_back(currentIndex);
            //   int x;
            //   std::cout << sh.Hash(cell.bmin.x, cell.bmax.y, cell.bmin.z) << " ";
            //   std::cin >> x;             
          }
        }
        currentPt.set(cell.bmin.x, cell.bmax.y, cell.bmax.z);
        if(sh.currentBox.contains(currentPt))
        {
          currentIndex = sh.Hash(cell.bmin.x, cell.bmax.y, cell.bmax.z);
          if(!sh.cellEmpty(currentIndex))
          {
            cellCollisions.push_back(currentIndex);
            // int x;
            // std::cout << sh.Hash(cell.bmin.x, cell.bmax.y, cell.bmax.z) << " ";
            // std::cin >> x;            
          }
        }
        currentPt.set(cell.bmax.x, cell.bmin.y, cell.bmin.z);
        if(sh.currentBox.contains(currentPt))
        {
          currentIndex = sh.Hash(cell.bmax.x, cell.bmin.y, cell.bmin.z);                                  
          if(!sh.cellEmpty(currentIndex))
          {
            cellCollisions.push_back(currentIndex);
            // int x;
            // std::cout << sh.Hash(cell.bmax.x, cell.bmin.y, cell.bmin.z) << " ";
            // std::cin >> x;
          }
        }
        currentPt.set(cell.bmax.x, cell.bmin.y, cell.bmax.z);
        if(sh.currentBox.contains(currentPt))
        {
          currentIndex = sh.Hash(cell.bmax.x, cell.bmin.y, cell.bmax.z);
          if(!sh.cellEmpty(currentIndex))
          {
            cellCollisions.push_back(currentIndex);
            // int x;
            // std::cout << sh.Hash(cell.bmax.x, cell.bmin.y, cell.bmax.z) << " ";
            // std::cin >> x;
          }
        }
        currentPt.set(cell.bmax.x, cell.bmax.y, cell.bmin.z);
        if(sh.currentBox.contains(currentPt))
        {
          currentIndex = sh.Hash(cell.bmax.x, cell.bmax.y, cell.bmin.z);
          if(!sh.cellEmpty(currentIndex))
          {
            cellCollisions.push_back(currentIndex);
            // int x;
            // std::cout << sh.Hash(cell.bmax.x, cell.bmax.y, cell.bmin.z) << " ";
            // std::cin >> x;
          }
        }
        if(sh.currentBox.contains(cell.bmax))
        {
          currentIndex = sh.Hash(cell.bmax);
          if(!sh.cellEmpty(currentIndex))
          {            
            cellCollisions.push_back(currentIndex);
            // int x;
            // std::cout << sh.Hash(cell.bmax) << " ";
            // std::cin >> x;
          }
        }
      }
    }
    ++it;
  }

  //std::cout << "Out of while" << "\n";
  //std::cout << cellCollisions.size() << "\n";
  //int BucketCounter = 0;
  int PtCounter = 0;
  if (!cellCollisions.empty())
  {
    //std::cout << "cellCollisions not empty\n";
    int i, j, k;
    Vector3 transPoint;
    //std::cout << "Into if" << "\n";  
    cellCollisions.sort();
    cellCollisions.unique();
    //std::cout << cellCollisions.size() << "\n";

    for(std::list<int>::iterator bucketIt = cellCollisions.begin(); bucketIt != cellCollisions.end(); bucketIt++){
      //std::cout << *bucketIt << "\n";
    }

    //std::cout << "Made unique" << "\n";
    for(std::list<int>::iterator bucketIt = cellCollisions.begin(); bucketIt != cellCollisions.end(); bucketIt++)
    {
      //std::cout << "Bucket: " << BucketCounter << "\n";
      //std::cout << "Took next bucket *****************************************" << "\n";
      // std::cout << *bucketIt << " ";
      // std::cout << sh.table[*bucketIt].size() << "\n";
      //if (sh.table[*bucketIt].size() < 10000){
        //std::cout << "Bucket has right size";
        for(std::vector<Vector3>::iterator ptIt = sh.table[*bucketIt].begin(); ptIt != sh.table[*bucketIt].end(); ptIt++)
        {
          //std::cout << PtCounter << "\n";
          //std::cout << "Getting index for point: " << ptIt->x << " " << ptIt->y << " " << ptIt->z << "\n";
          sh.currentTransform.mulPoint(*ptIt, transPoint);
          if(edt.bb.contains(transPoint))
          {
            edt.currentTransform.mulPointInverse(transPoint,transPoint);
            vg.GetIndex(transPoint,i,j,k);

            //std::cout << "Got index\n";

            //std::cout << "Index:" << i << " " << j << " " << k << "\n";
            //if(i>=0 && i<=vg.value.m && j>=0 && j<=vg.value.n && k>=0 && k<=vg.value.p)          
            if(vg.value(i,j,k) <= 0)
            {
              //Disregard list of points for now
              return true;
              // std::cout << "Index: " << i << " " << j << " " << k << "\n";
              // std::cout << "Point: " << ptIt->x << " " << ptIt->y << " " << ptIt->z << "\n";
              collisionPoints.push_back(transPoint);         
            }
          }
        }
       
      //std::cout << "Finished bucket" << "\n"; 
    }
  } 

  if(!collisionPoints.empty())
  {
    c = true;
    // for(std::vector<Meshing::Vector3>::iterator it = collisionPoints.begin(); it != collisionPoints.end(); it++)
    // {
    //   //std::cout << it->x << ", " << it->y << ", " << it->z << "\n";
    //   PtCounter++;
    // }
    // std::cout << PtCounter << "\n";
  }
  // std::cout << "done" << "\n";
  // std::cout.flush();

  return c;
}

//Function to collide a pointcloud and mesh by do a brute-force comparison
//  of the points in the cloud to the EDT
bool Collide(CollisionPointCloud& pc, EuclideanDistanceTransform& edt)
{
  Vector3 transPoint;
  int i,j,k;

  Meshing::VolumeGrid vg;
  vg.Resize(edt.distance->m, edt.distance->n, edt.distance->p);
  vg.value = *(edt.distance);
  vg.bb = edt.aabb; 
  
  for(std::vector<Vector3>::iterator it = pc.points.begin(); it!=pc.points.end(); it++)
  {
    pc.currentTransform.mul(*it, transPoint);
    if(edt.bb.contains(transPoint))
    {
      edt.currentTransform.mulPointInverse(transPoint,transPoint);
      vg.GetIndex(transPoint,i,j,k);  

      if(vg.value(i,j,k) <= 0)
      {
        return true;
      }
    }
  }
  return false;
}

} //namespace Geometry
  