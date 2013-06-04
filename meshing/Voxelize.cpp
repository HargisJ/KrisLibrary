#include "Voxelize.h"
#include "Rasterize.h"
#include "ClosestPoint.h"
#include "VolumeGrid.h"
#include <structs/FixedSizeHeap.h>
#include <structs/Heap.h>
#include <math/cast.h>
#include <math3d/Plane3D.h>
#include <math3d/Segment3D.h>
#include <geometry/primitives.h>
#include <math/random.h>
#include <Timer.h>
#include <list>
#include <set>
using namespace Geometry;

namespace Meshing {


void FitGridToMesh(int m,int n,int p,AABB3D& bb,const TriMesh& mesh)
{
  Assert(m>0 && n>0 && p>0);
  AABB3D mbb;
  mesh.GetAABB(mbb.bmin,mbb.bmax);
  Vector3 size=mbb.bmax-mbb.bmin;
  Vector3 cell;
  if(m > 2) cell.x = size.x/(m-2);
  else if(m == 2) cell.x = size.x;
  else if(m == 1) cell.x = size.x*2.0;
  if(n > 2) cell.y = size.y/(n-2);
  else if(n == 2) cell.y = size.y;
  else if(n == 1) cell.y = size.y*2.0;
  if(p > 2) cell.z = size.z/(p-2);
  else if(p == 2) cell.z = size.z;
  else if(p == 1) cell.z = size.z*2.0;
  size.x = cell.x*m;
  size.y = cell.y*n;
  size.z = cell.z*p;
  Vector3 center = (mbb.bmax+mbb.bmin)*0.5;
  bb.bmin = center - 0.5*size;
  bb.bmax = center + 0.5*size;
}

template <class T>
void FitGridToMesh(const Array3D<T>& cells,AABB3D& bb,const TriMesh& m)
{
  FitGridToMesh(cells.m,cells.n,cells.p,bb,m);
}

template <class T>
void GetGridCell(const Array3D<T>& cells,const AABB3D& bb,const IntTriple& index,AABB3D& cell)
{
  cell.bmin.x=bb.bmin.x + Real(index.a) / Real(cells.m) *(bb.bmax.x-bb.bmin.x);
  cell.bmin.y=bb.bmin.y + Real(index.b) / Real(cells.n) *(bb.bmax.y-bb.bmin.y);
  cell.bmin.z=bb.bmin.z + Real(index.c) / Real(cells.p) *(bb.bmax.z-bb.bmin.z);
  cell.bmax.x=bb.bmin.x + Real(index.a+1) / Real(cells.m) *(bb.bmax.x-bb.bmin.x);
  cell.bmax.y=bb.bmin.y + Real(index.b+1) / Real(cells.n) *(bb.bmax.y-bb.bmin.y);
  cell.bmax.z=bb.bmin.z + Real(index.c+1) / Real(cells.p) *(bb.bmax.z-bb.bmin.z);
}

template <class T>
void GetGridCellCenter(const Array3D<T>& cells,const AABB3D& bb,const IntTriple& index,Vector3& c)
{
  c.x=bb.bmin.x + (Real(index.a)+Half) / Real(cells.m) *(bb.bmax.x-bb.bmin.x);
  c.y=bb.bmin.y + (Real(index.b)+Half) / Real(cells.n) *(bb.bmax.y-bb.bmin.y);
  c.z=bb.bmin.z + (Real(index.c)+Half) / Real(cells.p) *(bb.bmax.z-bb.bmin.z);
}

bool QueryGrid(int m,int n,int p,const AABB3D& grid,const Vector3& query,IntTriple& index)
{
  Vector3 param;
  Vector3 size = grid.bmax-grid.bmin;
  param = (query - grid.bmin);
  if(param.x < 0 || param.x >= size.x) return false;
  if(param.y < 0 || param.y >= size.y) return false;
  if(param.z < 0 || param.z >= size.z) return false;
  param.x /= size.x;
  param.y /= size.y;
  param.z /= size.z;
  index.a = (int)iFloor(param.x*m);
  index.b = (int)iFloor(param.y*n);
  index.c = (int)iFloor(param.z*p);
  return true;
}


template <class T>
inline bool QueryGrid(const Array3D<T>& cells,const AABB3D& grid,const Vector3& query,IntTriple& index)
{
  return QueryGrid(cells.m,cells.n,cells.p,grid,query,index);
}

bool QueryGrid(int m,int n,int p,const AABB3D& grid,const AABB3D& query,IntTriple& low,IntTriple& high)
{
  if(!grid.intersects(query)) return false;
  Vector3 paramLow,paramHigh;
  Vector3 size = grid.bmax-grid.bmin;
  Assert(query.bmax.x >= query.bmin.x);
  Assert(query.bmax.y >= query.bmin.y);
  Assert(query.bmax.z >= query.bmin.z);
  paramLow = (query.bmin - grid.bmin);
  paramHigh = (query.bmax - grid.bmin);
  paramLow.x /= size.x;
  paramLow.y /= size.y;
  paramLow.z /= size.z;
  paramHigh.x /= size.x;
  paramHigh.y /= size.y;
  paramHigh.z /= size.z;
  low.a = (int)iFloor(paramLow.x*m);
  low.b = (int)iFloor(paramLow.y*n);
  low.c = (int)iFloor(paramLow.z*p);
  high.a = (int)iFloor(paramHigh.x*m);
  high.b = (int)iFloor(paramHigh.y*n);
  high.c = (int)iFloor(paramHigh.z*p);
  if(low.a < 0) low.a=0;
  if(low.b < 0) low.b=0;
  if(low.c < 0) low.c=0;
  if(high.a >= m) high.a=m-1;
  if(high.b >= n) high.b=n-1;
  if(high.c >= p) high.c=p-1;
  //maybe some numerical error could cause this
  if(low.a >= m) return false;
  if(low.b >= n) return false;
  if(low.c >= p) return false;
  if(high.a < 0) return false;
  if(high.b < 0) return false;
  if(high.c < 0) return false;
  return true;
}

template <class T>
inline bool QueryGrid(const Array3D<T>& cells,const AABB3D& grid,const AABB3D& query,IntTriple& low,IntTriple& high)
{
  return QueryGrid(cells.m,cells.n,cells.p,grid,query,low,high);
}

//estimates the fraction of the cell on the negative side of the plane
Real GridCellDensity(const AABB3D& cell,const Plane3D& p)
{
  Real dmin,dmax;
  p.distanceLimits(cell,dmin,dmax);
  if(dmin >= 0) return 0;
  if(dmax <= 0) return 1;
  Vector3 ptOnPlane = p.offset*p.normal;
  Vector3 planeNormal = p.normal;
  //scale and translate pt and plane so that cell is a [-1,1]^3 cube
  ptOnPlane -= (cell.bmin+cell.bmax)*0.5;
  Vector3 scale;
  scale.x = 2.0/(cell.bmax.x-cell.bmin.x);
  scale.y = 2.0/(cell.bmax.y-cell.bmin.y);
  scale.z = 2.0/(cell.bmax.z-cell.bmin.z);
  planeNormal.x /= scale.x;
  planeNormal.y /= scale.y;
  planeNormal.z /= scale.z;
  planeNormal.inplaceNormalize();
  ptOnPlane.x *= scale.x;
  ptOnPlane.y *= scale.y;
  ptOnPlane.z *= scale.z;
  //estimate by approximating with a sphere
  Real R = Sqrt(3.0);
  Real d = planeNormal.dot(ptOnPlane);
  if(d < -R || d > R) {
    cerr<<"Warning, numerical error in GridCellDensity"<<endl;
    cerr<<"   point "<<ptOnPlane<<endl;
    cerr<<"   d="<<d<<", R="<<R<<endl;
    getchar();
    if(d < -R) return 0;
    else return 1;
  }
  Assert(d >= -R && d <= R);
  Real R3 = R*R*R;
  Real V = 4.0/3.0*R3;
  Real Vd= (2.0/3.0*R3 + d*R*R - d*d*d/3.0);
  Assert(Vd >= 0 &&  Vd <= V);
  return Vd/V;
}



void GetTriangleBuckets(const TriMesh& m,const AABB3D& bb,Array3D<list<int> >& triangles)
{
  triangles.set(list<int>());
  Triangle3D tri;
  AABB3D query,cell;
  IntTriple hi,lo;
  for(size_t i=0;i<m.tris.size();i++) {
    m.GetTriangle(i,tri);
    query.setPoint(tri.a);
    query.expand(tri.b);
    query.expand(tri.c);
    bool q=QueryGrid(triangles,bb,query,lo,hi);
    if(!q) continue;
    //getchar();
    IntTriple index=lo;
    while(index.a <= hi.a) {
      while(index.b <= hi.b) {
	while(index.c <= hi.c) {
	  GetGridCell(triangles,bb,index,cell);
	  if(tri.intersects(cell)) {
	    triangles(index.a,index.b,index.c).push_back(i);
	  }
	  index.c++;
	}
	index.b++;
	index.c = lo.c;
      }
      index.a++;
      index.b = lo.b;
    }
  }
}

void GetSegmentCells(const Segment3D& s,vector<IntTriple>& cells)
{
  cells.resize(0);
  Vector3 d=s.b-s.a;
  IntTriple i;
  i.a = (int)iFloor(s.a.x);
  i.b = (int)iFloor(s.a.y);
  i.c = (int)iFloor(s.a.z);

  Vector3 cellCorner(i.a,i.b,i.c);
  Real invdx=1.0/d.x;
  Real invdy=1.0/d.y;
  Real invdz=1.0/d.z;
  Real param=0;  //goes from 0 to 1
  while(param < 1) {
    cells.push_back(i);
    //see which cell face is hit next
    int closest=0;  //1=+x,-1=-x,2=+y,-2=-y,3=+z,-3=-z
    param=Inf;
    if(d.x > 0) { //check +x face
      param = ((cellCorner.x+1.0)-s.a.x)*invdx;
      closest = 1;
    }
    else if(d.x < 0) { //check -x face
      param = (cellCorner.x-s.a.x)*invdx;
      closest = -1;
    }
    if(d.y > 0) { //check +y face
      if((cellCorner.y+1.0)-s.a.y < param*d.y) {
	param = ((cellCorner.y+1.0)-s.a.y)*invdy;
	closest = 2;
      }
    }
    else if(d.y < 0) { //check -y face
      if(cellCorner.y-s.a.y > param*d.y) {
	param = (cellCorner.y-s.a.y)*invdy;
	closest = -2;
      }
    }
    if(d.z > 0) { //check +z face
      if((cellCorner.z+1.0)-s.a.z < param*d.z)
	param = ((cellCorner.z+1.0)-s.a.z)*invdz; {
	closest = 3;
      }
    }
    else if(d.z < 0) { //check -y face
      if(cellCorner.z-s.a.z > param*d.z) {
	param = (cellCorner.z-s.a.z)*invdz;
	closest = -3;
      }
    }
    switch(closest) {
    case 0: param = 1; break;
    case 1: i.a++; cellCorner.x+=1.0; break;
    case -1: i.a--; cellCorner.x-=1.0; break;
    case 2: i.b++; cellCorner.y+=1.0; break;
    case -2: i.b--; cellCorner.y-=1.0; break;
    case 3: i.c++; cellCorner.z+=1.0; break;
    case -3: i.c--; cellCorner.z-=1.0; break;
    }
  }
}

void GetSegmentCells(const Segment3D& s,int m,int n,int p,const AABB3D& bb,vector<IntTriple>& cells)
{
  cells.resize(0);

  Real tmin,tmax;
  if(!s.intersects(bb,tmin,tmax)) return;
  Vector3 d = s.b - s.a;

  IntTriple i;
  Vector3 aparam = (s.a+tmin*d - bb.bmin);
  aparam.x /= bb.bmax.x-bb.bmin.x;
  aparam.y /= bb.bmax.y-bb.bmin.y;
  aparam.z /= bb.bmax.z-bb.bmin.z;
  i.a = (int)iFloor(aparam.x*m);
  i.b = (int)iFloor(aparam.y*n);
  i.c = (int)iFloor(aparam.z*p);
  if(i.a < 0) i.a=0;
  if(i.b < 0) i.b=0;
  if(i.c < 0) i.c=0;
  if(i.a >= m) i.a=m-1;
  if(i.b >= n) i.b=n-1;
  if(i.c >= p) i.c=p-1;
  
  Vector3 cellCorner;
  cellCorner.x = bb.bmin.x+Real(i.a)/Real(m)*(bb.bmax.x-bb.bmin.x);
  cellCorner.y = bb.bmin.y+Real(i.b)/Real(n)*(bb.bmax.y-bb.bmin.y);
  cellCorner.z = bb.bmin.z+Real(i.c)/Real(p)*(bb.bmax.z-bb.bmin.z);
  Vector3 cellSize = bb.bmax-bb.bmin;
  cellSize.x /= m;
  cellSize.y /= n;
  cellSize.z /= p;
  Real param=tmin;
  while(param < tmax) {
    cells.push_back(i);
    //see which cell face is hit next
    int closest=0;  //1=+x,-1=-x,2=+y,-2=-y,3=+z,-3=-z
    Real param=Inf;
    if(d.x > 0) { //check +x face
      param = ((cellCorner.x+cellSize.x)-s.a.x)/d.x;
      closest = 1;
    }
    else if(d.x < 0) { //check -x face
      param = (cellCorner.x-s.a.x)/d.x;
      closest = -1;
    }
    if(d.y > 0) { //check +y face
      if((cellCorner.y+cellSize.y)-s.a.y < param*d.y)
	param = ((cellCorner.y+cellSize.y)-s.a.y)/d.y;
      closest = 2;
    }
    else if(d.y < 0) { //check -y face
      if(cellCorner.y-s.a.y > param*d.y)
	param = (cellCorner.y-s.a.y)/d.y;
      closest = -2;
    }
    if(d.z > 0) { //check +z face
      if((cellCorner.z+cellSize.z)-s.a.z < param*d.z)
	param = ((cellCorner.z+cellSize.z)-s.a.z)/d.z;
      closest = 3;
    }
    else if(d.z < 0) { //check -y face
      if(cellCorner.z-s.a.z > param*d.z)
	param = (cellCorner.z-s.a.z)/d.z;
      closest = -3;
    }
    switch(closest) {
    case 0: param = 1; break;
    case 1: i.a++; cellCorner.x+=cellSize.x; break;
    case -1: i.a--; cellCorner.x-=cellSize.x; break;
    case 2: i.b++; cellCorner.y+=cellSize.y; break;
    case -2: i.b--; cellCorner.y-=cellSize.y; break;
    case 3: i.c++; cellCorner.z+=cellSize.z; break;
    case -3: i.c--; cellCorner.z-=cellSize.z; break;
    }
    if(i.a < 0 || i.a >= m) break;
    if(i.b < 0 || i.b >= n) break;
    if(i.c < 0 || i.c >= p) break;
  }
}


void GetTriangleCells(const Triangle3D& tri,vector<IntTriple>& cells)
{
  cells.resize(0);
  AABB3D query;
  query.setPoint(tri.a);
  query.expand(tri.b);
  query.expand(tri.c);
  IntTriple lo,hi;
  lo.a = (int)iFloor(query.bmin.x);
  lo.b = (int)iFloor(query.bmin.y);
  lo.c = (int)iFloor(query.bmin.z);
  hi.a = (int)iFloor(query.bmax.x);
  hi.b = (int)iFloor(query.bmax.y);
  hi.c = (int)iFloor(query.bmax.z);

  AABB3D cell;
  Vector3 cellCorner;
  Vector3 cellSize(1,1,1);
  cellCorner.x=lo.a;
  IntTriple index=lo;
  while(index.a <= hi.a) {
    index.b = lo.b;
    cellCorner.y=lo.b;
    while(index.b <= hi.b) {
      index.c = lo.c;
      cellCorner.z=lo.c;
      while(index.c <= hi.c) {
	cell.bmin = cellCorner;
	cell.bmax = cellCorner + cellSize;
	if(tri.intersects(cell)) cells.push_back(index);
	
	index.c++;
	cellCorner.z += cellSize.z;
      }
      index.b++;
      cellCorner.y += cellSize.y;
    }
    index.a++;
    cellCorner.x += cellSize.x;
  }
}

void GetTriangleCells(const Triangle3D& tri,int m,int n,int p,const AABB3D& bb,vector<IntTriple>& cells)
{
  cells.resize(0);
  AABB3D query;
  query.setPoint(tri.a);
  query.expand(tri.b);
  query.expand(tri.c);
  IntTriple lo,hi;
  bool q=QueryGrid(m,n,p,bb,query,lo,hi);
  if(!q) return;

  AABB3D cell;
  Vector3 cellCorner;
  Vector3 cellSize = bb.bmax-bb.bmin;
  cellSize.x /= m;
  cellSize.y /= n;
  cellSize.z /= p;
  cellCorner.x = bb.bmin.x+Real(lo.a)/Real(m)*(bb.bmax.x-bb.bmin.x);
  IntTriple index=lo;
  while(index.a <= hi.a) {
    index.b = lo.b;
    cellCorner.y = bb.bmin.y+Real(lo.b)/Real(n)*(bb.bmax.y-bb.bmin.y);
    while(index.b <= hi.b) {
      index.c = lo.c;
      cellCorner.z = bb.bmin.z+Real(lo.c)/Real(p)*(bb.bmax.z-bb.bmin.z);
      while(index.c <= hi.c) {
	cell.bmin = cellCorner;
	cell.bmax = cellCorner + cellSize;
	if(tri.intersects(cell)) cells.push_back(index);
	
	index.c++;
	cellCorner.z += cellSize.z;
      }
      index.b++;
      cellCorner.y += cellSize.y;
    }
    index.a++;
    cellCorner.x += cellSize.x;
  }

}

template <class T>
void RasterizeXYSegment(const Segment3D& s,int i,int j,const Array3D<T>& grid,const AABB3D& bb,vector<IntTriple>& cells)
{
  Real wa = (s.a.z-bb.bmin.z)/(bb.bmax.z-bb.bmin.z);
  Real wb = (s.b.z-bb.bmin.z)/(bb.bmax.z-bb.bmin.z);
  int kmin = (int)iFloor(wa*grid.p);
  int kmax = (int)iFloor(wb*grid.p);
  if(kmin > kmax) std::swap(kmin,kmax);
  for(int k=Max(0,kmin);k<Min(kmax+1,grid.p);k++) 
    cells.push_back(IntTriple(i,j,k));
}

void RasterizeXQuadrilateral(const Segment3D& s1,const Segment3D& s2,int i,int n,int p,vector<IntTriple>& cells)
{
  Triangle2D t;
  vector<IntPair> xcells;
  t.a.set(s1.a.y,s1.a.z);
  t.b.set(s1.b.y,s1.b.z);
  t.c.set(s2.a.y,s2.a.z);
  GetTriangleCells_Clipped(t,xcells,0,0,n,p);
  for(size_t c=0;c<xcells.size();c++) 
    cells.push_back(IntTriple(i,xcells[c].a,xcells[c].b));

  t.a.set(s1.b.y,s1.b.z);
  t.b.set(s2.b.y,s2.b.z);
  t.c.set(s2.a.y,s2.a.z);
  GetTriangleCells_Clipped(t,xcells,0,0,n,p);

  for(size_t c=0;c<xcells.size();c++)
    cells.push_back(IntTriple(i,xcells[c].a,xcells[c].b));

  //TODO: remove duplicates
}


void GetTriangleCells2(const Triangle3D& torig,int m,int n,int p,const AABB3D& bb,vector<IntTriple>& cells)
{
  Vector3 a,b,c;
  //rotate vertices so that a is max, b is next, then c
  //(normal order is unnecessary)
  //get a to the top
  if(Lexical3DOrder(torig.a,torig.b)) { //a < b
    if(Lexical3DOrder(torig.b,torig.c)) { //a < b < c
      a = torig.c;
      b = torig.b;
      c = torig.a;
    }
    else {  //a < b, c < b
      if(Lexical3DOrder(torig.a,torig.c)) { //a < c < b 
	a = torig.b;
	b = torig.c;
	c = torig.a;
      }
      else { //c < a < b
	a = torig.b;
	b = torig.a;
	c = torig.c;
      }
    }
  }
  else if(Lexical3DOrder(torig.a,torig.c)) { //b < a < c
    a = torig.c;
    b = torig.a;
    c = torig.b;
  }
  else {  //c < a, b < a
    if(Lexical3DOrder(torig.b,torig.c)) { //b < c < a
      a = torig.a;
      b = torig.c;
      c = torig.b;
    }
    else {  //c < b < a
      a = torig.a;
      b = torig.b;
      c = torig.c;
    }
  }

  //convert to standard frame
  a -= bb.bmin;
  b -= bb.bmin;
  c -= bb.bmin;
  a.x *= Real(m)/(bb.bmax.x-bb.bmin.x);
  b.x *= Real(m)/(bb.bmax.x-bb.bmin.x);
  c.x *= Real(m)/(bb.bmax.x-bb.bmin.x);
  a.y *= Real(n)/(bb.bmax.y-bb.bmin.y);
  b.y *= Real(n)/(bb.bmax.y-bb.bmin.y);
  c.y *= Real(n)/(bb.bmax.y-bb.bmin.y);
  a.z *= Real(p)/(bb.bmax.z-bb.bmin.z);
  b.z *= Real(p)/(bb.bmax.z-bb.bmin.z);
  c.z *= Real(p)/(bb.bmax.z-bb.bmin.z);

  int imax = (int)iFloor(a.x);
  int imid = (int)iFloor(b.x);
  int imin = (int)iFloor(c.x);
  Real uca_halfway;  //halfway point interpolating between segment c->a
  if(a.x == c.x) uca_halfway = 1;
  else uca_halfway = (b.x-c.x)/(a.x-c.x);  
  Real cbstep = 1.0/(b.x - c.x);
  Real castep = 1.0/(a.x - c.x);
  Real bastep = 1.0/(a.x - b.x);
  //c to b
  cells.resize(0);
  for(int i=Max(imin,0);i<Min(imid+1,m);i++) {
    //get uz range of slice at k
    Real xmin = Real(i);
    Real xmax = xmin + 1.0;
    //find the range of the triangle between xmin and xmax, voxelize it
    Segment3D segmin;
    Segment3D segmax;
    Assert(xmin <= b.x);
    Assert(xmax >= c.x);
    //these are the ranges of u parameters between xmin,xmax for cb,ca
    Real ucb1,ucb2,uca1,uca2;
    if(b.x == c.x) {
      ucb1 = 0;
      ucb2 = 1;
    }
    else {
      ucb1 = (xmin - c.x)*cbstep;
      ucb2 = ucb1 + cbstep;
      if(ucb1 < 0) ucb1 = 0;
      if(ucb2 > 1) ucb2 = 1;
    }
    if(a.x == c.x) {
      uca1 = 0;
      uca2 = uca_halfway;
    }
    else {
      uca1 = (xmin - c.x)*castep;
      uca2 = uca1 + castep;
      if(uca1 < 0) uca1 = 0;
      if(uca2 > uca_halfway) uca2 = uca_halfway;
    }
    segmin.a = (1-uca1)*c + uca1*a;
    segmin.b = (1-ucb1)*c + ucb1*b;

    segmax.a = (1-uca2)*c + uca2*a;
    segmax.b = (1-ucb2)*c + ucb2*b;

    //rasterize quadrilateral from segmin to segmax in y,z slice i
    RasterizeXQuadrilateral(segmin,segmax,i,n,p,cells);
  }

  //range from b to a
  for(int i=Max(imid,0);i<Min(imax+1,m);i++) {
    //get uz range of slice at k
    Real xmin = Real(i);
    Real xmax = xmin + 1.0;
    //find the range of the triangle between zmin and zmax, voxelize it
    Segment3D segmin;
    Segment3D segmax;
    Assert(xmin <= a.x);
    Assert(xmax >= b.x);
    //these are the ranges of u parameters between xmin,xmax for ba,ca
    Real uba1,uba2,uca1,uca2;
    if(a.x == b.x) {
      uba1 = 0;
      uba2 = 1;
    }
    else {
      uba1 = (xmin - b.x)*bastep;
      uba2 = uba1 + bastep;
      if(uba1 < 0) uba1 = 0;
      if(uba2 > 1) uba2 = 1;
    }
    if(a.x == c.x) {
      uca1 = uca_halfway;
      uca2 = 1;
    }
    else {
      uca1 = (xmin - c.x)*castep;
      uca2 = uca1+castep;
      if(uca1 < uca_halfway) uca1 = uca_halfway;
      if(uca2 > 1) uca2 = 1;
    }
    segmin.a = (1-uca1)*c + uca1*a;
    segmin.b = (1-uba1)*b + uba1*a;

    segmax.a = (1-uca2)*c + uca2*a;
    segmax.b = (1-uba2)*b + uba2*a;

    //rasterize quadrilateral from segmin to segmax in y,z slice i
    RasterizeXQuadrilateral(segmin,segmax,i,n,p,cells);
  }

  //Note: some cells may be revisited at imid! This may be inefficient
  //Also, RasterizeXQuadrilateral may revisit cells
  set<IntTriple> cellSet;
  for(size_t i=0;i<cells.size();i++)
    cellSet.insert(cells[i]);
  cells.resize(0);
  for(set<IntTriple>::const_iterator i=cellSet.begin();i!=cellSet.end();i++)
    cells.push_back(*i);

  /*
  vector<IntTriple> testCells;
  set<IntTriple> testSet;
  GetTriangleCells(torig,m,n,p,bb,testCells);
  for(size_t i=0;i<testCells.size();i++)
    testSet.insert(testCells[i]);
  for(size_t i=0;i<testCells.size();i++) {
    if(!cellSet.count(testCells[i])) {
      cout<<"Triangle "<<a<<" "<<b<<" "<<c<<endl;
      cerr<<"GetTriangleCells2 incorrect: doesnt contain "<<testCells[i]<<endl;
      getchar();
    }
  }
  for(size_t i=0;i<cells.size();i++) {
    if(!testSet.count(cells[i])) {
      cout<<"Triangle "<<a<<" "<<b<<" "<<c<<endl;
      cerr<<"GetTriangleCells2 incorrect: contains additional "<<cells[i]<<endl;
      getchar();
    }
  }
  */
}

template <class T>
void GetTriangleCells2(const Triangle3D& torig,const Array3D<T>& grid,const AABB3D& bb,vector<IntTriple>& cells)
{
  GetTriangleCells2(torig,grid.m,grid.n,grid.p,bb,cells);
}

void SurfaceOccupancyGrid(const TriMesh& m,Array3D<bool>& occupied,AABB3D& bb)
{
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(occupied,bb,m);
  occupied.set(false);
  Triangle3D tri;
  AABB3D query,cell;
  IntTriple hi,lo;
  for(size_t i=0;i<m.tris.size();i++) {
    m.GetTriangle(i,tri);
    query.setPoint(tri.a);
    query.expand(tri.b);
    query.expand(tri.c);
    bool q=QueryGrid(occupied,bb,query,lo,hi);
    if(!q) continue;
    //getchar();
    VolumeGridIterator<bool> it(occupied,bb);
    it.setRange(lo,hi);
    for(;!it.isDone();++it) {
      it.getCell(cell);
      if(tri.intersects(cell))
	occupied(it.index)=true;
    }
    /*
    IntTriple index=lo;
    while(index.a <= hi.a) {
      while(index.b <= hi.b) {
	while(index.c <= hi.c) {
	  if(!occupied(index)) {
	    GetGridCell(occupied,bb,index,cell);
	    if(tri.intersects(cell)) {
	      occupied(index)=true;
	    }
	  }
	  index.c++;
	}
	index.b++;
	index.c = lo.c;
      }
      index.a++;
      index.b = lo.b;
    }
    */
  }
}

void VolumeOccupancyGrid_FloodFill(const TriMesh& m,Array3D<bool>& occupied,AABB3D& bb,const IntTriple& seed,bool seedOccupied)
{
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(occupied,bb,m);

  Array3D<bool> surface(occupied.m,occupied.n,occupied.p);
  SurfaceOccupancyGrid(m,surface,bb);

  occupied.set(!seedOccupied);
  //set the surface cells to seedOccupied
  for(int i=0;i<surface.m;i++) 
    for(int j=0;j<surface.n;j++) 
      for(int k=0;k<surface.p;k++) 
	if(surface(i,j,k)) occupied(i,j,k) = seedOccupied;

  list<IntTriple> queue;
  queue.push_back(seed);
  while(!queue.empty()) {
    IntTriple index=queue.front();
    queue.pop_front();
    
    //go to neighbors
    index.a++;
    if(index.a < occupied.m && occupied(index.a,index.b,index.c) != seedOccupied) {
      occupied(index.a,index.b,index.c) = seedOccupied;
      queue.push_back(index);
    }
    index.a-=2;
    if(index.a >= 0 && occupied(index.a,index.b,index.c) != seedOccupied) {
      occupied(index.a,index.b,index.c) = seedOccupied;
      queue.push_back(index);
    }
    index.a++;

    index.b++;
    if(index.b < occupied.n && occupied(index.a,index.b,index.c) != seedOccupied) {
      occupied(index.a,index.b,index.c) = seedOccupied;
      queue.push_back(index);
    }
    index.b-=2;
    if(index.b >= 0 && occupied(index.a,index.b,index.c) != seedOccupied) {
      occupied(index.a,index.b,index.c) = seedOccupied;
      queue.push_back(index);
    }
    index.b++;

    index.c++;
    if(index.c < occupied.p && occupied(index.a,index.b,index.c) != seedOccupied) {
      occupied(index.a,index.b,index.c) = seedOccupied;
      queue.push_back(index);
    }
    index.c-=2;
    if(index.c >= 0 && occupied(index.a,index.b,index.c) != seedOccupied) {
      occupied(index.a,index.b,index.c) = seedOccupied;
      queue.push_back(index);
    }
    index.c++;
  }

  if(!seedOccupied) {  //need to set surface cells
    for(int i=0;i<surface.m;i++) 
      for(int j=0;j<surface.n;j++) 
	for(int k=0;k<surface.p;k++) 
	  if(surface(i,j,k)) occupied(i,j,k) = true;
  }
}

void DensityEstimate_FloodFill(const TriMeshWithTopology& m,Array3D<Real>& density,AABB3D& bb,const IntTriple& seed)
{
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(density,bb,m);

  Array3D<list<int> > buckets(density.m,density.n,density.p);
  GetTriangleBuckets(m,bb,buckets);

  density.set(1.0);
  //set the surface cells to 0.0
  for(int i=0;i<buckets.m;i++) 
    for(int j=0;j<buckets.n;j++) 
      for(int k=0;k<buckets.p;k++) 
	if(!buckets(i,j,k).empty()) density(i,j,k) = 0.0;

  list<IntTriple> queue;
  queue.push_back(seed);
  while(!queue.empty()) {
    IntTriple index=queue.front();
    queue.pop_front();
    
    //go to neighbors
    index.a++;
    if(index.a < density.m && density(index.a,index.b,index.c) != 0) {
      density(index.a,index.b,index.c) = 0;
      queue.push_back(index);
    }
    index.a-=2;
    if(index.a >= 0 && density(index.a,index.b,index.c) != 0) {
      density(index.a,index.b,index.c) = 0;
      queue.push_back(index);
    }
    index.a++;

    index.b++;
    if(index.b < density.n && density(index.a,index.b,index.c) != 0) {
      density(index.a,index.b,index.c) = 0;
      queue.push_back(index);
    }
    index.b-=2;
    if(index.b >= 0 && density(index.a,index.b,index.c) != 0) {
      density(index.a,index.b,index.c) = 0;
      queue.push_back(index);
    }
    index.b++;

    index.c++;
    if(index.c < density.p && density(index.a,index.b,index.c) != 0) {
      density(index.a,index.b,index.c) = 0;
      queue.push_back(index);
    }
    index.c-=2;
    if(index.c >= 0 && density(index.a,index.b,index.c) != 0) {
      density(index.a,index.b,index.c) = 0;
      queue.push_back(index);
    }
    index.c++;
  }

  //estimate the densities of the surface cells
  Triangle3D tri;
  for(int i=0;i<buckets.m;i++) 
    for(int j=0;j<buckets.n;j++) 
      for(int k=0;k<buckets.p;k++) {
	if(!buckets(i,j,k).empty()) {
	  const list<int>& b=buckets(i,j,k);
	  Real minDistance = Inf;
	  Vector3 minDir;
	  AABB3D cell;
	  Vector3 cellCenter;
	  GetGridCell(buckets,bb,IntTriple(i,j,k),cell);
	  cellCenter = 0.5*(cell.bmin+cell.bmax);
	  for(list<int>::const_iterator a=b.begin();a!=b.end();a++) {
	    m.GetTriangle(*a,tri);
	    Vector3 cp = tri.closestPoint(cellCenter);
	    Vector3 dir = cellCenter-cp;
	    Real dist = dir.norm();
	    if(dir.dot(tri.normal()) > 0) //on inside
	      dist = -dist;
	    if(dist < minDistance) {
	      minDistance = dist;
	      minDir = dir/dist;
	    }
	  }
	  Plane3D p;
	  p.offset = minDir.dot(cellCenter) - minDistance;
	  p.normal = minDir;
	  density(i,j,k) = GridCellDensity(cell,p);
	}
      }

}

void VolumeOccupancyGrid_CenterShooting(const TriMesh& m,Array3D<bool>& occupied,AABB3D& bb,int shootDirection)
{
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(occupied,bb,m);

  if(shootDirection != 0) FatalError("VolumeOccupancyGrid: can only do x direction now");

  Ray3D ray;
  ray.direction.setZero();
  ray.direction[shootDirection] = 1;

  Array3D<list<int> > triangleLists(occupied.m,occupied.n,occupied.p);
  GetTriangleBuckets(m,bb,triangleLists);
  
  Real xscale = (bb.bmax.x-bb.bmin.x)/Real(occupied.m);
  Triangle3D tri;
  for(int j=0;j<occupied.n;j++) {
    for(int k=0;k<occupied.p;k++) {
      GetGridCellCenter(occupied,bb,IntTriple(0,j,k),ray.source);
      ray.source.x = bb.bmin.x;
      //parity is 1 if it's outside, -1 if inside
      int parity = 1;
      int numHits = 0;
      //sweep a ray across x axis
      for(int i=0;i<occupied.m;i++) {
	Real startdistance = (Real(i))*xscale;
	Real enddistance = (Real(i+1))*xscale;
	Real centerdistance = (Real(i)+0.5)*xscale;
	const list<int>& triangles = triangleLists(i,j,k);
	Heap<Triangle3D,Real> entering,exiting;
	for(list<int>::const_iterator a=triangles.begin();a!=triangles.end();a++) {
	  m.GetTriangle(*a,tri);
	  Real t,u,v;
	  if(tri.rayIntersects(ray,&t,&u,&v)) {
	    if(t >= startdistance && t < enddistance) {
	      if(t > centerdistance) exiting.push(tri,-t);
	      else entering.push(tri,-t);
	    }
	  }
	}
	//take all the triangles from the entering cell face to the center
	while(!entering.empty()) {
	  tri = entering.top(); entering.pop();
	  Vector3 n=tri.normal();
	  if(n[shootDirection] < 0) { //entering into the triangle
	    if(parity > 0)   //on outside
	      parity = -parity;
	  }
	  else if(n[shootDirection] > 0) {  //exiting the triangle
	    if(parity < 0)   //on inside
	      parity = -parity;
	    else if(numHits==0) {  //first hit is on inside of triangle
	      for(int i2=0;i2<i;i2++)
		occupied(i2,j,k) = true;
	    }
	  }
	  numHits++;
	}

	//set the occupation grid
	occupied(i,j,k) = (parity < 0);

	//take all the triangles from the cell center to the exiting face
	while(!exiting.empty()) {
	  tri = exiting.top(); exiting.pop();
	  Vector3 n=tri.normal();
	  if(n[shootDirection] < 0) { //entering into the triangle
	    if(parity > 0)   //on outside
	      parity = -parity;
	  }
	  else if(n[shootDirection] > 0) {  //exiting the triangle
	    if(parity < 0)   //on inside
	      parity = -parity;
	    else if(numHits==0) {  //first hit is on inside of triangle
	      for(int i2=0;i2<=i;i2++)
		occupied(i2,j,k) = true;
	    }
	  }
	  numHits++;
	}
      }
    }
  }
}

void DensityEstimate_CenterShooting(const TriMesh& m,Array3D<Real>& density,AABB3D& bb,int shootDirection)
{
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(density,bb,m);

  if(shootDirection != 0) FatalError("DensityEstimate_RandomShooting: can only do x direction now");

  Ray3D ray;
  ray.direction.setZero();
  ray.direction[shootDirection] = 1;

  Array3D<list<int> > triangleLists(density.m,density.n,density.p);
  GetTriangleBuckets(m,bb,triangleLists);
  
  density.set(Zero);
  Real xscale = (bb.bmax.x-bb.bmin.x)/Real(density.m);
  Real xscaleInv = (Real)density.m/(bb.bmax.x-bb.bmin.x);
  Triangle3D tri;
  for(int j=0;j<density.n;j++) {
    for(int k=0;k<density.p;k++) {
      GetGridCellCenter(density,bb,IntTriple(0,j,k),ray.source);
      ray.source.x = bb.bmin.x;

      //parity is 1 if it's outside, -1 if inside
      int parity = 1;
      int numHits = 0;
      //sweep a ray across x axis
      for(int i=0;i<density.m;i++) {
	Real startdistance = (Real(i))*xscale;
	Real enddistance = (Real(i+1))*xscale;
	const list<int>& triangles = triangleLists(i,j,k);
	Heap<Triangle3D,Real> tris;
	for(list<int>::const_iterator a=triangles.begin();a!=triangles.end();a++) {
	  m.GetTriangle(*a,tri);
	  Real t,u,v;
	  if(tri.rayIntersects(ray,&t,&u,&v)) {
	    if(t >= startdistance && t < enddistance)
	      tris.push(tri,-t);
	  }
	}
	//integrate the length of segment intersected by the triangle interiors
	Real occupation=0;
	Real lastt = startdistance; //last intersection into the object
	while(!tris.empty()) {
	  tri = tris.top();
	  Real t=-tris.topPriority();
	  tris.pop();
	  
	  Vector3 n=tri.normal();
	  if(n[shootDirection] < 0) { //entering into the triangle
	    if(parity > 0) {  //on outside
	      parity = -parity;
	      lastt = t;
	    }
	  }
	  else if(n[shootDirection] > 0) {  //exiting the triangle
	    if(parity < 0) {  //on inside
	      parity = -parity;
	      occupation += t-lastt;
	    }
	    else if(numHits == 0) {  //first hit exited a triangle, assume all previous cells were occupied
	      for(int i2=0;i2<i;i2++)
		density(i2,j,k) += xscale;
	      //occupation += t-lastt;
	      parity = -parity;
	    }
	  }
	  numHits++;
	}
	if(parity < 0) //on the inside
	  occupation += enddistance - lastt;
	
	//set the occupation grid
	density(i,j,k) += occupation*xscaleInv;
	if(density(i,j,k) > 1.0) density(i,j,k) = 1.0;
      }
    }
  }
}


void DensityEstimate_RandomShooting(const TriMesh& m,Array3D<Real>& density,AABB3D& bb,int numSamples,int shootDirection)
{
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(density,bb,m);

  if(shootDirection != 0) FatalError("DensityEstimate_RandomShooting: can only do x direction now");

  Ray3D ray;
  ray.direction.setZero();
  ray.direction[shootDirection] = 1;

  Array3D<list<int> > triangleLists(density.m,density.n,density.p);
  GetTriangleBuckets(m,bb,triangleLists);
  
  density.set(Zero);
  Real xscale = (bb.bmax.x-bb.bmin.x)/Real(density.m);
  Triangle3D tri;
  for(int j=0;j<density.n;j++) {
    for(int k=0;k<density.p;k++) {
      AABB3D cell;
      GetGridCell(density,bb,IntTriple(0,j,k),cell);
      ray.source.x = bb.bmin.x;

      for(int sample=0;sample<numSamples;sample++) {
	ray.source.y=Rand(cell.bmin.y,cell.bmax.y);
	ray.source.z=Rand(cell.bmin.z,cell.bmax.z);

	//parity is 1 if it's outside, -1 if inside
	int parity = 1;
	int numHits = 0;
	//sweep a ray across x axis
	for(int i=0;i<density.m;i++) {
	  Real startdistance = (Real(i))*xscale;
	  Real enddistance = (Real(i+1))*xscale;
	  const list<int>& triangles = triangleLists(i,j,k);
	  Heap<Triangle3D,Real> tris;
	  for(list<int>::const_iterator a=triangles.begin();a!=triangles.end();a++) {
	    m.GetTriangle(*a,tri);
	    Real t,u,v;
	    if(tri.rayIntersects(ray,&t,&u,&v)) {
	      if(t >= startdistance && t < enddistance)
		tris.push(tri,-t);
	    }
	  }
	  //integrate the length of segment intersected by the triangle interiors
	  Real occupation=0;
	  Real lastt = startdistance; //last intersection into the object
	  while(!tris.empty()) {
	    tri = tris.top();
	    Real t=-tris.topPriority();
	    tris.pop();

	    Vector3 n=tri.normal();
	    if(n[shootDirection] < 0) { //entering into the triangle
	      if(parity > 0) {  //on outside
		parity = -parity;
		lastt = t;
	      }
	    }
	    else if(n[shootDirection] > 0) {  //exiting the triangle
	      if(parity < 0) {  //on inside
		parity = -parity;
		occupation += t-lastt;
	      }
	      else if(numHits == 0) {  //first hit exited a triangle, assume all previous cells were occupied
		for(int i2=0;i2<i;i2++)
		  density(i2,j,k) += xscale;
		//occupation += t-lastt;
		parity = -parity;
	      }
	    }
	    numHits++;
	  }
	  if(parity < 0) //on the inside
	    occupation += enddistance - lastt;

	  //set the occupation grid
	  density(i,j,k) += occupation;
	}
      }
    }
  }

  Real scale = Real(density.m)/(numSamples*(bb.bmax.x-bb.bmin.x));
  for(int i=0;i<density.m;i++) 
    for(int j=0;j<density.n;j++) 
      for(int k=0;k<density.p;k++) {
	density(i,j,k) *= scale;
	if(density(i,j,k) > 1) density(i,j,k) = 1;
      }
}

void SweepVisibilityGrid(const TriMesh& m,int direction,Array3D<bool>& visible,AABB3D& bb,bool singleSided)
{
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(visible,bb,m);

  SurfaceOccupancyGrid(m,visible,bb);
  if(singleSided) FatalError("Can only do double-sided visibility grid as of now");
  int dim = abs(direction)-1;
  int perp1 = (dim+1)%3, perp2 = (dim+2)%3;
  int dims[3]={visible.m,visible.n,visible.p};
  int index[3];
  for(index[perp1]=0;index[perp1]<dims[perp1];index[perp1]++) {
    for(index[perp2]=0;index[perp2]<dims[perp2];index[perp2]++) {
      bool hitsurface=false;
      if(direction < 0) {
        for(index[dim]=dims[dim]-1;index[dim]>=0;index[dim]--) {
	  if(hitsurface) {
	    visible(index[0],index[1],index[2])=false;
	  }
	  else {
	    hitsurface=visible(index[0],index[1],index[2]);
	    visible(index[0],index[1],index[2])=true;
	  }
	}
      }
      else {
        for(index[dim]=0;index[dim]<dims[dim];index[dim]++) {
	  if(hitsurface) {
	    visible(index[0],index[1],index[2])=false;
	  }
	  else {
	    hitsurface=visible(index[0],index[1],index[2]);
	    visible(index[0],index[1],index[2])=true;
	  }
	}
      }
    }
  }
}

enum FMMStatus { Accepted, Considered, Far };

struct TriangleClosestPointData
{
  //given the closest point is the given vertex, sets the dir member and
  //signedDistance
  void SetVertexDistance(const TriMeshWithTopology& mesh,int vertex,
			 const Vector3& queryPoint,const Vector3& closestPoint)
  {
    dir = queryPoint-closestPoint;
    Real distance = dir.norm();
    //check orientation at edge
    if(FuzzyZero(distance)) {  //right at the apex
      cout<<"Closest point is right at the apex of vertex "<<vertex<<endl;
      //average the adjacent normals
      Vector3 aveNormal(Zero);
      for(size_t i=0;i<mesh.incidentTris[vertex].size();i++) {
	int t=mesh.incidentTris[vertex][i];
	aveNormal += mesh.TriangleNormal(t);
      }
      Real nnorm = aveNormal.norm();
      if(FuzzyZero(nnorm)) { //norm is zero?!
	dir.setZero();
	signedDistance = distance;
	cout<<"Uhh... average normal is zero??"<<endl;
      }
      else {
	aveNormal /= nnorm;
	Real sign = Sign(dot(dir,aveNormal));
	dir = aveNormal*sign;
	signedDistance = sign*distance;
      }
      getchar();
    }
    else {
      //find the triangle that's the most perpendicular to the point
      Real maxAngle=-1;
      int maxAngleIndex=-1;
      for(size_t i=0;i<mesh.incidentTris[vertex].size();i++) {
	int t=mesh.incidentTris[vertex][i];
	Real angle = Abs(dot(mesh.TriangleNormal(t),dir));
	if(angle > maxAngle) {
	  maxAngle = angle;
	  maxAngleIndex = t;
	}
      }
      Assert(maxAngleIndex >= 0);
      Real sign = Sign(dot(mesh.TriangleNormal(maxAngleIndex),dir));
      dir /= distance;
      dir *= sign;
      signedDistance = sign*distance;
    }
  }

  //given the closest point is triangleIndex on the given edge, sets
  //the dir member and signedDistance
  void SetEdgeDistance(const TriMeshWithTopology& mesh,int edge,
		       const Vector3& queryPoint,const Vector3& closestPoint)
  {
    int adjIndex = mesh.triNeighbors[triangleIndex][edge];
    dir = queryPoint-closestPoint;
    Real distance = dir.norm();
    Vector3 normal = mesh.TriangleNormal(triangleIndex);

    if(adjIndex < 0) {
      printf("No adjacent triangle on edge\n");
      //query sign by triangle normal
      Real sign = Sign(dot(normal,dir));
      signedDistance = sign*distance;
      if(FuzzyZero(distance,Zero)) {
	dir = normal;
      }
      else {
	dir /= distance;
	dir *= sign;
      }
      return;
    }

    Vector3 adjNormal = mesh.TriangleNormal(adjIndex);
    
    //check orientation at edge
    if(FuzzyZero(distance)) {  //right at the apex
      printf("Cell center right on edge\n");
      //average the adjacent normals
      Vector3 aveNormal = normal + adjNormal;
      Real nnorm = aveNormal.norm();
      if(FuzzyZero(nnorm)) { //norm is zero?!
	cerr<<"Warning: closest point is right at the edge, and the normal of an edge is zero!"<<endl;
	getchar();
	dir.setZero();
	signedDistance = 0;
      }
      else {
	aveNormal /= nnorm;
	Real sign = Sign(dot(dir,aveNormal));
	dir = aveNormal;
	signedDistance = sign*distance;
      }
    }
    else {
      //find the triangle that's the most perpendicular to the point
      Real sign;
      if(Abs(dot(normal,dir)) > Abs(dot(adjNormal,dir))) 
	sign = Sign(dot(normal,dir));
      else
	sign = Sign(dot(adjNormal,dir));

      if(sign==0.0) sign=1.0;
      dir /= distance;
      dir *= sign;
      signedDistance = sign*distance;
    }
  }


  void Update(const TriMeshWithTopology& mesh,int index,const Vector3& queryPoint) {
    Vector3 point;
    //triangleIndex = ClosestPointDescent(mesh,queryPoint,index,point);
    triangleIndex = index;
    Calculate(mesh,triangleIndex,queryPoint);
  }

  void Calculate(const TriMeshWithTopology& mesh,int index,const Vector3& queryPoint) {
    Assert(!mesh.triNeighbors.empty());
    Assert(!mesh.incidentTris.empty());
    triangleIndex = index;

    Triangle3D tri;
    mesh.GetTriangle(triangleIndex,tri);
    Vector2 planeCoords = tri.closestPointCoords(queryPoint);
    Assert(planeCoords.x >= 0 && planeCoords.y >= 0);
    Assert(planeCoords.x+planeCoords.y <= One + Epsilon);
    Vector3 point = tri.planeCoordsToPoint(planeCoords);
    dir = queryPoint-point;
    Real distance = dir.norm();


    //figure out the sign of distance from the feature
    if(FuzzyZero(planeCoords.x)) {
      if(FuzzyZero(planeCoords.y)) {
	//closest at point a
	int vertex=mesh.tris[triangleIndex][0];
	SetVertexDistance(mesh,vertex,queryPoint,point);
      }
      else if(FuzzyEquals(planeCoords.y,1.0)) {
	//closest at point c
	int vertex=mesh.tris[triangleIndex][2];
	SetVertexDistance(mesh,vertex,queryPoint,point);
      }
      else {
	//closest at edge ac
	SetEdgeDistance(mesh,2,queryPoint,point);
      }
    }
    else if(FuzzyZero(planeCoords.y)) {
      if(FuzzyEquals(planeCoords.x,1.0)) {
	//closest at point b
	int vertex=mesh.tris[triangleIndex][1];
	SetVertexDistance(mesh,vertex,queryPoint,point);
      }
      else {
	//closest at edge ab
	SetEdgeDistance(mesh,0,queryPoint,point);
      }
    }
    else {
      if(FuzzyEquals(planeCoords.x+planeCoords.y,One)) {
	//closest at edge bc
	SetEdgeDistance(mesh,1,queryPoint,point);
      }
      else {
	//closest is triangle face

	//query sign by triangle normal
	Real sign = Sign(dot(tri.normal(),dir));
	if(sign == 0) sign = 1.0;
	signedDistance = sign*distance;
	if(FuzzyZero(distance,Zero)) {
	  dir = tri.normal();
	}
	else {
	  dir /= distance;
	  dir *= sign;
	}
      }
    }
  }

  int triangleIndex;     //indexed triangle
  Real signedDistance;   //signed distance to point
  Vector3 dir;           //signed distance gradient
};

void FastMarchingMethod(const TriMeshWithTopology& m,Array3D<Real>& distance,Array3D<Vector3>& gradient,AABB3D& bb,vector<IntTriple>& surfaceCells)
{
  int M=distance.m,N=distance.n,P=distance.p;
  if(gradient.m != M || gradient.n != N || gradient.p != P) gradient.resize(M,N,P);
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(distance,bb,m);

  Array3D<int> closestTriangle(M,N,P,-1);
  Array3D<FMMStatus> status(M,N,P,Far);
  //encode an index i,j,k as (i*N+j)*P+k
  FixedSizeHeap<Real> queue(M*N*P);

  IntTriple index,lo,hi;
  //fill in initial surface distances, get surface set
  cout<<"Filling initial surfaces"<<endl;
  distance.set(Inf);
  set<IntTriple> surfaceSet;
  Triangle3D tri;
  AABB3D query,cell;
  for(size_t i=0;i<m.tris.size();i++) {
    m.GetTriangle(i,tri);
    query.setPoint(tri.a);
    query.expand(tri.b);
    query.expand(tri.c);
    bool q=QueryGrid(distance,bb,query,lo,hi);
    if(!q) continue;
    /*
    index=lo;
    while(index.a <= hi.a) {
      while(index.b <= hi.b) {
	while(index.c <= hi.c) {
	  GetGridCell(distance,bb,index,cell);
    */
    VolumeGridIterator<Real> it(distance,bb);
    it.setRange(lo,hi);
    for(;!it.isDone();++it) {
      it.getCell(cell);
      index = it.getIndex();

	  if(tri.intersects(cell)) {
	    surfaceSet.insert(index);

	    Vector3 cellCenter = (cell.bmin+cell.bmax)*Half;
	    TriangleClosestPointData cp;
	    //TODO: use triangle neighborhood information
	    cp.Calculate(m,i,cellCenter);

	    if(Abs(cp.signedDistance) < Abs(distance(index))) {
	      distance(index) = cp.signedDistance;
	      gradient(index) = cp.dir;
	      closestTriangle(index) = (int)i;
	      int heapIndex = (index.a*N+index.b)*P+index.c;
	      queue.adjust(heapIndex,-Abs(cp.signedDistance));
	    }
	  }
	  /*
	  index.c++;
	}
	index.b++;
	index.c = lo.c;
      }
      index.a++;
      index.b = lo.b;
    }
	  */
    }
  }

  //cout<<"Converting surface cells"<<endl;
  //convert surface set to vector
  surfaceCells.resize(surfaceSet.size());
  int i=0;
  for(set<IntTriple>::const_iterator it=surfaceSet.begin();it!=surfaceSet.end();it++,i++) {
    surfaceCells[i]=*it;
  }

  cout<<"FMM iterations"<<endl;
  Vector3 cellSize = (bb.bmax-bb.bmin);
  cellSize.x /= M;
  cellSize.y /= N;
  cellSize.z /= P;
  //Timer timer,overallTimer;
  /*
  int numIters=0,numCPCalls=0,numHeapUpdates=0,maxHeapSize=0;
  double heapPopTime=0,heapUpdateTime=0,closestPointTime=0,overheadTime=0;
  */
  //process the cells in the queue
  vector<IntTriple> adj;
  adj.reserve(6);
  TriangleClosestPointData cp;
  while(!queue.empty()) {
    //maxHeapSize = Max(maxHeapSize,queue.size());
    //numIters++;
    //timer.Reset();
    int heapIndex = queue.top();
    queue.pop();
    //heapPopTime += timer.ElapsedTime();

    //timer.Reset();
    //decode heap index
    Assert(heapIndex >= 0);
    index.c = heapIndex % P;
    index.b = (heapIndex / P)%N;
    index.a = (heapIndex / (P*N));
    Assert(index.a < M);

    status(index.a,index.b,index.c)=Accepted;
    int lastClosestTri=closestTriangle(index.a,index.b,index.c);
    Assert(lastClosestTri >= 0 && lastClosestTri <(int)m.tris.size());

    //add all potentially adjacent neighbors
    adj.resize(0);
    if(index.a+1 < M) { adj.push_back(index); adj.back().a++; }
    if(index.a-1 >= 0) { adj.push_back(index); adj.back().a--; }
    if(index.b+1 < N) { adj.push_back(index); adj.back().b++; }
    if(index.b-1 >= 0) { adj.push_back(index); adj.back().b--; }
    if(index.c+1 < P) { adj.push_back(index); adj.back().c++; }
    if(index.c-1 >= 0) { adj.push_back(index); adj.back().c--; }
    //overheadTime += timer.ElapsedTime();
    for(size_t j=0;j<adj.size();j++) {
      //timer.Reset();
      index = adj[j];
      if(status(index.a,index.b,index.c)==Accepted) continue;
      status(index.a,index.b,index.c)=Considered;
      Vector3 cellCenter;
      GetGridCellCenter(distance,bb,index,cellCenter);
      //overheadTime += timer.ElapsedTime();

      //numCPCalls++;
      //timer.Reset();
      cp.signedDistance = Inf;
      cp.Update(m,lastClosestTri,cellCenter);
      //closestPointTime = timer.ElapsedTime();

      if(Abs(cp.signedDistance) < Abs(distance(index.a,index.b,index.c))) {
	//timer.Reset();
	distance(index.a,index.b,index.c) = cp.signedDistance;
	gradient(index.a,index.b,index.c) = cp.dir;
	closestTriangle(index.a,index.b,index.c) = cp.triangleIndex;
	//overheadTime += timer.ElapsedTime();
	//encode heap index
	//numHeapUpdates++;
	//timer.Reset();
	int heapIndex = (index.a*N+index.b)*P+index.c;
	queue.adjust(heapIndex,-Abs(cp.signedDistance));
	//heapUpdateTime += timer.ElapsedTime();
      }
    }
  }
  for(int i=0;i<status.m;i++)
    for(int j=0;j<status.n;j++)
      for(int k=0;k<status.p;k++)
	Assert(status(i,j,k)==Accepted);
  //printf("Overall time: %g\n",overallTimer.ElapsedTime());
  //printf("Overall: %g, heap pop: %g, heap update: %g, cp: %g, overhead %g\n",overallTimer.ElapsedTime(),heapPopTime,heapUpdateTime,closestPointTime,overheadTime);
  //printf("%d iterations, %d cp calls, %d heap updates, max heap size %d\n",numIters,numCPCalls,numHeapUpdates,maxHeapSize);
}

void FastMarchingMethod_Fill(const TriMeshWithTopology& m,Array3D<Real>& distance,Array3D<Vector3>& gradient,AABB3D& bb,vector<IntTriple>& surfaceCells)
{
  int M=distance.m,N=distance.n,P=distance.p;
  if(gradient.m != M || gradient.n != N || gradient.p != P) gradient.resize(M,N,P);
  if(bb.bmin.x > bb.bmax.x || bb.bmin.y > bb.bmax.y || bb.bmin.z > bb.bmax.z)
    FitGridToMesh(distance,bb,m);

  Array3D<bool> occupied(M,N,P),interior(M,N,P),exterior(M,N,P);
  VolumeOccupancyGrid_FloodFill(m,occupied,bb,IntTriple(0,0,0),0);
  interior = occupied;
  exterior = occupied;
  for(int i=0;i<M;i++)
    for(int j=0;j<N;j++)
      for(int k=0;k<P;k++) {
	if(!occupied(i,j,k)) {
	  //erode the interior mesh
	  if(i>0) interior(i-1,j,k)=false;
	  if(i+1<M) interior(i+1,j,k)=false;
	  if(j>0) interior(i,j-1,k)=false;
	  if(j+1<N) interior(i,j+1,k)=false;
	  if(k>0) interior(i,j,k-1)=false;
	  if(k+1<P) interior(i,j,k+1)=false;
	}
	else {
	  //grow the exterior mesh
	  if(i>0) exterior(i-1,j,k)=true;
	  if(i+1<M) exterior(i+1,j,k)=true;
	  if(j>0) exterior(i,j-1,k)=true;
	  if(j+1<N) exterior(i,j+1,k)=true;
	  if(k>0) exterior(i,j,k-1)=true;
	  if(k+1<P) exterior(i,j,k+1)=true;
	}
      } 

  Array3D<int> closestTriangle(M,N,P,-1);
  Array3D<FMMStatus> status(M,N,P,Far);
  //encode an index i,j,k as (i*N+j)*P+k
  FixedSizeHeap<Real> queue(M*N*P);

  IntTriple index,lo,hi;
  //fill in initial surface distances, get surface set
  cout<<"Filling initial surfaces"<<endl;
  distance.set(Inf);
  set<IntTriple> surfaceSet;
  Triangle3D tri;
  AABB3D query,cell;
  for(size_t i=0;i<m.tris.size();i++) {
    m.GetTriangle(i,tri);
    query.setPoint(tri.a);
    query.expand(tri.b);
    query.expand(tri.c);
    bool q=QueryGrid(distance,bb,query,lo,hi);
    if(!q) continue;
    /*
    index=lo;
    while(index.a <= hi.a) {
      while(index.b <= hi.b) {
	while(index.c <= hi.c) {
	  GetGridCell(distance,bb,index,cell);
    */
    VolumeGridIterator<Real> it(distance,bb);
    it.setRange(lo,hi);
    for(;!it.isDone();++it) {
      it.getCell(cell);
      index = it.getIndex();

      if(!interior(index.a,index.b,index.c) && tri.intersects(cell)) {
	    surfaceSet.insert(index);

	    Vector3 cellCenter = (cell.bmin+cell.bmax)*Half;
	    TriangleClosestPointData cp;
	    //TODO: use triangle neighborhood information
	    cp.Calculate(m,i,cellCenter);

	    if(Abs(cp.signedDistance) < Abs(distance(index))) {
	      distance(index) = cp.signedDistance;
	      gradient(index) = cp.dir;
	      closestTriangle(index) = (int)i;
	      int heapIndex = (index.a*N+index.b)*P+index.c;
	      queue.adjust(heapIndex,-Abs(cp.signedDistance));
	      status(index.a,index.b,index.c)=Accepted;
	    }
	  }
	  /*
	  index.c++;
	}
	index.b++;
	index.c = lo.c;
      }
      index.a++;
      index.b = lo.b;
    }
	  */
    }
  }

  //cout<<"Converting surface cells"<<endl;
  //convert surface set to vector
  surfaceCells.resize(surfaceSet.size());
  int i=0;
  for(set<IntTriple>::const_iterator it=surfaceSet.begin();it!=surfaceSet.end();it++,i++) {
    surfaceCells[i]=*it;
  }

  cout<<"FMM iterations"<<endl;
  Vector3 cellSize = (bb.bmax-bb.bmin);
  cellSize.x /= M;
  cellSize.y /= N;
  cellSize.z /= P;
  //Timer timer,overallTimer;
  /*
  int numIters=0,numCPCalls=0,numHeapUpdates=0,maxHeapSize=0;
  double heapPopTime=0,heapUpdateTime=0,closestPointTime=0,overheadTime=0;
  */
  //process the cells in the queue
  vector<IntTriple> adj;
  adj.reserve(6);
  TriangleClosestPointData cp;
  while(!queue.empty()) {
    //maxHeapSize = Max(maxHeapSize,queue.size());
    //numIters++;
    //timer.Reset();
    int heapIndex = queue.top();
    queue.pop();
    //heapPopTime += timer.ElapsedTime();

    //timer.Reset();
    //decode heap index
    Assert(heapIndex >= 0);
    index.c = heapIndex % P;
    index.b = (heapIndex / P)%N;
    index.a = (heapIndex / (P*N));
    Assert(index.a < M);

    status(index.a,index.b,index.c)=Accepted;
    int lastClosestTri=closestTriangle(index.a,index.b,index.c);
    Assert(lastClosestTri >= 0 && lastClosestTri <(int)m.tris.size());

    //add all potentially adjacent neighbors
    adj.resize(0);
    if(index.a+1 < M) { adj.push_back(index); adj.back().a++; }
    if(index.a-1 >= 0) { adj.push_back(index); adj.back().a--; }
    if(index.b+1 < N) { adj.push_back(index); adj.back().b++; }
    if(index.b-1 >= 0) { adj.push_back(index); adj.back().b--; }
    if(index.c+1 < P) { adj.push_back(index); adj.back().c++; }
    if(index.c-1 >= 0) { adj.push_back(index); adj.back().c--; }
    //overheadTime += timer.ElapsedTime();
    for(size_t j=0;j<adj.size();j++) {
      //timer.Reset();
      index = adj[j];
      if(status(index.a,index.b,index.c)==Accepted) continue;
      status(index.a,index.b,index.c)=Considered;
      Vector3 cellCenter;
      GetGridCellCenter(distance,bb,index,cellCenter);
      //overheadTime += timer.ElapsedTime();

      //numCPCalls++;
      //timer.Reset();
      cp.signedDistance = Inf;
      cp.Update(m,lastClosestTri,cellCenter);
      //closestPointTime = timer.ElapsedTime();

      if(cp.signedDistance < 0 && !exterior(index)) {
	//printf("Warning, exterior exploration led to negative signed distance?\n");
	continue;
      }
      else if(cp.signedDistance > 0 && interior(index)) {
	//printf("Warning, interior exploration led to positive signed distance?\n");
	continue;
      }
      if(Abs(cp.signedDistance) < Abs(distance(index.a,index.b,index.c))) {
	//timer.Reset();
	distance(index.a,index.b,index.c) = cp.signedDistance;
	gradient(index.a,index.b,index.c) = cp.dir;
	closestTriangle(index.a,index.b,index.c) = cp.triangleIndex;
	//overheadTime += timer.ElapsedTime();
	//encode heap index
	//numHeapUpdates++;
	//timer.Reset();
	int heapIndex = (index.a*N+index.b)*P+index.c;
	queue.adjust(heapIndex,-Abs(cp.signedDistance));
	//heapUpdateTime += timer.ElapsedTime();
      }
    }
  }

  for(int i=0;i<status.m;i++)
    for(int j=0;j<status.n;j++)
      for(int k=0;k<status.p;k++) {
	if(status(i,j,k)!=Accepted) {
	  //fill in with neighbors
	  if(interior(i,j,k)) {
	    Real minNeighbor = Inf;
	    if(i > 0 && status(i-1,j,k)==Accepted) minNeighbor=Min(minNeighbor,distance(i-1,j,k));
	    if(i+1 < M && status(i+1,j,k)==Accepted) minNeighbor=Min(minNeighbor,distance(i+1,j,k));
	    if(j > 0 && status(i,j-1,k)==Accepted) minNeighbor=Min(minNeighbor,distance(i,j-1,k));
	    if(j+1 < N && status(i,j+1,k)==Accepted) minNeighbor=Min(minNeighbor,distance(i,j+1,k));
	    if(k > 0 && status(i,j,k-1)==Accepted) minNeighbor=Min(minNeighbor,distance(i,j,k-1));
	    if(k+1 < P && status(i,j,k+1)==Accepted) minNeighbor=Min(minNeighbor,distance(i,j,k+1));
	    distance(i,j,k) = minNeighbor; 
	  }
	  else {
	    Real maxNeighbor = Inf;
	    if(i > 0 && status(i-1,j,k)==Accepted) maxNeighbor=Max(maxNeighbor,distance(i-1,j,k));
	    if(i+1 < M && status(i+1,j,k)==Accepted) maxNeighbor=Max(maxNeighbor,distance(i+1,j,k));
	    if(j > 0 && status(i,j-1,k)==Accepted) maxNeighbor=Max(maxNeighbor,distance(i,j-1,k));
	    if(j+1 < N && status(i,j+1,k)==Accepted) maxNeighbor=Max(maxNeighbor,distance(i,j+1,k));
	    if(k > 0 && status(i,j,k-1)==Accepted) maxNeighbor=Max(maxNeighbor,distance(i,j,k-1));
	    if(k+1 < P && status(i,j,k+1)==Accepted) maxNeighbor=Max(maxNeighbor,distance(i,j,k+1));
	    distance(i,j,k) = maxNeighbor; 
	  }
	}
      }
  //printf("Overall time: %g\n",overallTimer.ElapsedTime());
  //printf("Overall: %g, heap pop: %g, heap update: %g, cp: %g, overhead %g\n",overallTimer.ElapsedTime(),heapPopTime,heapUpdateTime,closestPointTime,overheadTime);
  //printf("%d iterations, %d cp calls, %d heap updates, max heap size %d\n",numIters,numCPCalls,numHeapUpdates,maxHeapSize);
}


void DensityEstimate_FMM(const TriMeshWithTopology& m,Array3D<Real>& density,AABB3D& bb)
{
  Array3D<Real> distance(density.m,density.n,density.p);
  Array3D<Vector3> gradient;
  vector<IntTriple> surfaceCells;
  FastMarchingMethod(m,distance,gradient,bb,surfaceCells);

  //get some consistency here
  Array3D<bool> occupied(density.m,density.n,density.p);
  VolumeOccupancyGrid_FloodFill(m,occupied,bb,IntTriple(0,0,0),false);
  Vector3 resolution = (bb.bmax-bb.bmin)*0.5;
  resolution.x /= density.m;
  resolution.y /= density.n;
  resolution.z /= density.p;
  for(int i=0;i<density.m;i++) 
    for(int j=0;j<density.n;j++) 
      for(int k=0;k<density.p;k++) {
	//can't rule on points close to the surface
	Vector3 surfpt = distance(i,j,k)*gradient(i,j,k);
	if(Abs(surfpt.x) < resolution.x && Abs(surfpt.y) < resolution.y && Abs(surfpt.z) < resolution.z) continue;
	if((distance(i,j,k) > 0) == occupied(i,j,k))
	  distance(i,j,k) *= -1.0;
      }

  DensityEstimate_FMM(distance,gradient,bb,density);
}

void DensityEstimate_FMM(const Array3D<Real>& distance,const Array3D<Vector3>& gradient,const AABB3D& bb,Array3D<Real>& density)
{
  Assert(distance.m == density.m && distance.n == density.n && distance.p == density.p);
  Assert(gradient.m == density.m && gradient.n == density.n && gradient.p == density.p);

  AABB3D cell;
  Plane3D p;
  for(int i=0;i<density.m;i++) {
    for(int j=0;j<density.n;j++) {
      for(int k=0;k<density.p;k++) {
	GetGridCell(density,bb,IntTriple(i,j,k),cell);
	p.normal = gradient(i,j,k);
	p.offset = dot(p.normal,(cell.bmin+cell.bmax)*0.5)-distance(i,j,k);
	density(i,j,k) = GridCellDensity(cell,p);
      }
    }
  }
}

} //namespace Meshing
