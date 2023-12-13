#include "precomp.h"
#include "assignment_1.h"
#include "stb_image_write.h"
#include "OBJ_Loader.h"
#include <iostream>
#include <string>

// THIS SOURCE FILE:
// Code for the article "How to Build a BVH", part 2: faster rays.
// This version improves ray traversal speed using ordered traversal
// and better split plane orientation / positions based on the SAH
// (Surface Area Heuristic).
// Feel free to copy this code to your own framework. Absolutely no
// rights are reserved. No responsibility is accepted either.
// For updates, follow me on twitter: @j_bikker.

TheApp* CreateApp() { return new BihApp(); }

// enable the use of SSE in the AABB intersection function
#define USE_SSE

//#define SCENE 0 //ROBOT
#define SCENE 1 //LANDSCAPE

#define METHOD 0 // BVH
//#define METHOD 1 // Grid
//#define METHOD 2 // BIH

#if METHOD == 0
	#define METHOD_NAME "BVH"
#elif METHOD == 1
	#define METHOD_NAME "Grid"
#elif METHOD == 2
	#define METHOD_NAME "BIH"
#endif

#if SCENE == 0
	#define N 12582
	#define OBJ_NAME "ROBOT"
#elif SCENE == 1
	#define N 98304/3
	#define OBJ_NAME "LANDSCAPE"
	#define OBJ_PATH "procedural ground.obj"
#endif

// minimal structs
struct Tri { float3 vertex0, vertex1, vertex2; float3 centroid; };
//TODO: Optimize struct with union
struct aabb
{
	float3 bmin = 1e30f, bmax = -1e30f;
	void grow( float3 p ) { bmin = fminf( bmin, p ); bmax = fmaxf( bmax, p ); }
	float area()
	{
		float3 e = bmax - bmin; // box extent
		return e.x * e.y + e.y * e.z + e.z * e.x;
	}
};
__declspec(align(64)) struct Ray
{
	Ray() { O4 = D4 = rD4 = _mm_set1_ps( 1 ); }
	union { struct { float3 O; float dummy1; }; __m128 O4; };
	union { struct { float3 D; float dummy2; }; __m128 D4; };
	union { struct { float3 rD; float dummy3; }; __m128 rD4; };
	float t = 1e30f;
};

// Method-specific structs
#if METHOD == 0
struct BVHNode
{
	union { struct { float3 aabbMin; uint leftFirst; }; __m128 aabbMin4; };
	union { struct { float3 aabbMax; uint triCount; }; __m128 aabbMax4; };
	bool isLeaf() { return triCount > 0; }
};
#elif METHOD == 1
#define FLOAT_MAX  3.402823466e+38
#define FLOAT_MIN  1.175494351e-38
#if SCENE == 0
#define GRID_SIZE  64
#elif SCENE == 1
#define GRID_SIZE  128
#endif

struct GridCell {
public:
	vector<int> triangles;
};
struct Grid {
public:
	GridCell grid[GRID_SIZE][GRID_SIZE][GRID_SIZE];
	float3 min = float3(FLOAT_MAX); //min coordinates on the grid
	float3 max = float3(FLOAT_MIN); //max coordinates on the grid
	float3 origin = float3(0);
	float3 size = float3(0.0f);
	float3 cellSize = float3(0.0f);
};

struct IntersectionData {
public:
	int triIdx; //triangle index of intersection with ray
	float t; //Ray's t parameter value at time of intersection
};
#elif METHOD == 2
struct BihNode {
	int type; // Lowest bits: axis (00, 01, 10) or leaf(11)
	int leftIdx; // Node only
	union {
		float clip[2];
		int triangle;
	};
};
#endif

// application data
Tri tri[N];

#if METHOD == 0
uint triIdx[N];
BVHNode* bvhNode = 0;
uint rootNodeIdx = 0, nodesUsed = 2;
#elif METHOD == 1
Grid grid;
#elif METHOD == 2
uint triIdx[N];
BihNode bihNode[N * 2];
uint bihRootIdx = 0, nodesUsed = 2;
aabb totalBoundary;
#endif

// Analysis data
int traversal_steps[SCRHEIGHT * SCRWIDTH];
int triangle_intersects[SCRHEIGHT * SCRWIDTH];
uint8_t intrs_data[SCRWIDTH * SCRHEIGHT * 3];
int tri_intrs_count;
int min_traversal_steps;
int min_triangle_intersections;
int max_traversal_steps;
int max_triangle_intersects;

// functions
void IntersectTri( Ray& ray, const Tri& tri )
{
	tri_intrs_count++;
	const float3 edge1 = tri.vertex1 - tri.vertex0;
	const float3 edge2 = tri.vertex2 - tri.vertex0;
	const float3 h = cross( ray.D, edge2 );
	const float a = dot( edge1, h );
	if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
	const float f = 1 / a;
	const float3 s = ray.O - tri.vertex0;
	const float u = f * dot( s, h );
	if (u < 0 || u > 1) return;
	const float3 q = cross( s, edge1 );
	const float v = f * dot( ray.D, q );
	if (v < 0 || u + v > 1) return;
	const float t = f * dot( edge2, q );
	if (t > 0.0001f) ray.t = min( ray.t, t );
}

inline float IntersectAABB( const Ray& ray, const float3 bmin, const float3 bmax )
{
	float tx1 = (bmin.x - ray.O.x) * ray.rD.x, tx2 = (bmax.x - ray.O.x) * ray.rD.x;
	float tmin = min( tx1, tx2 ), tmax = max( tx1, tx2 );
	float ty1 = (bmin.y - ray.O.y) * ray.rD.y, ty2 = (bmax.y - ray.O.y) * ray.rD.y;
	tmin = max( tmin, min( ty1, ty2 ) ), tmax = min( tmax, max( ty1, ty2 ) );
	float tz1 = (bmin.z - ray.O.z) * ray.rD.z, tz2 = (bmax.z - ray.O.z) * ray.rD.z;
	tmin = max( tmin, min( tz1, tz2 ) ), tmax = min( tmax, max( tz1, tz2 ) );
	if (tmax >= tmin && tmin < ray.t && tmax > 0) return tmin; else return 1e30f;
}

float IntersectAABB_SSE( const Ray& ray, const __m128& bmin4, const __m128& bmax4 )
{
	static __m128 mask4 = _mm_cmpeq_ps( _mm_setzero_ps(), _mm_set_ps( 1, 0, 0, 0 ) );
	__m128 t1 = _mm_mul_ps( _mm_sub_ps( _mm_and_ps( bmin4, mask4 ), ray.O4 ), ray.rD4 );
	__m128 t2 = _mm_mul_ps( _mm_sub_ps( _mm_and_ps( bmax4, mask4 ), ray.O4 ), ray.rD4 );
	__m128 vmax4 = _mm_max_ps( t1, t2 ), vmin4 = _mm_min_ps( t1, t2 );
	float tmax = min( vmax4.m128_f32[0], min( vmax4.m128_f32[1], vmax4.m128_f32[2] ) );
	float tmin = max( vmin4.m128_f32[0], max( vmin4.m128_f32[1], vmin4.m128_f32[2] ) );
	if (tmax >= tmin && tmin < ray.t && tmax > 0) return tmin; else return 1e30f;
}

// Method specific functions
#if METHOD == 0
// forward declarations
void Subdivide(uint nodeIdx);
void UpdateNodeBounds(uint nodeIdx);

int IntersectBVH(Ray& ray)
{
	BVHNode* node = &bvhNode[rootNodeIdx], * stack[64];
	uint stackPtr = 0;
	int counter = 0;
	while (1)
	{
		counter += 1;
		if (node->isLeaf())
		{
			for (uint i = 0; i < node->triCount; i++)
				IntersectTri(ray, tri[triIdx[node->leftFirst + i]]);
			if (stackPtr == 0) break; else node = stack[--stackPtr];
			continue;
		}
		BVHNode* child1 = &bvhNode[node->leftFirst];
		BVHNode* child2 = &bvhNode[node->leftFirst + 1];
#ifdef USE_SSE
		float dist1 = IntersectAABB_SSE(ray, child1->aabbMin4, child1->aabbMax4);
		float dist2 = IntersectAABB_SSE(ray, child2->aabbMin4, child2->aabbMax4);
#else
		float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax);
		float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax);
#endif
		if (dist1 > dist2) { swap(dist1, dist2); swap(child1, child2); }
		if (dist1 == 1e30f)
		{
			if (stackPtr == 0) break; else node = stack[--stackPtr];
		}
		else
		{
			node = child1;
			if (dist2 != 1e30f) stack[stackPtr++] = child2;
		}
	}
	return counter;
}

void BuildBVH()
{
	// create the BVH node pool
	bvhNode = (BVHNode*)_aligned_malloc(sizeof(BVHNode) * N * 2, 64);
	// populate triangle index array
	for (int i = 0; i < N; i++) triIdx[i] = i;
	// calculate triangle centroids for partitioning
	for (int i = 0; i < N; i++)
		tri[i].centroid = (tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.3333f;
	// assign all triangles to root node
	BVHNode& root = bvhNode[rootNodeIdx];
	root.leftFirst = 0, root.triCount = N;
	UpdateNodeBounds(rootNodeIdx);
	// subdivide recursively
	Timer t;
	Subdivide(rootNodeIdx);
	printf("BVH (%i nodes) constructed in %.2fms.\n", nodesUsed, t.elapsed() * 1000);
}

void UpdateNodeBounds(uint nodeIdx)
{
	BVHNode& node = bvhNode[nodeIdx];
	node.aabbMin = float3(1e30f);
	node.aabbMax = float3(-1e30f);
	for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
	{
		uint leafTriIdx = triIdx[first + i];
		Tri& leafTri = tri[leafTriIdx];
		node.aabbMin = fminf(node.aabbMin, leafTri.vertex0);
		node.aabbMin = fminf(node.aabbMin, leafTri.vertex1);
		node.aabbMin = fminf(node.aabbMin, leafTri.vertex2);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex0);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex1);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex2);
	}
}

float EvaluateSAH(BVHNode& node, int axis, float pos)
{
	// determine triangle counts and bounds for this split candidate
	aabb leftBox, rightBox;
	int leftCount = 0, rightCount = 0;
	for (uint i = 0; i < node.triCount; i++)
	{
		Tri& triangle = tri[triIdx[node.leftFirst + i]];
		if (triangle.centroid[axis] < pos)
		{
			leftCount++;
			leftBox.grow(triangle.vertex0);
			leftBox.grow(triangle.vertex1);
			leftBox.grow(triangle.vertex2);
		}
		else
		{
			rightCount++;
			rightBox.grow(triangle.vertex0);
			rightBox.grow(triangle.vertex1);
			rightBox.grow(triangle.vertex2);
		}
	}
	float cost = leftCount * leftBox.area() + rightCount * rightBox.area();
	return cost > 0 ? cost : 1e30f;
}

void Subdivide(uint nodeIdx)
{
	// terminate recursion
	BVHNode& node = bvhNode[nodeIdx];
	// determine split axis using SAH
	int bestAxis = -1;
	float bestPos = 0, bestCost = 1e30f;
	for (int axis = 0; axis < 3; axis++) for (uint i = 0; i < node.triCount; i++)
	{
		Tri& triangle = tri[triIdx[node.leftFirst + i]];
		float candidatePos = triangle.centroid[axis];
		float cost = EvaluateSAH(node, axis, candidatePos);
		if (cost < bestCost)
			bestPos = candidatePos, bestAxis = axis, bestCost = cost;
	}
	int axis = bestAxis;
	float splitPos = bestPos;
	float3 e = node.aabbMax - node.aabbMin; // extent of parent
	float parentArea = e.x * e.y + e.y * e.z + e.z * e.x;
	float parentCost = node.triCount * parentArea;
	if (bestCost >= parentCost) return;
	// in-place partition
	int i = node.leftFirst;
	int j = i + node.triCount - 1;
	while (i <= j)
	{
		if (tri[triIdx[i]].centroid[axis] < splitPos)
			i++;
		else
			swap(triIdx[i], triIdx[j--]);
	}
	// abort split if one of the sides is empty
	int leftCount = i - node.leftFirst;
	if (leftCount == 0 || leftCount == node.triCount) return;
	// create child nodes
	int leftChildIdx = nodesUsed++;
	int rightChildIdx = nodesUsed++;
	bvhNode[leftChildIdx].leftFirst = node.leftFirst;
	bvhNode[leftChildIdx].triCount = leftCount;
	bvhNode[rightChildIdx].leftFirst = i;
	bvhNode[rightChildIdx].triCount = node.triCount - leftCount;
	node.leftFirst = leftChildIdx;
	node.triCount = 0;
	UpdateNodeBounds(leftChildIdx);
	UpdateNodeBounds(rightChildIdx);
	// recurse
	Subdivide(leftChildIdx);
	Subdivide(rightChildIdx);
}
#elif METHOD == 1
bool CheckGridCell(int x, int y, int z, Ray& ray, vector<IntersectionData>& intersected_triangles) {
	//if(camera_position == 3) cout << x << "," << y << "," << z << endl;
	GridCell gc = grid.grid[x][y][z];
	vector<int> triangles = gc.triangles;
	float t = ray.t;
	for (int i = 0; i < triangles.size(); i++) {
		IntersectTri(ray, tri[triangles[i]]);
		if (ray.t < 1e30f) {
			//cout << ray.t << endl;
			t = min(t, ray.t);
			/*	IntersectionData id;
				id.t = ray.t;
				id.triIdx = triangles[i];
				intersected_triangles.push_back(id);*/
		}
	}

	if (ray.t < 1e30f) {
		ray.t = t;
		return true;
	}
	return false;
}

bool FindRayGridIntersections(Ray& ray, float3& intersectionpoint, float3& exitpoint) {
	float tNear, tFar;

	float3 neart, fart;

	// Check for parallel rays that will pass the grid
	if (ray.D.x == 0) if (ray.O.x < grid.min.x || ray.O.x > grid.max.x) return false;
	if (ray.D.y == 0) if (ray.O.y < grid.min.y || ray.O.y > grid.max.y) return false;
	if (ray.D.z == 0) if (ray.O.z < grid.min.z || ray.O.z > grid.max.z) return false; 

	// compute t at  both (enter and exit) intersection points for each axis
	neart.x = (grid.min.x - ray.O.x) / ray.D.x;
	fart.x = (grid.max.x - ray.O.x) / ray.D.x;
	neart.y = (grid.min.y - ray.O.y) / ray.D.y;
	fart.y = (grid.max.y - ray.O.y) / ray.D.y;
	neart.z = (grid.min.z - ray.O.z) / ray.D.z;
	fart.z = (grid.max.z - ray.O.z) / ray.D.z;

	// find which one is the entrance point and which is the exit point
	tNear = max({ min(neart.x, fart.x), min(neart.y, fart.y), min(neart.z, fart.z) });
	tFar = min({ max(neart.x, fart.x), max(neart.y, fart.y), max(neart.z, fart.z) });

	// Check if there is any valid intersection
	if (tNear > tFar || tFar < 0) {
		return false; 
	}

	// Calculate both intersection points
	intersectionpoint = ray.O + tNear * ray.D;
	exitpoint = ray.O + tFar * ray.D;

	//if (abs(exitpoint.x - intersectionpoint.x) > 1.0f || abs(exitpoint.y - intersectionpoint.y) > 1.0f || abs(exitpoint.z - intersectionpoint.z) > 1.0f) {
	//	std::cout << intersectionpoint.x << ", " << intersectionpoint.y << ", " << intersectionpoint.z << std::endl;
	//	std::cout << exitpoint.x << ", " << exitpoint.y << ", " << exitpoint.z << std::endl;
	//}
	return true;
}

int IntersectGrid(Ray& ray) {
	int traversalSteps = 0;
	float3 intersectionpoint, exitpoint;
	if (!FindRayGridIntersections(ray, intersectionpoint, exitpoint)) return traversalSteps;

	int dx = abs(exitpoint.x - intersectionpoint.x) / grid.cellSize.x;
	int dy = abs(exitpoint.y - intersectionpoint.y) / grid.cellSize.y;
	int dz = abs(exitpoint.z - intersectionpoint.z) / grid.cellSize.z;

	int stepX = intersectionpoint.x < exitpoint.x ? 1 : -1;
	int stepY = intersectionpoint.y < exitpoint.y ? 1 : -1;
	int stepZ = intersectionpoint.z < exitpoint.z ? 1 : -1;
	float hypotenuse = sqrt(dx * dx + dy * dy + dz * dz);
	float tMaxX = hypotenuse * (grid.cellSize.x / 2.0f) / dx;
	float tMaxY = hypotenuse * (grid.cellSize.y / 2.0f) / dy;
	float tMaxZ = hypotenuse * (grid.cellSize.z / 2.0f) / dz;
	float tDeltaX = hypotenuse / dx;
	float tDeltaY = hypotenuse / dy;
	float tDeltaZ = hypotenuse / dz;
	//float tDeltaX = dx / abs(ray.D.x);
	//float tDeltaY = dy / abs(ray.D.y);
	//float tDeltaZ = dz / abs(ray.D.z);

	//int cellIndex[] = { 0, 0, 0 };
	int x0 = floor((intersectionpoint.x - grid.min.x) / grid.cellSize.x);
	x0 = clamp(x0, 0, GRID_SIZE - 1);

	//cellIndex[0] = clamp(cellIndex[0], 0, GRID_SIZE - 1);
	int y0 = floor((intersectionpoint.y - grid.min.y) / grid.cellSize.y);
	y0 = clamp(y0, 0, GRID_SIZE - 1);

	//cellIndex[1] = clamp(cellIndex[1], 0, GRID_SIZE - 1);
	int z0 = floor((intersectionpoint.z - grid.min.z) / grid.cellSize.z);
	z0 = clamp(z0, 0, GRID_SIZE - 1);

	//cellIndex[2] = clamp(cellIndex[2], 0, GRID_SIZE - 1);

	vector<IntersectionData> intersected_triangles;
	float t;
	while (1) {
		traversalSteps++;
		tri_intrs_count++;
		t = min({ tMaxX, tMaxY, tMaxZ });
		if (CheckGridCell(x0, y0, z0, ray, intersected_triangles)) {
			if (ray.t < t) return traversalSteps;
		}
		//CheckGridCell(x0, y0, z0, ray, intersected_triangles);
		if (tMaxX < tMaxY) {
			if (tMaxX < tMaxZ) {
				x0 = x0 + stepX;
				tMaxX = tMaxX + tDeltaX;
			}
			else if (tMaxX > tMaxZ) {
				z0 = z0 + stepZ;
				tMaxZ = tMaxZ + tDeltaZ;
			}
			else {
				x0 = x0 + stepX;
				tMaxX = tMaxX + tDeltaX;
				z0 = z0 + stepZ;
				tMaxZ = tMaxZ + tDeltaZ;
			}
		}
		else if (tMaxX > tMaxY) {
			if (tMaxY < tMaxZ) {
				y0 = y0 + stepY;
				tMaxY = tMaxY + tDeltaY;
			}
			else if (tMaxY > tMaxZ) {
				z0 = z0 + stepZ;
				tMaxZ = tMaxZ + tDeltaZ;
			}
			else {
				y0 = y0 + stepY;
				tMaxY = tMaxY + tDeltaY;
				z0 = z0 + stepZ;
				tMaxZ = tMaxZ + tDeltaZ;

			}
		}
		else {
			if (tMaxY < tMaxZ) {
				y0 = y0 + stepY;
				tMaxY = tMaxY + tDeltaY;
				x0 = x0 + stepX;
				tMaxX = tMaxX + tDeltaX;
			}
			else if (tMaxY > tMaxZ) {
				z0 = z0 + stepZ;
				tMaxZ = tMaxZ + tDeltaZ;
			}
			else {
				x0 = x0 + stepX;
				tMaxX = tMaxX + tDeltaX;
				y0 = y0 + stepY;
				tMaxY = tMaxY + tDeltaY;
				z0 = z0 + stepZ;
				tMaxZ = tMaxZ + tDeltaZ;

			}
		}
		if (x0 >= GRID_SIZE || x0 < 0 || y0 >= GRID_SIZE || y0 < 0 || z0 >= GRID_SIZE || z0 < 0) return traversalSteps;

		//std::cout << x0 << ", " << y0 << ", " << z0 << std::endl;
	}
}


void BuildGrid() {
	grid.size.x = (grid.max.x - grid.min.x);
	grid.size.y = (grid.max.y - grid.min.y);
	grid.size.z = (grid.max.z - grid.min.z);
	grid.origin.x = grid.min.x + (grid.size.x / 2.0f);
	grid.origin.y = grid.min.y + (grid.size.y / 2.0f);
	grid.origin.z = grid.min.z + (grid.size.z / 2.0f);
	grid.cellSize.x = grid.size.x / GRID_SIZE;
	grid.cellSize.y = grid.size.y / GRID_SIZE;
	grid.cellSize.z = grid.size.z / GRID_SIZE;
	for (int i = 0; i < N; i++) {
		float xmin = min(tri[i].vertex0.x, min(tri[i].vertex1.x, tri[i].vertex2.x));
		float xmax = max(tri[i].vertex0.x, max(tri[i].vertex1.x, tri[i].vertex2.x));
		float ymin = min(tri[i].vertex0.y, min(tri[i].vertex1.y, tri[i].vertex2.y));
		float ymax = max(tri[i].vertex0.y, max(tri[i].vertex1.y, tri[i].vertex2.y));
		float zmin = min(tri[i].vertex0.z, min(tri[i].vertex1.z, tri[i].vertex2.z));
		float zmax = max(tri[i].vertex0.z, max(tri[i].vertex1.z, tri[i].vertex2.z));

		//	std::cout << grid.cellSize.x << ", " << grid.cellSize.y << ", " << grid.cellSize.z << std::endl;


		//std::cout << "i" << i << std::endl;
		//std::cout << xmin << ", " << xmax << ", " << std::endl;
		//std::cout << ymin << ", " << ymax << ", " << std::endl;
		//std::cout << zmin << ", " << zmax << ", " << std::endl;
		int xstart = floor(((xmin - grid.min.x) / grid.cellSize.x));
		int ystart = floor((ymin - grid.min.y) / grid.cellSize.y);
		int zstart = floor((zmin - grid.min.z) / grid.cellSize.z);
		xstart = clamp(xstart - 1, 0, GRID_SIZE - 1);
		ystart = clamp(ystart - 1, 0, GRID_SIZE - 1);
		zstart = clamp(zstart - 1, 0, GRID_SIZE - 1);
		int xend = floor(((xmax - grid.min.x) / grid.cellSize.x));
		int yend = floor((ymax - grid.min.y) / grid.cellSize.y);
		int zend = floor((zmax - grid.min.z) / grid.cellSize.z);
		xend = clamp(xend + 1, 0, GRID_SIZE - 1);
		yend = clamp(yend + 1, 0, GRID_SIZE - 1);
		zend = clamp(zend + 1, 0, GRID_SIZE - 1);
		int xsize = xend - xstart;
		int ysize = yend - ystart;
		int zsize = zend - zstart;
		xsize = clamp(xsize, 0, GRID_SIZE - 1);
		ysize = clamp(ysize, 0, GRID_SIZE - 1);
		zsize = clamp(zsize, 0, GRID_SIZE - 1);


		for (int x = xstart; x <= xstart + xsize; x++) {
			//std::cout << "x" << x << std::endl;
			for (int y = ystart; y <= ystart + ysize; y++) {
				//std::cout << "y" << y << std::endl;
				for (int z = zstart; z <= zstart + zsize; z++) {
					grid.grid[x][y][z].triangles.push_back(i);
				}
			}
}
	}
}

void CheckBounds(float3 v0, float3 v1, float3 v2) {
	float minx = min({ v0.x, v1.x, v2.x }), maxx = max({ v0.x, v1.x, v2.x });
	float miny = min({ v0.y, v1.y, v2.y }), maxy = max({ v0.y, v1.y, v2.y });
	float minz = min({ v0.z, v1.z, v2.z }), maxz = max({ v0.z, v1.z, v2.z });

	if (minx < grid.min.x) grid.min.x = minx;
	if (miny < grid.min.y) grid.min.y = miny;
	if (minz < grid.min.z) grid.min.z = minz;
	if (maxx > grid.max.x) grid.max.x = maxx;
	if (maxy > grid.max.y) grid.max.y = maxy;
	if (maxz > grid.max.z) grid.max.z = maxz;
}
#elif METHOD == 2
void Subdivide(int nodeId, aabb boundary, int count) {
	BihNode& node = bihNode[nodeId];
	if (count <= 1) {
		node.type = 3;
		return;
	}
	int longestAxis = 0;
	float longestSide = 0.0;
	for (int a = 0; a < 3; a++) {
		float length = boundary.bmax[a] - boundary.bmin[a];
		if (length > longestSide) {
			longestAxis = a;
			longestSide = length;
		}
	}
	float splitPlane = boundary.bmin[longestAxis] + (boundary.bmax[longestAxis] - boundary.bmin[longestAxis]) * 0.5;
	int i = node.triangle;
	int j = i + count - 1;
	// TODO: Fix object sorting
	while (i <= j) {
		if (tri[triIdx[i]].centroid[longestAxis] < splitPlane)
			i++;
		else
			swap(triIdx[i], triIdx[j--]);
	}
	int leftCount = i - node.triangle;
	if (leftCount == 0) {
		// All nodes are in right child
		boundary.bmin[longestAxis] = splitPlane;
		Subdivide(nodeId, boundary, count);
		return;
	}
	else if (leftCount == count) {
		// All nodes are in left child
		boundary.bmax[longestAxis] = splitPlane;
		Subdivide(nodeId, boundary, count);
		return;
	}
	// TODO: Squeeze boundary (probably just do a max from boundary.bmin[a]
	float leftThresh = boundary.bmin[longestAxis];
	for (int k = node.triangle; k < i; k++) {
		Tri triangle = tri[triIdx[k]];
		leftThresh = max(
			max(
				leftThresh,
				triangle.vertex0[longestAxis]
			),
			max(
				triangle.vertex1[longestAxis],
				triangle.vertex2[longestAxis]
			)
		);
	}
	float rightThresh = boundary.bmax[longestAxis];
	for (int k = i; k < node.triangle + count; k++) {
		Tri triangle = tri[triIdx[k]];
		rightThresh = min(
			min(
				rightThresh,
				triangle.vertex0[longestAxis]
			),
			min(
				triangle.vertex1[longestAxis],
				triangle.vertex2[longestAxis]
			)
		);
	}
	// Calculate right threshold
	int leftChildIdx = nodesUsed;
	nodesUsed += 2;
	bihNode[leftChildIdx].triangle = node.triangle;
	bihNode[leftChildIdx+1].triangle = i;
	node.type = longestAxis;
	node.leftIdx = leftChildIdx;
	node.clip[0] = leftThresh;
	node.clip[1] = rightThresh;
	aabb leftBoundary = boundary;
	leftBoundary.bmax[longestAxis] = leftThresh;
	aabb rightBoundary = boundary;
	rightBoundary.bmin[longestAxis] = rightThresh;
	Subdivide(leftChildIdx, leftBoundary, leftCount);
	Subdivide(leftChildIdx+1, rightBoundary, count - leftCount);
}

void BuildBih() {
	// populate triangle index array
	for (int i = 0; i < N; i++) triIdx[i] = i;
	// calculate triangle centroids for partitioning
	for (int i = 0; i < N; i++)
		tri[i].centroid = (tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.3333f;
	// Calculate boundary for space
	for (uint i = 0; i < N; i++) {
		uint index = triIdx[i];
		Tri& triangle = tri[index];
		totalBoundary.grow(triangle.vertex0);
		totalBoundary.grow(triangle.vertex1);
		totalBoundary.grow(triangle.vertex2);
	}
	// Subdivide recursively
	Timer t;
	BihNode& root = bihNode[bihRootIdx];
	// Set root parameters
	root.triangle = 0;
	Subdivide(bihRootIdx, totalBoundary, N);
	printf("BIH constructed in %.2fms.\n", t.elapsed() * 1000);
}

std::string DumpBih(int i) {
	BihNode node = bihNode[i];
	std::string base = "{";
	if (node.type == 3) {
		base = base + "start=" + std::to_string(node.triangle);
		base = base + ", triIdx=" + std::to_string(triIdx[node.triangle]);
	}
	else {
		base = base + "lclip=" + std::to_string(node.clip[0]);
		base = base + ", rclip=" + std::to_string(node.clip[1]);
		switch (node.type)
		{
		case 0:
			base = base + ", type=X";
			break;
		case 1:
			base = base + ", type=Y";
			break;
		case 2:
			base = base + ", type=Z";
			break;
		default:
			base = base + ", type=" + std::to_string(node.type);
			break;
		}
		if (node.leftIdx != 0) {
			base = base + ", left=" + DumpBih(node.leftIdx);
			base = base + ", right=" + DumpBih(node.leftIdx + 1);
		}
	}
	return base + "}";
}

int IntersectBih(Ray& ray, aabb rootBoundary, int nodeIdx) {
	uint nodeStack[64];
	aabb boundaryStack[64];
	int stackPtr = 0;

	nodeStack[stackPtr] = bihRootIdx;
	boundaryStack[stackPtr] = rootBoundary;

	int counter = 0;

	while (stackPtr >= 0) {
		counter++;
		uint nodeI = nodeStack[stackPtr];
		BihNode node = bihNode[nodeI];
		aabb boundary = boundaryStack[stackPtr--];

		if (node.type == 3) {
			Tri triangle = tri[triIdx[node.triangle]];
			IntersectTri(ray, triangle);
			if (ray.t < 1e30f) {
				return counter;
			}
			else {
				continue;
			}
		}

		aabb leftBound = boundary;
		leftBound.bmax[node.type] = node.clip[0];
		aabb rightBound = boundary;
		rightBound.bmin[node.type] = node.clip[1];

		// TODO: Intersections can be easier probably
		float tLeft = IntersectAABB(ray, leftBound.bmin, leftBound.bmax);
		float tRight = IntersectAABB(ray, rightBound.bmin, rightBound.bmax);
	
		// Cases
		if (tLeft == 1e30f && tRight == 1e30f) {
			// No hit
			continue;
		} else if (tLeft != 1e30f && tRight == 1e30f) {
			// Only left hit
			nodeStack[++stackPtr] = node.leftIdx;
			boundaryStack[stackPtr] = leftBound;
		}
		else if (tLeft == 1e30f && tRight != 1e30f) {
			// Only right hit
			nodeStack[++stackPtr] = node.leftIdx + 1;
			boundaryStack[stackPtr] = rightBound;
		}
		else {
			// Double hit
			if (ray.D[node.type] > 0) {
				// Check left node first
				nodeStack[++stackPtr] = node.leftIdx + 1;
				boundaryStack[stackPtr] = rightBound;

				nodeStack[++stackPtr] = node.leftIdx;
				boundaryStack[stackPtr] = leftBound;
			}
			else {
				// Check right node first
				nodeStack[++stackPtr] = node.leftIdx;
				boundaryStack[stackPtr] = leftBound;

				nodeStack[++stackPtr] = node.leftIdx + 1;
				boundaryStack[stackPtr] = rightBound;

			}
		}
		
	}
	return counter;
}
#endif

//#define RECORD_TIME
#define RECORD_RAYS

#ifdef RECORD_TIME
const int CAMERA_FRAMES = 10;
float times[CAMERA_FRAMES];
#endif
#ifdef RECORD_RAYS
const int CAMERA_FRAMES = 1;
#endif

int camera_position = 0;
int counter = 0;
float total_time_ms;

void BihApp::Init()
{
#if SCENE == 0
	FILE* file = fopen("assets/unity.tri", "r");
	float a, b, c, d, e, f, g, h, i;
	for (int t = 0; t < N; t++)
	{
		fscanf(file, "%f %f %f %f %f %f %f %f %f\n",
			&a, &b, &c, &d, &e, &f, &g, &h, &i);
		tri[t].vertex0 = float3(a, b, c);
		tri[t].vertex1 = float3(d, e, f);
		tri[t].vertex2 = float3(g, h, i);
		#if METHOD == 1
		CheckBounds(tri[t].vertex0, tri[t].vertex1, tri[t].vertex2);
		#endif

	}
	fclose(file);
#elif SCENE == 1
	objl::Loader loader;
	loader.LoadFile(OBJ_PATH);
	vector<objl::Vertex> vertices = loader.LoadedMeshes[0].Vertices;
	vector<uint> indices = loader.LoadedMeshes[0].Indices;

	//printf("Size vertices: %d, Size indices: %d\n", vertices.size(), indices.size());
	// aabb (-32.857990, -11.336102, -32.637764 : 32.926079, 10.171690, 32.333691)
	//Size vertices: 65536, Size indices: 98304

	for (int i = 0, triidx = 0; i < N * 3; i += 3, triidx++)
	{
		//cout << vertices[i].Position.X << ", " << vertices[i].Position.Y << ", " << vertices[i].Position.Z << endl;
		tri[triidx].vertex0 = float3(vertices[indices[i]].Position.X, vertices[indices[i]].Position.Y, vertices[indices[i]].Position.Z);
		tri[triidx].vertex1 = float3(vertices[indices[i + 1]].Position.X, vertices[indices[i + 1]].Position.Y, vertices[indices[i + 1]].Position.Z);
		tri[triidx].vertex2 = float3(vertices[indices[i + 2]].Position.X, vertices[indices[i + 2]].Position.Y, vertices[indices[i + 2]].Position.Z);
#if METHOD == 1
		CheckBounds(tri[triidx].vertex0, tri[triidx].vertex1, tri[triidx].vertex2);
#endif
	}
#endif
	BuildBVH();
#if METHOD == 0
#elif METHOD == 1
	BuildGrid();
#elif METHOD == 2
	BuildBih();
	printf("aabb (%f, %f, %f : %f, %f, %f)\n", totalBoundary.bmin[0], totalBoundary.bmin[1], totalBoundary.bmin[2], totalBoundary.bmax[0], totalBoundary.bmax[1], totalBoundary.bmax[2]);
	//std::cout << DumpBih(0) << "\n";
#endif

#ifdef RECORD_TIME
	string filename;
	filename += OBJ_NAME;
	filename += "_";
	filename += METHOD_NAME;
	filename += "_time.csv";
	timefile.open(filename.c_str());
#endif

#ifdef RECORD_RAYS
	string filename;
	filename += OBJ_NAME;
	filename += "_";
	filename += METHOD_NAME;
	filename += "_rays.csv";
	rayfile.open(filename.c_str());
	rayfile << "camera position,ray_i,traversal steps,triangle intersections";
#endif
}

const int POSITIONS = 4;
#if SCENE == 0
float3 cameraPositions[POSITIONS] = {
	float3(-1.5f, -0.2f, -3.5f),

	float3(-4.5f, -0.2f, 0.0f),

	float3(-1.5f, -0.2f, 3.5f),

	float3(1.5f, -0.2f, 0.0f),
};

float3  cameraPoints[3 * POSITIONS] = {
	float3(-2.5f, 0.8f, -1.5f),
	float3(-0.5f, 0.8f, -1.5f),
	float3(-2.5f, -1.2f, -1.5f),

	float3(-2.5f, 0.8f, 1.0f),
	float3(-2.5f, 0.8f, -1.0f),
	float3(-2.5f, -1.2f, 1.0f),

	float3(-0.5f, 0.8f, 1.5f),
	float3(-2.5f, 0.8f, 1.5f),
	float3(-0.5f, -1.2f, 1.5f),

	float3(-0.5f, 0.8f, -1.0f),
	float3(-0.5f, 0.8f, 1.0f),
	float3(-0.5f, -1.2f, -1.0f),
};
#elif SCENE == 1
float3 cameraPositions[POSITIONS] = {
	float3(0.f, 2.f, -2.f),

	float3(-2.f, 2.f, 0.f),

	float3(0.f, 2.f, 2.f),

	float3(2.f, 2.f, 0.f),
};

float3  cameraPoints[3 * POSITIONS] = {
	cameraPositions[0] + float3(-2.f, 0.f, sqrt(2.0f)),
	cameraPositions[0] + float3(2.f, 0.f, sqrt(2.0f)),
	cameraPositions[0] + float3(-2.f, -sqrt(2.f), 0.f),

	cameraPositions[1] + float3(sqrt(2.0f), 0.f, 2.f),
	cameraPositions[1] + float3(sqrt(2.0f), 0.f, -2.f),
	cameraPositions[1] + float3(0.f, -sqrt(2.f), 2.f),

	cameraPositions[2] + float3(2.f, 0.f, -sqrt(2.0f)),
	cameraPositions[2] + float3(-2.f, 0.f, -sqrt(2.0f)),
	cameraPositions[2] + float3(2.f, -sqrt(2.f), 0.f),

	cameraPositions[3] + float3(-sqrt(2.0f), 0.f, -2.f),
	cameraPositions[3] + float3(-sqrt(2.0f), 0.f, 2.f),
	cameraPositions[3] + float3(0.f, -sqrt(2.f), -2.f),
};
#endif


void BihApp::Tick( float deltaTime )
{
	min_traversal_steps = MAXINT;
	min_triangle_intersections = MAXINT;
	max_traversal_steps = 0;
	max_triangle_intersects = 0;
	//camera_position = 3;
	// draw the scene
	screen->Clear(0);
	// define the corners of the screen in worldspace
	float3 p0 = cameraPoints[3 * camera_position];
	float3 p1 = cameraPoints[3 * camera_position + 1];
	float3 p2 = cameraPoints[3 * camera_position + 2];
	Ray ray;
	Timer t;
	// analysis
	int ind = 0;

	// render tiles of pixels
	for (int y = 0; y < SCRHEIGHT; y++) for (int x = 0; x < SCRWIDTH; x++)
	{
		tri_intrs_count = 0;
		// calculate the position of a pixel on the screen in worldspace
		float3 pixelPos = p0 + (p1 - p0) * (x / (float)SCRWIDTH) + (p2 - p0) * (y / (float)SCRHEIGHT);
		// define the ray in worldspace
		ray.O = cameraPositions[camera_position];
		ray.D = normalize(pixelPos - ray.O), ray.t = 1e30f;
		ray.rD = float3(1 / ray.D.x, 1 / ray.D.y, 1 / ray.D.z);
#if METHOD == 0
		int traversals = IntersectBVH(ray);
#elif METHOD == 1
		int traversals = IntersectGrid(ray);
#elif METHOD == 2
		int traversals = IntersectBih(ray, totalBoundary, 0);
#endif
		uint c = 500 - (int)(ray.t * 42);
		if (ray.t < 1e30f) screen->Plot(x, y, c * 0x10101);

		traversal_steps[ind] = traversals;
		triangle_intersects[ind] = tri_intrs_count;
		min_traversal_steps = min(min_traversal_steps, traversals);
		min_triangle_intersections = min(min_triangle_intersections, tri_intrs_count);
		max_traversal_steps = max(max_traversal_steps, traversals);
		max_triangle_intersects = max(max_triangle_intersects, tri_intrs_count);

#ifdef RECORD_RAYS
		rayfile << camera_position << "," << ind << "," << traversals << "," << tri_intrs_count << "\n";
#endif
		ind += 1;
	}
	float elapsed = t.elapsed() * 1000;
#ifdef RECORD_TIME
	times[counter] = elapsed;
#endif
	printf( "tracing time: %.2fms (%5.2fK rays/s)\n", elapsed, sqr( 630 ) / elapsed );

	counter = (counter + 1) % CAMERA_FRAMES;
	if (counter == 0) {
#ifdef RECORD_TIME
		timefile << "Position " << camera_position << ": ";
		for (int i = 0; i < CAMERA_FRAMES-1; i++) {
			timefile << times[i] << ", ";
		}
		timefile << times[CAMERA_FRAMES - 1] << "\n";
#endif

		std::string intersect_filename;
		intersect_filename += OBJ_NAME;
		intersect_filename += "_";
		intersect_filename += METHOD_NAME;
		intersect_filename += "_CAM";
		intersect_filename += to_string(camera_position);
		intersect_filename += "_intersects.png";
		for (int i = 0; i < SCRHEIGHT * SCRWIDTH; i++) {
			intrs_data[3 * i] = 0;
			intrs_data[3 * i + 1] = (int)((255.0f / (float)max_triangle_intersects) * triangle_intersects[i]);
			intrs_data[3 * i + 2] = 0;
		}
		stbi_write_png(intersect_filename.c_str(), SCRWIDTH, SCRHEIGHT, 3, intrs_data, (int)(sizeof(uint8_t) * SCRWIDTH * 3));

		std::string traversal_filename;
		traversal_filename += OBJ_NAME;
		traversal_filename += "_";
		traversal_filename += METHOD_NAME;
		traversal_filename += "_CAM";
		traversal_filename += to_string(camera_position);
		traversal_filename += "_traversals.png";
		for (int i = 0; i < SCRHEIGHT * SCRWIDTH; i++) {
			intrs_data[3 * i] = (int)((255.0f / (float)max_traversal_steps) * traversal_steps[i]);
			intrs_data[3 * i + 1] = 0;
			intrs_data[3 * i + 2] = 0;
		}
		stbi_write_png(traversal_filename.c_str(), SCRWIDTH, SCRHEIGHT, 3, intrs_data, (int)(sizeof(uint8_t) * SCRWIDTH * 3));


#if defined(RECORD_TIME) || defined(RECORD_RAYS)
		if (camera_position + 1 == POSITIONS) {
			timefile.close();
			rayfile.close();
			exit(-1);
		}
		else {
			camera_position += 1;
		}
#else
		camera_position = (camera_position + 1) % POSITIONS;
#endif
	}
}

// EOF