#include "precomp.h"
#include "basics.h"

// THIS SOURCE FILE:
// Code for the article "How to Build a BVH", part 1: basics. Link:
// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics
// This is bare-bones BVH construction and traversal code, running in
// a minimalistic framework.
// Feel free to copy this code to your own framework. Absolutely no
// rights are reserved. No responsibility is accepted either.
// For updates, follow me on twitter: @j_bikker.

TheApp* CreateApp() { return new BasicBVHApp(); }

// triangle count
#define N	16
#define FLOAT_MAX  3.402823466e+38
#define FLOAT_MIN  1.175494351e-38
#define GRID_SIZE  4

#define GRID_ACC 1

// forward declarations
void Subdivide(uint nodeIdx);
void UpdateNodeBounds(uint nodeIdx);

// minimal structs
struct Tri { float3 vertex0, vertex1, vertex2; float3 centroid; };
__declspec(align(32)) struct BVHNode
{
	float3 aabbMin, aabbMax;
	uint leftFirst, triCount;
	bool isLeaf() { return triCount > 0; }
};
struct Ray { float3 O, D; float t = 1e30f; };

// application data grid
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

Grid grid;

// application data bvh
Tri tri[N];
uint triIdx[N];
BVHNode bvhNode[N * 2];
uint rootNodeIdx = 0, nodesUsed = 1;

// functions

void IntersectTri(Ray& ray, const Tri& tri)
{
	const float3 edge1 = tri.vertex1 - tri.vertex0;
	const float3 edge2 = tri.vertex2 - tri.vertex0;
	const float3 h = cross(ray.D, edge2);
	const float a = dot(edge1, h);
	if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
	const float f = 1 / a;
	const float3 s = ray.O - tri.vertex0;
	const float u = f * dot(s, h);
	if (u < 0 || u > 1) return;
	const float3 q = cross(s, edge1);
	const float v = f * dot(ray.D, q);
	if (v < 0 || u + v > 1) return;
	const float t = f * dot(edge2, q);
	if (t > 0.0001f) ray.t = min(ray.t, t);
}

bool IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax)
{
	float tx1 = (bmin.x - ray.O.x) / ray.D.x, tx2 = (bmax.x - ray.O.x) / ray.D.x;
	float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
	float ty1 = (bmin.y - ray.O.y) / ray.D.y, ty2 = (bmax.y - ray.O.y) / ray.D.y;
	tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
	float tz1 = (bmin.z - ray.O.z) / ray.D.z, tz2 = (bmax.z - ray.O.z) / ray.D.z;
	tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
	return tmax >= tmin && tmin < ray.t && tmax > 0;
}

void IntersectBVH(Ray& ray, const uint nodeIdx)
{
	BVHNode& node = bvhNode[nodeIdx];
	if (!IntersectAABB(ray, node.aabbMin, node.aabbMax)) return;
	if (node.isLeaf())
	{
		for (uint i = 0; i < node.triCount; i++)
			IntersectTri(ray, tri[triIdx[node.leftFirst + i]]);
	}
	else
	{
		IntersectBVH(ray, node.leftFirst);
		IntersectBVH(ray, node.leftFirst + 1);
	}
}

void BuildBVH()
{
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
	Subdivide(rootNodeIdx);
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
		node.aabbMin = fminf(node.aabbMin, leafTri.vertex0),
			node.aabbMin = fminf(node.aabbMin, leafTri.vertex1),
			node.aabbMin = fminf(node.aabbMin, leafTri.vertex2),
			node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex0),
			node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex1),
			node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex2);
	}
}

void Subdivide(uint nodeIdx)
{
	// terminate recursion
	BVHNode& node = bvhNode[nodeIdx];
	if (node.triCount <= 2) return;
	// determine split axis and position
	float3 extent = node.aabbMax - node.aabbMin;
	int axis = 0;
	if (extent.y > extent.x) axis = 1;
	if (extent.z > extent[axis]) axis = 2;
	float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
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

void CheckGridCell(int x, int y, int z, Ray& ray) {
	GridCell gc = grid.grid[x][y][z];
	vector<int> triangles = gc.triangles;

	for (int i = 0; i < triangles.size(); i++) {
		IntersectTri(ray, tri[triangles[i]]);
		if (ray.t < 1e30f) {
			//std::cout << "SLOEP" << std::endl;
			return; 
		}
	}
}

void IntersectGird(Ray& ray) {
	float tMin = (grid.min.x - ray.O.x) / ray.D.x;
	float tMax = (grid.max.x - ray.O.x) / ray.D.x;

	if (tMin > tMax) std::swap(tMin, tMax);

	float txMin = tMin;

	float tyMin = (grid.min.y - ray.O.y) / ray.D.y;
	float tyMax = (grid.max.y - ray.O.y) / ray.D.y;

	if (tyMin > tyMax) std::swap(tyMin, tyMax);

	if ((tMin > tyMax) || (tyMin > tMax)) {
		return; // No intersection with the XY faces
	}

	if (tyMin > tMin) tMin = tyMin;
	if (tyMax < tMax) tMax = tyMax;

	float tzMin = (grid.min.z - ray.O.z) / ray.D.z;
	float tzMax = (grid.max.z - ray.O.z) / ray.D.z;

	if (tzMin > tzMax) std::swap(tzMin, tzMax);

	if ((tMin > tzMax) || (tzMin > tMax)) {
		return; // No intersection with the box
	}

	if (tzMin > tMin) tMin = tzMin;
	// if (tzMax < tMax) tMax = tzMax; // Uncomment if you want to use tMax for something

	float3 intersectionPoint = float3(0);
	intersectionPoint.x = ray.O.x + tMin * ray.D.x;
	intersectionPoint.y = ray.O.y + tMin * ray.D.y;
	intersectionPoint.z = ray.O.z + tMin * ray.D.z;

	int cellIndex[] = { 0, 0, 0 };
	cellIndex[0] = floor((intersectionPoint.x - grid.min.x) / grid.cellSize.x);
	cellIndex[1] = floor((intersectionPoint.y - grid.min.y) / grid.cellSize.y);
	cellIndex[2] = floor((intersectionPoint.z - grid.min.z) / grid.cellSize.z);

	CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);

	float3 deltaT = float3(0);
	//QUESTION isn't dividing by the ray direction unsafe? since it might be 0
	deltaT.x = ray.D.x < 0 ? -grid.size.x / ray.D.x : grid.size.x / ray.D.x;
	deltaT.y = ray.D.y < 0 ? -grid.size.y / ray.D.y : grid.size.y / ray.D.y;
	deltaT.z = ray.D.z < 0 ? -grid.size.z / ray.D.z : grid.size.z / ray.D.z;
	float tx = deltaT.x;
	float ty = deltaT.y;
	float tz = deltaT.z;
	//float tx = txMin == tMin ? deltaT.x : 0.0f;
	//float ty = tyMin == tMin ? deltaT.y : 0.0f;
	//float tz = tzMin == tMin ? deltaT.z : 0.0f;
	//ray.t = 0.0f;
	//return;

	int counter = 0;
	float t = 0;
	//std::cout << "AAA" << std::endl;
	while (1) {
		//std::cout << cellIndex[0] << ", " << cellIndex[1] << ", " << cellIndex[2] << std::endl;
		counter++;
		if (tx < ty && tx < tz) {
			//std::cout << "x" << std::endl;
			t = tx;
			tx += deltaT.x;
			cellIndex[0] = ray.D.x < 0 ? cellIndex[0] - 1 : cellIndex[0] + 1;
			if (cellIndex[0] >= 0 && cellIndex[0] < GRID_SIZE) CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);
			else break;
		}
		else if (ty < tx && ty < tz) {
			//std::cout << "y" << std::endl;
			t = ty;
			ty += deltaT.y;
			cellIndex[1] = ray.D.y < 0 ? cellIndex[1] - 1 : cellIndex[1] + 1;
			if (cellIndex[1] >= 0 && cellIndex[1] < GRID_SIZE) CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);
			else break;
		}
		else {
			//std::cout << "z" << std::endl;
			t = tz;
			tz += deltaT.z;
			cellIndex[2] = ray.D.z < 0 ? cellIndex[2] - 1 : cellIndex[2] + 1;
			if (cellIndex[2] > -1 && cellIndex[2] < GRID_SIZE)CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);
			else break;
		}
		if (counter == 100) break;
	}
	
}

void IntersectGrid(Ray& ray) {
	float tx, ty, tz;
	int cellIndex[] = { 0, 0, 0 };
	float3 mint = float3((grid.min.x - ray.O.x) / ray.D.x,
		(grid.min.y - ray.O.y) / ray.D.y,
		(grid.min.z - ray.O.z) / ray.D.z);
	float3 maxt = float3((grid.max.x - ray.O.x) / ray.D.x,
		(grid.max.y - ray.O.y) / ray.D.y,
		(grid.max.z - ray.O.z) / ray.D.z);
	float tmin = max(max(mint.x, mint.y), mint.z);
	float tmax = min(min(maxt.x, maxt.y), maxt.z);
	//std::cout << tmin << std::endl;
	//std::cout << tmax << std::endl;

	if (tmin > tmax) {
		//std::cout << "SLOEP" << std::endl;
		return; //no intersection
	}
	float3 intrspoint = ray.O + (tmin * ray.D);
	cellIndex[0] = floor((intrspoint.x - grid.min.x) / grid.cellSize.x);
	cellIndex[1] = floor((intrspoint.y - grid.min.y) / grid.cellSize.y);
	cellIndex[2] = floor((intrspoint.z - grid.min.z) / grid.cellSize.z);
	ray.t = 0.0f;


	
	std::cout << "A" << std::endl;
	if (cellIndex[0] >= 0 && cellIndex[0] <= GRID_SIZE - 1) {
		if (cellIndex[1] >= 0 && cellIndex[1] <= GRID_SIZE - 1) {
			if (cellIndex[2] >= 0 && cellIndex[2] <= GRID_SIZE - 1) {
				CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);
			}
		}
	}

	float3 deltaT = float3(0);
	//QUESTION isn't dividing by the ray direction unsafe? since it might be 0
	deltaT.x = ray.D.x < 0 ? -grid.size.x * (1.0 / ray.D.x) : grid.size.x * (1.0 / ray.D.x);
	deltaT.y = ray.D.y < 0 ? -grid.size.y * (1.0 / ray.D.y) : grid.size.y * (1.0 / ray.D.y);
	deltaT.z = ray.D.z < 0 ? -grid.size.z * (1.0 / ray.D.z) : grid.size.z * (1.0 / ray.D.z);

	tx = deltaT.x, ty = deltaT.y, tz = deltaT.z;

	

	/*std::cout << ray.D.x << std::endl;
	std::cout << ray.D.y << std::endl;
	std::cout << ray.D.z << std::endl;*/


	//float3 rayOriginGrid = ray.O - grid.min;
	//float3 oCell = rayOriginGrid / grid.cellSize;
	//std::cout << oCell.x << ", " << oCell.y << ", " << oCell.z << std::endl;
	//tx = (grid.cellSize.x - ray.O.x) / ray.D.x;
	//ty = (grid.cellSize.y - ray.O.y) / ray.D.y;
	//tz = (grid.cellSize.z - ray.O.z) / ray.D.z;
	/*if (ray.D.x < 0) {
		deltaT.x = -grid.cellSize.x / ray.D.x;
		tx = floor(oCell.x) * grid.cellSize.x - ray.O.x / ray.D.x;
	}
	else {
		deltaT.x = -grid.cellSize.x / ray.D.x;
		tx = (floor(oCell.x) + 1.0f) * grid.cellSize.x - ray.O.x / ray.D.x;
	}
	if (ray.D.y < 0) {
		deltaT.y = -grid.cellSize.y / ray.D.y;
		ty = floor(oCell.y) * grid.cellSize.y - ray.O.y / ray.D.y;
	}
	else {
		deltaT.y = -grid.cellSize.y / ray.D.y;
		ty = (floor(oCell.y) + 1.0f) * grid.cellSize.y - ray.O.y / ray.D.y;
	}
	if (ray.D.z < 0) {
		deltaT.z = -grid.cellSize.z / ray.D.z;
		tz = floor(oCell.z) * grid.cellSize.z - ray.O.z / ray.D.z;
	}
	else {
		deltaT.z = -grid.cellSize.z / ray.D.z;
		tz = (floor(oCell.z) + 1.0f) * grid.cellSize.z - ray.O.z / ray.D.z;
	}*/

	//deltaT = { grid.cellSize.x / abs(ray.D.x), grid.cellSize.y / abs(ray.D.y), grid.cellSize.z / abs(ray.D.z) };
	//float tx = (((floor(oCell.x) + 1.0f) * grid.cellSize.x) - ray.O.x) / ray.D.x;
	//float ty = (((floor(oCell.y) + 1.0f) * grid.cellSize.y) - ray.O.y) / ray.D.y;
	//float tz = (((floor(oCell.z) + 1.0f) * grid.cellSize.z) - ray.O.z) / ray.D.z;
	//float t = 0;
	//int cellIndex[] = { floor(tx / grid.size.x), floor(ty / grid.size.y), floor(tz / grid.size.z) };
	//std::cout << tx << ", " << ty << ", " << tz << std::endl;
	//std::cout << cellIndex[0] << ", " << cellIndex[1] << ", " << cellIndex[2] << std::endl;

	//int counter = 0;
	//float t = 0;
	//std::cout << "AAA" << std::endl;
	//while (1) {
		//std::cout << cellIndex[0] << ", " << cellIndex[1] << ", " << cellIndex[2] << std::endl;
		//counter++;
		//if (tx < ty && tx < tz) {
		//	//std::cout << "x" << std::endl;
		//	t = tx;
		//	tx += deltaT.x;
		//	cellIndex[0] = ray.D.x < 0 ? cellIndex[0] - 1 : cellIndex[0] + 1;
		//	if (cellIndex[0] >= 0 && cellIndex[0] < GRID_SIZE) CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);
		//	else break;
		//}
		//else if (ty < tx && ty < tz) {
		//	//std::cout << "y" << std::endl;
		//	t = ty;
		//	ty += deltaT.y;
		//	cellIndex[1] = ray.D.y < 0 ? cellIndex[1] - 1 : cellIndex[1] + 1;
		//	if (cellIndex[1] >= 0 && cellIndex[1] < GRID_SIZE) CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);
		//	else break;
		//}
		//else {
			//std::cout << "z" << std::endl;	
			//t = tz;
			//tz += deltaT.z;
			//cellIndex[2] = ray.D.z < 0 ? cellIndex[2] - 1 : cellIndex[2] + 1;
			//if (cellIndex[2] > -1 && cellIndex[2] < GRID_SIZE)CheckGridCell(cellIndex[0], cellIndex[1], cellIndex[2], ray);
			//else break;
		//}
		//if (counter == 100) break;
	//}
	//std::cout << cellIndex[0] << ", " << cellIndex[1] << ", " << cellIndex[2] << std::endl;
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
		float xmin = min({ tri[i].vertex0.x, tri[i].vertex1.x, tri[i].vertex2.x });
		float xmax = max({ tri[i].vertex0.x, tri[i].vertex1.x, tri[i].vertex2.x });
		float ymin = min({ tri[i].vertex0.y, tri[i].vertex1.y, tri[i].vertex2.y });
		float ymax = max({ tri[i].vertex0.y, tri[i].vertex1.y, tri[i].vertex2.y });
		float zmin = min({ tri[i].vertex0.z, tri[i].vertex1.z, tri[i].vertex2.z });
		float zmax = max({ tri[i].vertex0.z, tri[i].vertex1.z, tri[i].vertex2.z });

		//	std::cout << grid.cellSize.x << ", " << grid.cellSize.y << ", " << grid.cellSize.z << std::endl;


		for (float x = xmin - grid.min.x; x <= xmax - grid.min.x; x += grid.cellSize.x) {
			for (float y = ymin - grid.min.y; y <= ymax - grid.min.y; y += grid.cellSize.y) {
				for (float z = zmin - grid.min.z; z <= zmax - grid.min.z; z += grid.cellSize.z) {
					int xi = floor(x / grid.cellSize.x);
					int yi = floor(y / grid.cellSize.y);
					int zi = floor(z / grid.cellSize.z);
					grid.grid[xi][yi][zi].triangles.push_back(i);
				}
			}
		}
	}

	for (int x = 0; x < GRID_SIZE; x++) {
		for (int y = 0; y < GRID_SIZE; y++) {
			for (int z = 0; z < GRID_SIZE; z++) {
				vector<int> triidxs = grid.grid[x][y][z].triangles;
				for (int i = 0; i < triidxs.size(); i++) {
					std::cout << triidxs.size() << std::endl;
				}
			}
		}
	}
}

void CheckBounds(float3 v0, float3 v1, float3 v2) {
	float minx = min(min(v0.x, v1.x), v2.x), maxx = max(max(v0.x, v1.x), v2.x);
	float miny = min(min(v0.y, v1.y), v2.y), maxy = max(max(v0.y, v1.y), v2.y);
	float minz = min(min(v0.z, v1.z), v2.z), maxz = max(max(v0.z, v1.z), v2.z);

	if (minx < grid.min.x) grid.min.x = minx;
	if (miny < grid.min.y) grid.min.y = miny;
	if (minz < grid.min.z) grid.min.z = minz;
	if (maxx > grid.max.x) grid.max.x = maxx;
	if (maxy > grid.max.y) grid.max.y = maxy;
	if (maxz > grid.max.z) grid.max.z = maxz;

	//if (v0.x < grid.min.x) grid.min.x = v0.x; if (v1.x < grid.min.x) grid.min.x = v1.x; if (v2.x < grid.min.x) grid.min.x = v2.x;
	//if (v0.x > grid.max.x) grid.max.x = v0.x; if (v1.x > grid.max.x) grid.max.x = v1.x; if (v2.x > grid.max.x) grid.max.x = v2.x;
	//if (v0.y < grid.min.y) grid.min.y = v0.y; if (v1.y < grid.min.y) grid.min.y = v1.y; if (v2.y < grid.min.y) grid.min.y = v2.y;
	//if (v0.y > grid.max.y) grid.max.y = v0.y; if (v1.y > grid.max.y) grid.max.y = v1.y; if (v2.y > grid.max.y) grid.max.y = v2.y;
	//if (v0.z < grid.min.z) grid.min.z = v0.z; if (v1.z < grid.min.z) grid.min.z = v1.z; if (v2.z < grid.min.z) grid.min.z = v2.z;
	//if (v0.z > grid.max.z) grid.max.z = v0.z; if (v1.z > grid.max.z) grid.max.z = v1.z; if (v2.z > grid.max.z) grid.max.z = v2.z;
}

void BasicBVHApp::Init()
{
	// intialize a scene with N random triangles
	for (int i = 0; i < N; i++)
	{
		float3 r0 = float3(RandomFloat(), RandomFloat(), RandomFloat());
		float3 r1 = float3(RandomFloat(), RandomFloat(), RandomFloat());
		float3 r2 = float3(RandomFloat(), RandomFloat(), RandomFloat());
		tri[i].vertex0 = r0 * 9 - float3(5);
		tri[i].vertex1 = tri[i].vertex0 + r1, tri[i].vertex2 = tri[i].vertex0 + r2;
		//std::cout << tri[i].vertex0.x << ", " << tri[i].vertex0.y << ", " << tri[i].vertex0.z << std::endl;
		CheckBounds(tri[i].vertex0, tri[i].vertex1, tri[i].vertex2);
	}

	std::cout << grid.min.x << ", " << grid.min.y << ", " << grid.min.z << std::endl;
	std::cout << grid.max.x << ", " << grid.max.y << ", " << grid.max.z << std::endl;
#if GRID_ACC:
	BuildGrid();
#else
	// construct the BVH
	BuildBVH();
#endif
}

void TickBVH(Tmpl8::Surface* screen) {
	// draw the scene
	screen->Clear(0);
	// define the corners of the screen in worldspace
	float3 p0(-1, 1, -15), p1(1, 1, -15), p2(-1, -1, -15);
	Ray ray;
	Timer t;
	for (int y = 0; y < SCRHEIGHT; y++) for (int x = 0; x < SCRWIDTH; x++)
	{
		// calculate the position of a pixel on the screen in worldspace
		float3 pixelPos = p0 + (p1 - p0) * (x / (float)SCRWIDTH) + (p2 - p0) * (y / (float)SCRHEIGHT);
		// define the ray in worldspace
		ray.O = float3(0, 0, -18);
		ray.D = normalize(pixelPos - ray.O);
		// initially the ray has an 'infinite length'
		ray.t = 1e30f;
#if 0
		for (int i = 0; i < N; i++) IntersectTri(ray, tri[i]);
#else
		IntersectBVH(ray, rootNodeIdx);
#endif
		if (ray.t < 1e30f) screen->Plot(x, y, 0xffffff);
	}
	float elapsed = t.elapsed() * 1000;
	printf("tracing time: %.2fms (%5.2fK rays/s)\n", elapsed, sqr(630) / elapsed);
}

void TickGrid(Tmpl8::Surface* screen) {
	// draw the scene
	screen->Clear(0);
	// define the corners of the screen in worldspace
	float3 p0(-1, 1, -15), p1(1, 1, -15), p2(-1, -1, -15);
	Ray ray;
	Timer t;
	std::cout << "start" << std::endl;
	//int ct = 0;
	int count = 0;
	for (int y = 0; y < SCRHEIGHT; y += 1) for (int x = 0; x < SCRWIDTH; x += 1)
	{
		// calculate the position of a pixel on the screen in worldspace
		float3 pixelPos = p0 + (p1 - p0) * (x / (float)SCRWIDTH) + (p2 - p0) * (y / (float)SCRHEIGHT);
		// define the ray in worldspace
		ray.O = float3(0, 0, -16);
		ray.D = normalize(pixelPos - ray.O);
		// initially the ray has an 'infinite length'
		ray.t = 1e30f;

		//for (int i = 0; i < N; i++) IntersectTri(ray, tri[i]);
		IntersectGird(ray);
		//for (int x = 0; x < GRID_SIZE; x++) {
		//	//std::cout << x << std::endl;
		//	for (int y = 0; y < GRID_SIZE; y++) {
		//		for (int z = 0; z < GRID_SIZE; z++) {
		//			CheckGridCell(x, y, z, ray);
		//			//vector<int> triidxs = grid.grid[x][y][z].triangles;
		//			//for (int i = 0; i < triidxs.size(); i++) {
		//			//	
		//			//}
		//		}
		//	}
		//}

		if (ray.t < 1e30f) screen->Plot(x, y, 0xffffff);
	}
	std::cout << count << std::endl;
	float elapsed = t.elapsed() * 1000;
	printf("tracing time: %.2fms (%5.2fK rays/s)\n", elapsed, sqr(630) / elapsed);
}

void BasicBVHApp::Tick(float deltaTime)
{
#if GRID_ACC:
	TickGrid(screen);
#else
	TickBVH(screen);
#endif
}

// EOF