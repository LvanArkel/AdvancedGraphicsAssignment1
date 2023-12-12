#include "precomp.h"
#include "basics.h"
#include "stb_image_write.h"
#include "OBJ_Loader.h"

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
#define N	12582	
#define FLOAT_MAX  3.402823466e+38
#define FLOAT_MIN  1.175494351e-38
#define GRID_SIZE  64

#define OBJ_PATH "data/scene1/suzanne2.obj"

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

struct IntersectionData {
public:
	int triIdx; //triangle index of intersection with ray
	float t; //Ray's t parameter value at time of intersection
};

Grid grid;

// application data bvh
Tri tri[N];
uint triIdx[N];
BVHNode bvhNode[N * 2];
uint rootNodeIdx = 0, nodesUsed = 1;

const int POSITIONS = 4;
const int CAMERA_FRAMES = 5;
int camera_position = 0;
int counter = 0;

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

//analyzation data
uint8_t intrs_data[SCRHEIGHT * SCRWIDTH * 3]; //intersection count data for image generation
int tri_intrs_count = 0; //number of times per ray a triangle intersection check was done
// functions

void IntersectTri(Ray& ray, const Tri& tri)
{
	//tri_intrs_count++;
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

	// Check for parallel rays
	if (ray.D.x == 0) if (ray.O.x < grid.min.x || ray.O.x > grid.max.x) return false; // No intersection
	if (ray.D.y == 0) if (ray.O.y < grid.min.y || ray.O.y > grid.max.y) return false; // No intersection
	if (ray.D.z == 0) if (ray.O.z < grid.min.z || ray.O.z > grid.max.z) return false; // No intersection

	// Calculate t at intersection for each axis
	neart.x = (grid.min.x - ray.O.x) / ray.D.x;
	fart.x = (grid.max.x - ray.O.x) / ray.D.x;
	neart.y = (grid.min.y - ray.O.y) / ray.D.y;
	fart.y = (grid.max.y - ray.O.y) / ray.D.y;
	neart.z = (grid.min.z - ray.O.z) / ray.D.z;
	fart.z = (grid.max.z - ray.O.z) / ray.D.z;

	// pick the closest and the furthest
	tNear = max({ min(neart.x, fart.x), min(neart.y, fart.y), min(neart.z, fart.z) });
	tFar = min({ max(neart.x, fart.x), max(neart.y, fart.y), max(neart.z, fart.z) });

	// Check if there is any valid intersection
	if (tNear > tFar || tFar < 0) {
		return false; // No intersection
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

void IntersectGrid(Ray& ray) {
	float3 intersectionpoint, exitpoint;
	if (!FindRayGridIntersections(ray, intersectionpoint, exitpoint)) return;

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

		tri_intrs_count++;
		t = min({ tMaxX, tMaxY, tMaxZ });
		if (CheckGridCell(x0, y0, z0, ray, intersected_triangles)) {
			if (ray.t < t) return;
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
		if (x0 >= GRID_SIZE || x0 < 0 || y0 >= GRID_SIZE || y0 < 0 || z0 >= GRID_SIZE || z0 < 0) return;

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

void BasicBVHApp::Init()
{
	//objl::Loader loader;
	//loader.LoadFile(OBJ_PATH);
	//vector<objl::Vertex> vertices = loader.LoadedMeshes[0].Vertices;

	//int triidx = 0;
	////for (int i = 0; i < N * 3; i += 3)
	//for (int i = 0; i < N * 3; i += 3)
	//{
	//	//cout << vertices[i].Position.X << ", " << vertices[i].Position.Y << ", " << vertices[i].Position.Z << endl;
	//	tri[triidx].vertex0 = float3(vertices[i].Position.X, vertices[i].Position.Y, vertices[i].Position.Z);
	//	tri[triidx].vertex1 = float3(vertices[i + 1].Position.X, vertices[i + 1].Position.Y, vertices[i + 1].Position.Z);
	//	tri[triidx].vertex2 = float3(vertices[i + 2].Position.X, vertices[i + 2].Position.Y, vertices[i + 2].Position.Z);

	//	CheckBounds(tri[triidx].vertex0, tri[triidx].vertex1, tri[triidx].vertex2);

	//	triidx++;
	//}

	//std::cout << "grid min; " << grid.min.x << ", " << grid.min.y << ", " << grid.min.z << std::endl;
	//std::cout << "grid max; " << grid.max.x << ", " << grid.max.y << ", " << grid.max.z << std::endl;
	//BuildGrid();

	FILE* file = fopen("assets/unity.tri", "r");
	float a, b, c, d, e, f, g, h, i;
	for (int t = 0; t < N; t++)
	{
		fscanf(file, "%f %f %f %f %f %f %f %f %f\n",
			&a, &b, &c, &d, &e, &f, &g, &h, &i);
		tri[t].vertex0 = float3(a, b, c);
		tri[t].vertex1 = float3(d, e, f);
		tri[t].vertex2 = float3(g, h, i);
		CheckBounds(tri[t].vertex0, tri[t].vertex1, tri[t].vertex2);

	}
	fclose(file);
	BuildGrid();
}


void TickGrid(Tmpl8::Surface* screen) {
	if (counter % CAMERA_FRAMES == 0) {
		camera_position = (camera_position + 1) % POSITIONS;
	}
	//camera_position = 1;
	counter += 1;
	// draw the scene
	screen->Clear(0);
	// define the corners of the screen in worldspace
	float3 p0 = cameraPoints[3 * camera_position];
	float3 p1 = cameraPoints[3 * camera_position + 1];
	float3 p2 = cameraPoints[3 * camera_position + 2];
	Ray ray;
	Timer t;


	int count = 0;
	int ind = 0;


	for (int y = 0; y < SCRHEIGHT; y += 1) for (int x = 0; x < SCRWIDTH; x += 1)
	{
		// calculate the position of a pixel on the screen in worldspace
		float3 pixelPos = p0 + (p1 - p0) * (x / (float)SCRWIDTH) + (p2 - p0) * (y / (float)SCRHEIGHT);
		// define the ray in worldspace
		ray.O = cameraPositions[camera_position];
		ray.D = normalize(pixelPos - ray.O);
		// initially the ray has an 'infinite length'
		ray.t = 1e30f;

		//for (int i = 0; i < N; i++) IntersectTri(ray, tri[i]);

		IntersectGrid(ray);

		//intrs_data[ind++] = 0, intrs_data[ind++] = (int)((255.0f / (N * 0.1f)) * tri_intrs_count), intrs_data[ind++] = 0;
		intrs_data[ind++] = (int)((255.0f / (N * 0.005f)) * tri_intrs_count), intrs_data[ind++] = 0, intrs_data[ind++] = 0;
		tri_intrs_count = 0;
		uint c = 500 - (int)(ray.t * 42);
		if (ray.t < 1e30f) screen->Plot(x, y, c * 0x10101);
	}

	std::cout << count << std::endl;
	float elapsed = t.elapsed() * 1000;
	printf("tracing time: %.2fms (%5.2fK rays/s)\n", elapsed, sqr(630) / elapsed);
	stbi_write_png("intsection-count.png", SCRWIDTH, SCRHEIGHT, 3, intrs_data, (int)(sizeof(uint8_t) * SCRWIDTH * 3));
}

void BasicBVHApp::Tick(float deltaTime)
{
	TickGrid(screen);
}

// EOF