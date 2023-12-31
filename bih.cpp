#include "precomp.h"
#include "bih.h"
#include "stb_image_write.h"
#include <format>
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

// triangle count
#define N	12582 // hardcoded for the unity vehicle mesh

// forward declarations
void Subdivide( uint nodeIdx );
void UpdateNodeBounds( uint nodeIdx );

// minimal structs
struct Tri { float3 vertex0, vertex1, vertex2; float3 centroid; };
struct BihNode {
	int type; // Lowest bits: axis (00, 01, 10) or leaf(11)
	int leftIdx; // Node only
	union {
		float clip[2];
		int triangle;
	};
};
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

// application data
Tri tri[N];
uint triIdx[N];
BihNode bihNode[N * 2];
uint bihRootIdx = 0, nodesUsed = 2;
aabb totalBoundary;

// Analysis data
int traversal_steps[SCRHEIGHT * SCRWIDTH];
int triangle_intersects[SCRHEIGHT * SCRWIDTH];
int tri_intrs_count;
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

void BihApp::Init()
{
	FILE* file = fopen( "assets/unity.tri", "r" );
	float a, b, c, d, e, f, g, h, i;
	for (int t = 0; t < N; t++)
	{
		fscanf( file, "%f %f %f %f %f %f %f %f %f\n",
			&a, &b, &c, &d, &e, &f, &g, &h, &i );
		tri[t].vertex0 = float3( a, b, c );
		tri[t].vertex1 = float3( d, e, f );
		tri[t].vertex2 = float3( g, h, i );
	}
	fclose( file );
	// intialize a scene with N random triangles
	/*for (int i = 0; i < N; i++)
	{
		float3 r0 = float3(RandomFloat(), RandomFloat(), RandomFloat());
		float3 r1 = float3(RandomFloat(), RandomFloat(), RandomFloat());
		float3 r2 = float3(RandomFloat(), RandomFloat(), RandomFloat());
		tri[i].vertex0 = r0 * 9 - float3(5);
		tri[i].vertex1 = tri[i].vertex0 + r1, tri[i].vertex2 = tri[i].vertex0 + r2;
	}*/
	/*tri[0] = Tri{ 
		float3(17.0, 29.0, 0.0),
		float3(93.0, 47.0, 0.0),
		float3(32.0, 92.0, 0.0),
	};
	tri[1] = Tri{
		float3(255.0, 258.0, 0.0),
		float3(380.0, 316.0, 0.0),
		float3(302.0, 543.0, 0.0),
	};
	tri[2] = Tri{
		float3(380.0, 316.0, 0.0),
		float3(255.0, 258.0, 0.0),
		float3(348.0, 166.0, 0.0),
	};
	tri[3] = Tri{
		float3(33.0, 378.0, 0.0),
		float3(229.0, 436.0, 0.0),
		float3(154.0, 572.0, 0.0),
	};
	tri[4] = Tri{
		float3(32.0, 92.0, 0.0),
		float3(93.0, 47.0, 0.0),
		float3(19.0, 271.0, 0.0),
	};
	tri[5] = Tri{
		float3(380.0, 316.0, 0.0),
		float3(348.0, 166.0, 0.0),
		float3(452.0, 241.0, 0.0),
	};*/
	// construct the BVH
	BuildBih();
	//std::cout << DumpBih(0) << "\n";
}

const int POSITIONS = 4;
const int CAMERA_FRAMES = 20;
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


void BihApp::Tick( float deltaTime )
{
	max_traversal_steps = 0;
	max_triangle_intersects = 0;
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
		int traversals = IntersectBih(ray, totalBoundary, 0);
		uint c = 500 - (int)(ray.t * 42);
		if (ray.t < 1e30f) screen->Plot(x, y, c * 0x10101);

		traversal_steps[ind++] = traversals;
		triangle_intersects[ind++] = tri_intrs_count;
		max_traversal_steps = max(max_traversal_steps, traversals);
		max_triangle_intersects = max(max_triangle_intersects, tri_intrs_count);
	}
	float elapsed = t.elapsed() * 1000;
	printf( "tracing time: %.2fms (%5.2fK rays/s)\n", elapsed, sqr( 630 ) / elapsed );

	/*if (counter % CAMERA_FRAMES == 1) {
		std::string filename;
		filename += "intersection-count";
		filename += to_string(camera_position);
		filename += ".png";
		stbi_write_png(filename.c_str(), SCRWIDTH, SCRHEIGHT, 3, intrs_data, (int)(sizeof(uint8_t) * SCRWIDTH * 3));
	}*/
}

// EOF