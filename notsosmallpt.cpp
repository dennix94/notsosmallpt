#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <random>
#define M_PI 3.1415926535897932384626433832795
#define M_1_PI 0.31830988618379067153776752674503
#define EXPLICIT 1
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);
double uni_rand(){
	return distribution(generator);
}

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
	double x, y, z;                  // position, also color (r,g,b)
	Vec(double x_ = 0, double y_ = 0, double z_ = 0){ x = x_; y = y_; z = z_; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x*b, y*b, z*b); }
	Vec mult(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec& norm(){ return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
	double dot(const Vec &b) const { return x*b.x + y*b.y + z*b.z; } // cross:
	Vec operator%(Vec&b){ return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
};

struct Ray {
	Vec o, d;
	Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Sphere {
	double rad;       // radius
	Vec p, e, c;      // position, emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
	double intersect(const Ray &r) const { // returns distance, 0 if nohit
		Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, eps = 1e-4, b = op.dot(r.d), det = b*b - op.dot(op) + rad*rad;
		if (det<0) return 0; else det = sqrt(det);
		return (t = b - det)>eps ? t : ((t = b + det) > eps ? t : 0);
	}
};

Sphere spheres[] = {//Scene: radius, position, emission, color, material
	Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),//Left
	Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF),//Rght
	Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),//Back
	Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),//Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Botm
	Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Top
	Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1)*.999, SPEC),//Mirr
	Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1)*.999, REFR),//Glas
	//Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF) //Lite
	Sphere(1.5, Vec(50, 81.6 - 16.5, 81.6), Vec(600, 600, 600), Vec(), DIFF),//Lite
};

inline double clamp(double x){ return x < 0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

inline bool intersect(const Ray &r, double &t, int &id){
	double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
	for (int i = int(n); i--;) if ((d = spheres[i].intersect(r)) && d < t){ t = d; id = i; }
	return t < inf;
}

Vec radiance(const Ray &r, int depth){
	double t;                               // distance to intersection
	int id = 0;                               // id of intersected object
	if (!intersect(r, t, id)) return Vec(); // if miss, return black
	const Sphere &obj = spheres[id];        // the hit object
	Vec x = r.o + r.d*t, n = (x - obj.p).norm(), nl = n.dot(r.d)<0 ? n : n*-1, f = obj.c;
	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
	if (depth > 100) return obj.e; // *** Added to prevent stack overflow
	if (++depth>5) if (uni_rand()<p) f = f*(1 / p); else return obj.e; //R.R.
	if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
		double r1 = 2 * M_PI*uni_rand(), r2 = uni_rand(), r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x)>.1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w%u;
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();
		return obj.e + f.mult(radiance(Ray(x, d), depth));
	}
	else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
		return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth));
	Ray reflRay(x, r.d - n * 2 * n.dot(r.d));     // Ideal dielectric REFRACTION
	bool into = n.dot(nl)>0;                // Ray from outside going in?
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
	if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0)    // Total internal reflection
		return obj.e + f.mult(radiance(reflRay, depth));
	Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
	double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
	return obj.e + f.mult(depth>2 ? (uni_rand()<P ?   // Russian roulette
		radiance(reflRay, depth)*RP : radiance(Ray(x, tdir), depth)*TP) :
		radiance(reflRay, depth)*Re + radiance(Ray(x, tdir), depth)*Tr);
}

Vec Radiance_ori(const Ray &r, int depth){
	//Vec x; //intersection point
	int obj_hit = -1;
	double t_min = 1e20;
	double t = 1e20;
	bool flag = false;
	//intersect the ray to the objects in the scene
	double n_elem = sizeof(spheres) / sizeof(Sphere);
	for (int i = 0; i < int(n_elem); i++){
		t = spheres[i].intersect(r);
		if (t < t_min && t != 0){
			obj_hit = i;
			t_min = t;
			flag = true;
		}
	}
	if (!flag)
		return Vec();
	Vec x = r.o + r.d*t_min;
	int id = obj_hit;
	const Sphere &obj = spheres[obj_hit];        // the hit object
	//x = r.o + r.d*t;
    Vec n = (x - obj.p).norm(), f = obj.c;
	bool into = n.dot(r.d) > 0 ? false : true;              // Ray from outside going in?
	Vec nl = into ? n : n*-1;

	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
	if (depth > 100) return obj.e; // *** Added to prevent stack overflow
	if (++depth>5) if (uni_rand()<p) f = f*(1 / p); else return obj.e; //R.R.
	if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
		double r1 = 2 * M_PI*uni_rand();
		double r2 = uni_rand();
		double r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w%u;
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();   //PBRT: P.775
		return spheres[id].e + f.mult(Radiance_ori(Ray(x, d), depth));
	}
	else if (obj.refl == SPEC){           // Ideal SPECULAR reflection
		//1.mirror, reflection only.
		Vec ray_dir = r.d - n * 2 * (r.d.dot(n));
		return spheres[id].e + f.mult(Radiance_ori(Ray(x, ray_dir), depth));
	}
	double nc = 1, nt = 1.5;
	double nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl);
	Ray reflection_ray(x, r.d - n * 2 * (r.d.dot(n)));
	// total internal reflection
	// cosine of reflection angle
	double cos_sq = 1 - nnt*nnt*(1 - ddn*ddn);
	if (cos_sq < 0)
		return spheres[id].e + f.mult(Radiance_ori(reflection_ray, depth));
	// reflection and refration
	// we already have the reflection ray.
	//refration
	// solve for the direction of refration ray.
	// tangent and normal.
	//Vec tangent = (r.d - n*ddn)*nnt;
	//Vec normal  = nl*-1 * sqrt(cos_sq);
	Vec tdir = (r.d*nnt - nl*(ddn*nnt + sqrt(cos_sq))).norm();

	//fresnel refraction
	double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1- (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0)*c*c*c*c*c;
	double Tr = 1 - Re; // transmittion?
	double P_Re = .25 + .5*Re;
	double P_Tr = 1 - P_Re;
	double RP = Re / P_Re;
	double TP = Tr / P_Tr;

	return spheres[id].e + f.mult(depth>2 ? (uni_rand() < P_Re ?   // Russian roulette
		Radiance_ori(reflection_ray, depth)*RP : Radiance_ori(Ray(x, tdir), depth)*TP) :
		Radiance_ori(reflection_ray, depth)*Re + Radiance_ori(Ray(x, tdir), depth)*Tr);
}

// input: ray, scene (objects), depth (number of bounces)
Vec Radiance(Ray r, int depth, int E = 1){

	Vec x; //intersection point
	int obj_hit = -1;
	double t_min = 1e20;
	double t = 1e20;
	bool flag = false;

	//intersect the ray to the objects in the scene
	double n_elem = sizeof(spheres) / sizeof(Sphere);
	for (int i = 0; i < int(n_elem); i++){
		t = spheres[i].intersect(r);
		if (t < t_min && t != 0){
			obj_hit = i;
			t_min = t;
			flag = true;
		}
	}
	if (!flag)
		return Vec();
	x = r.o + r.d*t_min;


	if (depth > 100) // too many bounces.
		return spheres[obj_hit].e;
	Vec f = spheres[obj_hit].c;
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
	if (++depth > 5) if (uni_rand() < p) f = f*(1 / p); else return spheres[obj_hit].e*E; //R.R.

	Vec n = (x - spheres[obj_hit].p).norm(); //normal of intersection point x. DO NOT FORGET TO NORMALIZE IT
	double nc = 1, nt = 1.5;
	//from outside?
	bool outside = n.dot(r.d) > 0 ? false : true;
	Vec nl = outside ? n : n*-1;
	//switch between different surface material.
	if (spheres[obj_hit].refl == SPEC){
		//1.mirror, reflection only.
		Vec ray_dir = r.d - n * 2 * (r.d.dot(n));
		return spheres[obj_hit].e + f.mult(Radiance(Ray(x, ray_dir), depth));
	}
	else if (spheres[obj_hit].refl == REFR){
		//2.glass, reflection and refration.
		double nnt = outside ? nc / nt : nt / nc, ddn = r.d.dot(nl);
		Ray reflection_ray(x, r.d - n * 2 * (r.d.dot(n)));
		// total internal reflection
		// cosine of reflection angle
		double cos_sq = 1 - nnt*nnt*(1 - ddn*ddn);
		if (cos_sq < 0)
			return spheres[obj_hit].e + f.mult(Radiance(reflection_ray, depth));
		// reflection and refration
		// we already have the reflection ray.
		//refration
		// solve for the direction of refration ray.
		// tangent and normal.
		//Vec tangent = (r.d - n*ddn)*nnt;
		//Vec normal  = nl*-1 * sqrt(cos_sq);
		Vec tdir = (r.d*nnt - nl*(ddn*nnt + sqrt(cos_sq))).norm();

		//fresnel refraction
		double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1- (outside ? -ddn : tdir.dot(n));
		double Re = R0 + (1 - R0)*c*c*c*c*c;
		double Tr = 1 - Re; // transmittion?
		double P_Re = .25 + .5*Re;
		double P_Tr = 1 - P_Re;
		double RP = Re / P_Re;
		double TP = Tr / P_Tr;

		return spheres[obj_hit].e + f.mult(depth>2 ? (uni_rand() < P_Re ?   // Russian roulette
			Radiance(reflection_ray, depth)*RP : Radiance(Ray(x, tdir), depth)*TP) :
			Radiance(reflection_ray, depth)*Re + Radiance(Ray(x, tdir), depth)*Tr);
	}
	else if (spheres[obj_hit].refl == DIFF){
		//3.matte, diffusion.
		double r1 = 2 * M_PI*uni_rand();
		double r2 = uni_rand();
		double r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w%u;
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();   //PBRT: P.775

#if EXPLICIT
		// Loop over any lights
		Vec e;
		for (int i = 0; i<n_elem; i++){
			const Sphere &s = spheres[i];
			if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0) continue; // skip non-lights
			int obj_hit_light = -1;
			Vec sw = s.p - x, su = ((fabs(sw.x)>.1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw%su;
			double cos_a_max = sqrt(1 - s.rad*s.rad / (x - s.p).dot(x - s.p));
			double eps1 = uni_rand(), eps2 = uni_rand();
			double cos_a = 1 - eps1 + eps1*cos_a_max;
			double sin_a = sqrt(1 - cos_a*cos_a);
			double phi = 2 * M_PI*eps2;
			Vec l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
			l.norm();
			if (intersect(Ray(x, l), t, obj_hit_light) && obj_hit_light == i){  // shadow ray
				double omega = 2 * M_PI*(1 - cos_a_max);
				e = e + f.mult(s.e*l.dot(nl)*omega)*M_1_PI;  // 1/pi for brdf
			}
		}
        	return spheres[obj_hit].e*E + e + f.mult(Radiance(Ray(x, d), depth, 0));
#else
		return spheres[obj_hit].e + f.mult(Radiance(Ray(x, d), depth));

#endif
	}


}

int main(int argc, char *argv[]){

	int defalt_sample = 1024;

	//int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : defalt_sample; // # samples
	int w = 100, h = 100, samps = argc == 2 ? atoi(argv[1]) / 4 : defalt_sample; // # samples

	Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
	Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, r, *c = new Vec[w*h];

#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
	for (int y = 0; y < h; y++){                       // Loop over image rows
		// *** Commented out for Visual Studio, fprintf is not thread-safe
		//fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100.*y / (h - 1));
		for (unsigned short x = 0; x < w; x++)   // Loop cols
			for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++)     // 2x2 subpixel rows
				for (int sx = 0; sx < 2; sx++, r = Vec()){        // 2x2 subpixel cols
					for (int s = 0; s < samps; s++){
						double r1 = 2 * uni_rand(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double r2 = 2 * uni_rand(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						Vec d = cx*(((sx + .5 + dx) / 2 + x) / w - .5) +
							cy*(((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
						r = r + Radiance(Ray(cam.o + d * 140, d.norm()), 0)*(1. / samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior
					c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
				}
	}
	FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i < w*h; i++)
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}