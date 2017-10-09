#pragma once
#include "coating.h"
#include "geometry.h"



struct Landing {
	double xl, yl, zl;
	unsigned int land_id = UINT_MAX;

	Landing() { land_id = UINT_MAX; }
	void reset() { land_id = UINT_MAX; }

};




void deposit_cpu(coating::particle* p_arg, Triangle* tri_arg, int MAXNUM_arg, Landing* l_arg) {

	for (int tid = 0; tid < MAXNUM_arg; tid++) {

		double nx = sin(p_arg->omega_p)*cos(p_arg->theta_p);
		double ny = sin(p_arg->omega_p)*sin(p_arg->theta_p);
		double nz = -cos(p_arg->omega_p);


		//if (tri_arg[tid].nx*nx + tri_arg[tid].ny*ny + tri_arg[tid].nz*nz < 0){
		//	continue; // the injected particle does not land on this plane
		//}
		// now find intersection of particle with the triangle plane

		// xp=x0+nx*t, yp=y0+ny*t, zp=z0+nz*t
		// a'x+b'y+c'z=a*'x'0+b'*y'0+c'*z'0

		double num = tri_arg[tid].nx*(tri_arg[tid].p0.x - p_arg->x_p) + tri_arg[tid].ny*(tri_arg[tid].p0.y - p_arg->y_p) + tri_arg[tid].nz*(tri_arg[tid].p0.z - p_arg->z_p);
		double denum = tri_arg[tid].nx*nx + tri_arg[tid].ny*ny + tri_arg[tid].nz*nz;

		double t = num / denum;

		double xpp = nx*t + p_arg->x_p;
		double ypp = ny*t + p_arg->y_p;
		double zpp = nz*t + p_arg->z_p;

		double x_p = coating::gun::z_g*tan(p_arg->omega_p)*cos(p_arg->theta_p) + coating::gun::x_g;
		double y_p = coating::gun::z_g*tan(p_arg->omega_p)*sin(p_arg->theta_p) + coating::gun::y_g;
		double z_p = -coating::gun::z_g;

		Point O(xpp, ypp, zpp);

		Triangle OP0P1(O, tri_arg[tid].p0, tri_arg[tid].p1);
		Triangle OP0P2(O, tri_arg[tid].p0, tri_arg[tid].p2);
		Triangle OP1P2(O, tri_arg[tid].p1, tri_arg[tid].p2);

		//cout << OP0P1.area() + OP0P2.area() + OP1P2.area() <<"   "<< tri_arg[tid].area() << endl;

		//if ((OP0P1.area() + OP0P2.area() + OP1P2.area() < tri_arg[tid].area() + 1e-6) && (OP0P1.area() + OP0P2.area() + OP1P2.area() > tri_arg[tid].area() - 1e-6)){
		if (abs((OP0P1.area() + OP0P2.area() + OP1P2.area()) - tri_arg[tid].area()) < 1e-5) {
			l_arg->xl = xpp;
			l_arg->yl = ypp;
			l_arg->zl = zpp;
			l_arg->land_id = tid;

			p_arg->impact();
			p_arg->dmax_p = p_arg->spreadfactor()*p_arg->D_p;
			p_arg->h_p = 2 * p_arg->D_p / (3.0*p_arg->spreadfactor()*p_arg->spreadfactor());
			tri_arg[tid].h = tri_arg[tid].h + (p_arg->h_p*(0.25*PI*p_arg->D_p*p_arg->D_p) / tri_arg[tid].area());


			//cout << "FOUND A MATCH!" << endl;
			//cout << xpp << "   " << ypp << "   " << zpp<<endl;

			//return;
		}






	}
}



	//struct Landing {
	//	double xl, yl, zl;
	//	unsigned int land_id = UINT_MAX;

	//	Landing() { land_id = UINT_MAX; }
	//	void reset() { land_id = UINT_MAX; }

	//};

bool isTriangle_inside_splat(coating::particle& p_arg,Triangle* t_arg) {
	// this function returns true if any points of a given triangle resides inside the splat diameter

	if ((t_arg->p0.x - p_arg.x_p)*(t_arg->p0.x - p_arg.x_p) +
		(t_arg->p0.y - p_arg.y_p)*(t_arg->p0.y - p_arg.y_p) +
		(t_arg->p0.z - p_arg.z_p)*(t_arg->p0.z - p_arg.z_p)
		< (p_arg.spreadfactor()*p_arg.D_p*p_arg.spreadfactor()*p_arg.D_p)) return true;

	if ((t_arg->p1.x - p_arg.x_p)*(t_arg->p1.x - p_arg.x_p) +
		(t_arg->p1.y - p_arg.y_p)*(t_arg->p1.y - p_arg.y_p) +
		(t_arg->p1.z - p_arg.z_p)*(t_arg->p1.z - p_arg.z_p)
		< (p_arg.spreadfactor()*p_arg.D_p*p_arg.spreadfactor()*p_arg.D_p)) return true;

	if ((t_arg->p2.x - p_arg.x_p)*(t_arg->p2.x - p_arg.x_p) +
		(t_arg->p2.y - p_arg.y_p)*(t_arg->p2.y - p_arg.y_p) +
		(t_arg->p2.z - p_arg.z_p)*(t_arg->p2.z - p_arg.z_p)
		< (p_arg.spreadfactor()*p_arg.D_p*p_arg.spreadfactor()*p_arg.D_p)) return true;

	return false;
}



void deposit_particle(coating::particle& p_arg, vector<Triangle>& tri_arg, unsigned int triId) {








}

