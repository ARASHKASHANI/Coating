#include "geometry.h"
#include "coating.h"


using namespace std;
using std::default_random_engine;
using std::normal_distribution;


#define CHECK_GPU(call) \
{ \
const cudaError_t error = call; \
	if (error != cudaSuccess) \
						{ \
		printf("\n\n");\
		printf("Error occured @ \n");\
		printf("FILE= %s   :   LINE= %d, \n", __FILE__, __LINE__); \
		printf("reason: %s\n\n", cudaGetErrorString(error)); \
		printf("\n");\
		getchar();\
		exit(1);\
					}\
}\


#define CHECK_GPU_KERNEL(error) \
{ \
	if (error != cudaSuccess) \
								{ \
		printf("\n\n");\
		printf("Error occured @ \n");\
		printf("FILE= %s   :   LINE= %d, \n", __FILE__, __LINE__); \
		printf("reason: %s\n\n",  cudaGetErrorString(error)); \
		printf("\n");\
		getchar();\
		exit(1);\
						}\
}\

cudaError_t error = cudaSuccess;

struct Landing {
	double xl, yl, zl;
	unsigned int land_id = UINT_MAX;

	/*__host__ __device__*/  Landing() { land_id = UINT_MAX; }
	/*__host__ __device__*/ void reset() { land_id = UINT_MAX; }

};

struct particleprop {
	double rho_p;
	double cp_p;
	double mu_p; // viscosity of particle
	double k_p; // heat conduction coefficient of particle
	double sigma_p; // surface tension coefficient of particle
	double Hf_p; // latent heat of partcile
	double cangle_p; // liquid-solid contact angle


};


__global__ void deposit(coating::particle* p_arg, Triangle* tri_arg, Point* P_arg, int MAXNUM_arg, double Tsub_arg, particleprop parprop_arg, Landing* l_arg) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid < MAXNUM_arg) {


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

		/*double x_p = coating::gun::z_g*tan(p_arg->omega_p)*cos(p_arg->theta_p) + coating::gun::x_g;
		double y_p = coating::gun::z_g*tan(p_arg->omega_p)*sin(p_arg->theta_p) + coating::gun::y_g;
		double z_p = -coating::gun::z_g;*/

		Point O(xpp, ypp, zpp);

		Triangle OP0P1(O, tri_arg[tid].p0, tri_arg[tid].p1);
		Triangle OP0P2(O, tri_arg[tid].p0, tri_arg[tid].p2);
		Triangle OP1P2(O, tri_arg[tid].p1, tri_arg[tid].p2);

		//cout << OP0P1.area() + OP0P2.area() + OP1P2.area() <<"   "<< tri_arg[tid].area() << endl;

		//if ((OP0P1.area() + OP0P2.area() + OP1P2.area() < tri_arg[tid].area() + 1e-6) && (OP0P1.area() + OP0P2.area() + OP1P2.area() > tri_arg[tid].area() - 1e-6)){
		
		if (tri_arg[tid].nx*nx + tri_arg[tid].ny*ny + tri_arg[tid].nz*nz < 0) {
			//if (abs((OP0P1.area() + OP0P2.area() + OP1P2.area()) - tri_arg[tid].area()) < 1e-5) {
			if (!( (OP0P1.area() + OP0P2.area() + OP1P2.area()) > tri_arg[tid].area() ) ) {
				l_arg->xl = xpp;
				l_arg->yl = ypp;
				l_arg->zl = zpp;
				l_arg->land_id = tid;

				//cout << "FOUND A MATCH!" << endl;
				//cout << xpp << "   " << ypp << "   " << zpp<<endl;


				p_arg->V_z = p_arg->V_p*cos(p_arg->omega_p);
				p_arg->V_x = p_arg->V_p*sin(p_arg->omega_p)*cos(p_arg->theta_p);
				p_arg->V_y = p_arg->V_p*sin(p_arg->omega_p)*sin(p_arg->theta_p);

				double We = parprop_arg.rho_p*p_arg->V_z*p_arg->V_z*p_arg->D_p / parprop_arg.sigma_p;
				double Re = parprop_arg.rho_p*p_arg->V_z*p_arg->D_p / parprop_arg.mu_p;
				double Pe = p_arg->V_z* p_arg->D_p *  parprop_arg.rho_p* parprop_arg.cp_p / parprop_arg.k_p;
				double St = parprop_arg.cp_p*(p_arg->T_p - Tsub_arg) / parprop_arg.Hf_p;

				double num = We + 12.0;
				double denum = 3.0*(1.0 - cos(parprop_arg.cangle_p)) + (4.0*We / sqrt(Re)) + We*sqrt(0.75*St / Pe);

				double spreadfactor = sqrt(num / denum);

				p_arg->dmax_p = spreadfactor*p_arg->D_p;
				p_arg->h_p = 2 * p_arg->D_p / (3.0*spreadfactor*spreadfactor);

				double PI = 3.14159265;

				tri_arg[tid].h = tri_arg[tid].h + (p_arg->h_p*(0.25*PI*p_arg->D_p*p_arg->D_p) / tri_arg[tid].area());

				P_arg[tri_arg[tid].pid[0]].h = tri_arg[tid].h;
				P_arg[tri_arg[tid].pid[1]].h = tri_arg[tid].h;
				P_arg[tri_arg[tid].pid[2]].h = tri_arg[tid].h;

				return;
			}

		}


	}

}



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
		if (abs((OP0P1.area() + OP0P2.area() + OP1P2.area()) - tri_arg[tid].area())<1e-5) {
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



void write_file_triangles(Triangle* t_arg, Point* p_arg, unsigned int gi_arg) {
	ofstream f;
	stringstream filename;
	filename << "";
	filename << "output--" << gi_arg << ".dat";
	//filename << "output.dat";

	f.open(filename.str());

	f << "TITLE = \"Thickness plot\"" << endl;
	f << "VARIABLES = \"X\", \"Y\", \"Z\", \"h\" \"Z'\" " << endl;
	f << "ZONE N = " << NUM_POINTS << ", E = " << NUM_TRIANGLES << ", F = FEPOINT,  ET=TRIANGLE " << endl;

	/*for (unsigned int j = 0; j < NY; j++) {
	for (unsigned int i = 0; i < NX; i++) {

	unsigned int ij = i + NX*j;*/
	for (unsigned int i = 0; i < NUM_POINTS; i++) {

		f << p_arg[i].x << "   " << p_arg[i].y << "   " << p_arg[i].z << "   " << p_arg[i].h << "   " << p_arg[i].z+ p_arg[i].h << endl;
		//p_arg[t_arg[i].pid[0]].x  << "   " << p_arg[t_arg[i].pid[0]].y  << "   " << p_arg[t_arg[i].pid[0]].z  << "   " << t_arg[i].h << endl;
		//f << p_arg[t_arg[i].pid[1]].x  << "   " << p_arg[t_arg[i].pid[1]].y  << "   " << p_arg[t_arg[i].pid[1]].z << "   "  << t_arg[i].h << endl;
		//f << p_arg[t_arg[i].pid[2]].x  << "   " << p_arg[t_arg[i].pid[2]].y  << "   " << p_arg[t_arg[i].pid[2]].z << "   "  << t_arg[i].h << endl;
	}

	for (unsigned int i = 0; i < NUM_TRIANGLES; i++) {
		f << t_arg[i].pid[0] + 1 << "   " << t_arg[i].pid[1] + 1 << "   " << t_arg[i].pid[2] + 1 << endl;

	}

	f.flush();
	f.close();

}


using std::cout;
using namespace coating;



int main() {


	pt_data mydata;

	
	find_avglength("Cyl.stl");
	/*AVG_L = 0.0159931;
	NUM_TRIANGLES = 198899;
	NUM_POINTS_RAW = 596697;*/

	mydata = from_stl("Cyl.stl");
	Triangle* t = mydata.t;
	Point* P = mydata.p;

	find_neighbors(t, P, NUM_TRIANGLES);
	triangle_report(t, P, NUM_TRIANGLES);

	

	write_file_triangles(t, P, 50000);
	//cout << endl << "Press Enter ...";
	//getchar();

	//Triangle* t;
	//ifstream f("triangles.txt");
	//string s;
	//stringstream ss;


	//NUM_TRIANGLES = 198899;
	//f >> NUM_TRIANGLES;
	//f >> s;

	//getline(f, s);
	//ss << s;
	//ss >> NUM_TRIANGLES;
	//ss.clear();


	//t = new Triangle[NUM_TRIANGLES];
	//for (int i = 0; i < NUM_TRIANGLES; i++) {

	//	//f >> s;



	//	getline(f, s);

	//	ss << s;
	//	ss >> i >> t[i].nx >> t[i].ny >> t[i].nz >>
	//		t[i].p0.x >> t[i].p0.y >> t[i].p0.z >>
	//		t[i].p1.x >> t[i].p1.y >> t[i].p1.z >>
	//		t[i].p2.x >> t[i].p2.y >> t[i].p2.z >>
	//		t[i].neighbor[0] >> t[i].neighbor[1] >> t[i].neighbor[2];
	//	t[i].h = 0.0;
	//	ss.clear();
	//	//cout << "\r" << i;
	//}
	//f.close();
	//cout << endl;



	Triangle* t_dev;
	CHECK_GPU(cudaMalloc((void **)&t_dev, NUM_TRIANGLES * sizeof(Triangle)));
	CHECK_GPU(cudaMemcpy(t_dev, t, NUM_TRIANGLES * sizeof(Triangle), cudaMemcpyHostToDevice));


	Point* P_dev;
	CHECK_GPU(cudaMalloc((void **)&P_dev, NUM_POINTS * sizeof(Point)));
	CHECK_GPU(cudaMemcpy(P_dev, P, NUM_POINTS * sizeof(Point), cudaMemcpyHostToDevice));

	coating::particle* p_dev;
	CHECK_GPU(cudaMalloc((void **)&p_dev, 1 * sizeof(coating::particle)));

	Landing l;
	l.xl = 0.0;
	l.yl = 0.0;
	l.zl = 0.0;


	Landing* l_dev;
	CHECK_GPU(cudaMalloc((void **)&l_dev, 1 * sizeof(Landing)));
	CHECK_GPU(cudaMemcpy(l_dev, &l, 1 * sizeof(Landing), cudaMemcpyHostToDevice));



	// settings parameters
	coating::settings::Dmean = 50.0;
	coating::settings::Dstd = 10.0;

	coating::settings::Vmean = 60.0;
	coating::settings::Vstd = 5.7;

	coating::settings::Tmean = 1609.0;
	coating::settings::Tstd = 219.0;

	coating::settings::Omean = 0.0;
	coating::settings::Ostd = 5.0;



	// substrate parameters

	coating::Substrate::XMIN = 0.0;
	coating::Substrate::XMAX = 0.005;
	coating::Substrate::YMIN = 0.0;
	coating::Substrate::YMAX = 0.005;
	NX = 250;
	NY = 250;
	coating::Substrate::DELX = (coating::Substrate::XMAX - coating::Substrate::XMIN) / NX;
	coating::Substrate::DELY = (coating::Substrate::XMAX - coating::Substrate::XMIN) / NY;
	coating::Substrate::T_sub = 200.0;


	//coating::Substrate::h = new double[NX*NY];
	//h = new double[NX*NY];
	//
	//	for (unsigned int i = 0; i < NX*NY; i++) { /*coating::Substrate::h[i] = 0.0; */ h[i] = 0.0; }



	// gun parameters

	coating::gun::x_g = 0.0;
	coating::gun::y_g = 0.0;
	coating::gun::z_g = 10.0;

	coating::gun::u_g = 0.0;
	coating::gun::v_g = 0.0;
	coating::gun::w_g = 0.0;

	coating::gun::mfrate = 0.126*0.001;
	coating::gun::injtime = 0.005;


	// particle properties

	coating::particle::rho_p = 3000.0;;
	coating::particle::cp_p = 1300.0;;
	coating::particle::mu_p = 0.175E-4;; // viscosity of particle
	coating::particle::k_p = 60.0;; // heat conduction coefficient of particle
	coating::particle::sigma_p = 0.69;; // surface tension coefficient of particle
	coating::particle::Hf_p = 0.1075E7;; // latent heat of partcile
	coating::particle::cangle_p = 0.0;; // liquid-solid contact angle

	particleprop parprop;
	parprop.rho_p = particle::rho_p;
	parprop.cp_p = particle::cp_p;
	parprop.mu_p = particle::mu_p;
	parprop.k_p = particle::k_p;
	parprop.sigma_p = particle::sigma_p;
	parprop.Hf_p = particle::Hf_p;
	parprop.cangle_p = particle::cangle_p;


	// --------------------------------------------------------------------------------------------------------------- //



	unsigned int numtime = 1000;
	double dt = coating::gun::injtime / numtime;
	double time = 0.0;
	double mass_total = coating::gun::injtime*coating::gun::mfrate;
	coating::particle p;
	p.reset();
	CHECK_GPU(cudaMemcpy(p_dev, &p, 1 * sizeof(coating::particle), cudaMemcpyHostToDevice));


	//deposit_cpu(&p, t, NUM_TRIANGLES, &l);
	deposit << < int(NUM_TRIANGLES / 1024) + 1, 1024 >> > (p_dev, t_dev, P_dev, NUM_TRIANGLES, Substrate::T_sub, parprop, l_dev);
	error = cudaGetLastError();
	CHECK_GPU_KERNEL(error);



	CHECK_GPU(cudaMemcpy(&l, l_dev, 1 * sizeof(Landing), cudaMemcpyDeviceToHost));


	//cout << t[l.land_id].area() << endl;


	double PI = 3.14159265;

	gi = 0;

	//while (time<coating::gun::injtime) {


	while (gi<numtime) {

		double mass_dt_init = coating::gun::mfrate*dt;
		double mass_dt = coating::gun::mfrate*dt;
		double mass_dt_old = mass_dt;

		gii = 0;
		tStart = cpuSecond();

		while (mass_dt >(PI / 6.0)*(1E-6)*(1E-6)*(1.E-6)) { // while remaining mass is larger than mass of 1 micron particle

			p.reset();
			CHECK_GPU(cudaMemcpy(p_dev, &p, 1 * sizeof(coating::particle), cudaMemcpyHostToDevice));


			//deposit_cpu(&p, t, NUM_TRIANGLES, &l);
			deposit << < int(NUM_TRIANGLES / 1024) + 1, 1024 >> > (p_dev, t_dev, P_dev, NUM_TRIANGLES, Substrate::T_sub, parprop, l_dev);
			error = cudaGetLastError();
			CHECK_GPU_KERNEL(error);

			CHECK_GPU(cudaMemcpy(&l, l_dev, 1 * sizeof(Landing), cudaMemcpyDeviceToHost));


			mass_dt -= (PI / 6.0)*p.D_p*p.D_p*p.D_p;


			gii++;
			//cout << ii << "   " << mass_dt<< endl;

			//cout << p.D_p << "   " << p.T_p << "   " << p.theta_p << std::endl;

			//_getch();
			if (!(gii % 10000)) {
				//cout << "*-----------------------------------*" << endl;
				cout << "\rparticle#: " << gii << "   mass/mass0%= " << (mass_dt_old / mass_dt_init)*100.0;
				//cout << "*-----------------------------------*" << endl;
			}

			if (mass_dt> (PI / 6.0)*(1E-6)*(1E-6)*(1.E-6)) {
				mass_dt_old = mass_dt;
			}


		}


		cout << "\rparticle#: " << gii << "   mass/mass0%= " << (mass_dt_old / mass_dt_init)*100.0;
		gi++;
		cout << endl << "Ntime= " << gi << " time= " << time << "  / endtime= " << coating::gun::injtime << endl;
		tStop = cpuSecond();
		cout << "elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
		cout << "*-------------------------------------------*" << endl;


		//write_file(coating::Substrate::h, NX, NY);
		//write_file(h, NX, NY);
		if (!(gi % 100) || gi == 1) {
			//write_file(NX, NY);


			CHECK_GPU(cudaMemcpy(t, t_dev, NUM_TRIANGLES * sizeof(Triangle), cudaMemcpyDeviceToHost));
			CHECK_GPU(cudaMemcpy(P, P_dev, NUM_POINTS * sizeof(Point), cudaMemcpyDeviceToHost));
			write_file_triangles(t, P, gi);
		}

		//cout << "Press a key ..." << endl;
		//_getch();
		//Sleep(10 * 1000);

		time += dt;
	}

	//if (gi==numtime) write_file(NX, NY);
	//write_file(NX, NY);

	if (gi == numtime) {
		CHECK_GPU(cudaMemcpy(t, t_dev, NUM_TRIANGLES * sizeof(Triangle), cudaMemcpyDeviceToHost));
		CHECK_GPU(cudaMemcpy(P, P_dev, NUM_POINTS * sizeof(Point), cudaMemcpyDeviceToHost));
		write_file_triangles(t, P, gi);
	}
	CHECK_GPU(cudaMemcpy(t, t_dev, NUM_TRIANGLES * sizeof(Triangle), cudaMemcpyDeviceToHost));
	CHECK_GPU(cudaMemcpy(P, P_dev, NUM_POINTS * sizeof(Point), cudaMemcpyDeviceToHost));
	write_file_triangles(t, P, gi);

	//cout << gi << endl;
	//_getch();

	//delete[] coating::Substrate::h;
	delete[] h;


	cudaFree(t_dev);
	cudaFree(p_dev);
	cudaFree(l_dev);


	_getch();

	return 0;
}






