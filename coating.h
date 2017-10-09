#pragma once

#include<time.h>

#include <random>
#include <conio.h>

#include <string>
#include <sstream>
#include <fstream>



double* h = new double[250 * 250];
//double h[1000 * 1000];

unsigned int NX = 250;
unsigned int NY = 250;


unsigned int gi = 0;
unsigned int gii = 0;


// Constants
double PI = 3.14159265;

double gx = 0.0;
double gy = 0.0;
double gz = 9.81;

std::default_random_engine generator((unsigned int)time(0));


namespace coating {

	//std::default_random_engine generator((unsigned int)time(0));


	class settings {

	public:

		static double delt;
		static double mass_t; // mass injected per timestep

		// mean and standard variation of variables
		static double Vmean, Vstd;
		static double Dmean, Dstd;
		static double Tmean, Tstd;
		static double Omean, Ostd; // dispersion angle


	};

	double settings::Vmean = 0.0;
	double settings::Vstd = 0.0;
	double settings::Dmean = 0.0;
	double settings::Dstd = 0.0;
	double settings::Tmean = 0.0;
	double settings::Tstd = 0.0;
	double settings::Omean = 0.0;
	double settings::Ostd = 0.0;

	//-------------------------------------------------------------------------------------------------------------------------//

	class gun {

	public:

		static double x_g;
		static double y_g;
		static double z_g;

		static double u_g;
		static double v_g;
		static double w_g;

		static double mfrate; // mass flow rate
		static double injtime; // total injection time;


	};

	double gun::x_g = 0.0;
	double gun::y_g = 0.0;
	double gun::z_g = 0.0;

	double gun::u_g = 0.0;
	double gun::v_g = 0.0;
	double gun::w_g = 0.0;

	double gun::mfrate = 0.0;
	double gun::injtime = 0.0;


	//-------------------------------------------------------------------------------------------------------------------------//


	class Substrate {

	public:

		//static double* h;

		// Domain/Substrate parameters
		static double XMIN, XMAX;
		static double YMIN, YMAX;
		static double DELX, DELY;
		//static unsigned int NX, NY;

		static double T_sub;

		Substrate() {

		}

		~Substrate() {

		}


	};

	//double* Substrate::h = NULL;

	double Substrate::XMIN = 0.0;
	double Substrate::XMAX = 0.0;
	double Substrate::YMIN = 0.0;
	double Substrate::YMAX = 0.0;
	double Substrate::DELX = 0.0;
	double Substrate::DELY = 0.0;
	//unsigned int Substrate::NX = 0;
	//unsigned int Substrate::NY = 0;

	double Substrate::T_sub = 0.0;


	//-------------------------------------------------------------------------------------------------------------------------//

	class particle {
	public:

		// particle parameters
		double D_p;
		double V_p; // particle velocity magnitude
		double V_x, V_y, V_z;
		double V_normal; // normal impact velocity
		double T_p;
		double x_p, y_p, z_p;
		double dmax_p; // diameter of cylinder
		double h_p; // height of cylinder

		double omega_p; // particle disperse angle
		double theta_p;  // particle azimuthal angle

		// bounding box parameters
		unsigned int ib_min;
		unsigned int ib_max;
		unsigned int jb_min;
		unsigned int jb_max;


		// particle properties
		static double rho_p;
		static double cp_p;
		static double mu_p; // viscosity of particle
		static double k_p; // heat conduction coefficient of particle
		static double sigma_p; // surface tension coefficient of particle
		static double Hf_p; // latent heat of partcile
		static double cangle_p; // liquid-solid contact angle

		//inline double Re() { return rho_p*abs(V_z)*D_p / mu_p; }
		inline double Re() { return rho_p*abs(V_normal)*D_p / mu_p; }
		//inline double We() { return rho_p*V_z*V_z*D_p / sigma_p; }
		inline double We() { return rho_p*V_normal*V_normal*D_p / sigma_p; }
		//inline double Pe() { return abs(V_z)*D_p * rho_p*cp_p / k_p; }
		inline double Pe() { return abs(V_normal)*D_p * rho_p*cp_p / k_p; }
		inline double St() { return cp_p*(T_p - Substrate::T_sub) / Hf_p; }

		inline double spreadfactor() {

			double num = We() + 12.0;
			double denum = 3.0*(1.0 - cos(cangle_p)) + (4.0*We() / sqrt(Re())) + We()*sqrt(0.75*St() / Pe());
			return sqrt(num / denum);
		}

		inline void find_boundingbox() {

			ib_min = int((x_p - 0.5*dmax_p) / Substrate::DELX);
			ib_max = int((x_p + 0.5*dmax_p) / Substrate::DELX);

			jb_min = int((y_p - 0.5*dmax_p) / Substrate::DELY);
			jb_max = int((y_p + 0.5*dmax_p) / Substrate::DELY);


		}

		inline void impact() {

			V_z = -V_p*cos(omega_p);
			V_x = V_p*sin(omega_p)*cos(theta_p);
			V_y = V_p*sin(omega_p)*sin(theta_p);

			//x_p = gun::z_g*tan(omega_p)*cos(theta_p) + gun::x_g;
			//y_p = gun::z_g*tan(omega_p)*sin(theta_p) + gun::y_g;



		}

		particle() {
			//reset();
			x_p = 0.0;
			y_p = 0.0;
			z_p = gun::z_g;


			std::normal_distribution<double> normdist_Dp(settings::Dmean, settings::Dstd);
			std::normal_distribution<double> normdist_Vp(settings::Vmean, settings::Vstd);
			std::normal_distribution<double> normdist_Tp(settings::Tmean, settings::Tstd);
			std::normal_distribution<double> normdist_Op(settings::Omean, settings::Ostd);


			std::uniform_real_distribution<double> unif(0.0, 360.0);

			D_p = normdist_Dp(generator);
			while (D_p < 1.0) { D_p = normdist_Dp(generator); }
			//D_p = normdist_Dp(generator)*1e-6;
			D_p = D_p*1e-6;

			V_p = normdist_Vp(generator);

			T_p = normdist_Tp(generator);
			while (T_p  < 0.0) { T_p = normdist_Tp(generator); }
			//T_p = normdist_Tp(generator);

			omega_p = normdist_Op(generator)*(PI / 180.0);
			theta_p = unif(generator)*(PI / 180.0);


			V_z = -V_p*cos(omega_p);
			V_x = V_p*sin(omega_p)*cos(theta_p);
			V_y = V_p*sin(omega_p)*sin(theta_p);
		}


		inline void deposit() {

			//if (x_p< coating::Substrate::XMIN || x_p> coating::Substrate::XMAX || y_p< coating::Substrate::YMIN || y_p> coating::Substrate::YMAX) return;
			if ((ib_max> NX - 1) || (jb_max> NY - 1)) return;

			for (unsigned int j = jb_min; j <= jb_max; j++) {
				for (unsigned int i = ib_min; i <= ib_max; i++) {

					unsigned int ij = i + NX*j;
					//coating::Substrate::h[ij] = coating::Substrate::h[ij] + (h_p*(0.25*PI*D_p*D_p) / ((ib_max - ib_min + 1)*(jb_max - jb_min + 1)*coating::Substrate::DELX*coating::Substrate::DELY));
					h[ij] = h[ij] + (h_p*(0.25*PI*D_p*D_p) / ((ib_max - ib_min + 1)*(jb_max - jb_min + 1)*coating::Substrate::DELX*coating::Substrate::DELY));

				}
			}

		}


		inline void reset() {

			x_p = 0.0;
			y_p = 0.0;
			z_p = gun::z_g;


			std::normal_distribution<double> normdist_Dp(settings::Dmean, settings::Dstd);
			std::normal_distribution<double> normdist_Vp(settings::Vmean, settings::Vstd);
			std::normal_distribution<double> normdist_Tp(settings::Tmean, settings::Tstd);
			std::normal_distribution<double> normdist_Op(settings::Omean, settings::Ostd);


			std::uniform_real_distribution<double> unif(0.0, 360.0);

			D_p = normdist_Dp(generator);
			while (D_p < 1.0) { D_p = normdist_Dp(generator); }
			//D_p = normdist_Dp(generator)*1e-6;
			D_p = D_p*1e-6;

			V_p = normdist_Vp(generator);

			T_p = normdist_Tp(generator);
			while (T_p  < 0.0) { T_p = normdist_Tp(generator); }
			//T_p = normdist_Tp(generator);

			omega_p = normdist_Op(generator)*(PI / 180.0);
			theta_p = unif(generator)*(PI / 180.0);

			V_z = -V_p*cos(omega_p);
			V_x = V_p*sin(omega_p)*cos(theta_p);
			V_y = V_p*sin(omega_p)*sin(theta_p);

			//impact();

			//spreadfactor();
			//dmax_p = spreadfactor()*D_p;
			//h_p = 2 * D_p / (3.0*spreadfactor()*spreadfactor());

			//find_boundingbox();

			//deposit();






		}





	};


	double particle::rho_p = 0.0;;
	double particle::cp_p = 0.0;;
	double particle::mu_p = 0.0;; // viscosity of particle
	double particle::k_p = 0.0;; // heat conduction coefficient of particle
	double particle::sigma_p = 0.0;; // surface tension coefficient of particle
	double particle::Hf_p = 0.0;; // latent heat of partcile
	double particle::cangle_p = 0.0;; // liquid-solid contact angle






}


