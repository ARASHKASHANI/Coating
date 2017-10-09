#include "geometry.h"
#include "coating.h"
#include "deposit.h"
//#include "kdtree.h"
#include "octree.h"

using namespace std;
using std::default_random_engine;
using std::normal_distribution;



bool isTri_inVector(unsigned int id, vector<unsigned int> &v){

	if (v.empty()) return false;
	else {

		if (std::find(v.begin(), v.end(), id) != v.end()) {
			return true;

		}

		return false;

	}
}



struct imp{

public:
	unsigned int id;
	Point im;

	imp(unsigned int i, Point p) : id(i), im(p) {}
	imp() { id = UINT_MAX; im.x = 0.0; im.y = 0.0; im.z = 0.0; }
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




using std::cout;
using namespace coating;

pt_data DATA;


double x, y, z;
double t;


imp brute_find(Ray& r, pt_data& D) {

	vector<pair<double, unsigned int>> CandidTri;

	Point pimpact;

	for (unsigned int i = 0; i < D.t.size(); i++) {
		double nx = r.n.vx;
		double ny = r.n.vy;
		double nz = r.n.vz;

		//double num = D.t[i].nx*(D.t[i].p0.x - r.p.x) +
		//	D.t[i].ny*(D.t[i].p0.y - r.p.y) + D.t[i].nz*(D.t[i].p0.z - r.p.z);
		//double denum = D.t[i].nx*nx + D.t[i].ny*ny + D.t[i].nz*nz;

		double p0x = D.p[D.t[i].pid[0]].first.x;
		double p0y = D.p[D.t[i].pid[0]].first.y;
		double p0z = D.p[D.t[i].pid[0]].first.z;

		double num = D.t[i].nx*(p0x - r.p.x) +
			D.t[i].ny*(p0y - r.p.y) + D.t[i].nz*(p0z - r.p.z);
		double denum = D.t[i].nx*nx + D.t[i].ny*ny + D.t[i].nz*nz;

		


		

		double t = 0;
		if (denum != 0) {
			t = num / denum;
		}
		// if (denum==0.0) then what



		double xpp = nx*t + r.p.x;
		double ypp = ny*t + r.p.y;
		double zpp = nz*t + r.p.z;

		pimpact.x = xpp;
		pimpact.y = ypp;
		pimpact.z = zpp;


		Point O(xpp, ypp, zpp);

		// now check if this point is inside the triangle

		Triangle OP0P1(O, D.t[i].p0, D.t[i].p1);
		Triangle OP0P2(O, D.t[i].p0, D.t[i].p2);
		Triangle OP1P2(O, D.t[i].p1, D.t[i].p2);

		//double d;
		//if (d = abs((OP0P1.area2() + OP0P2.area2() + OP1P2.area2()) - D.t[i].area2()) < 1e-5) {
			//distance_pair.push_back(pair<double, unsigned int>(d, numTri));
		//d = ((OP0P1.area() + OP0P2.area() + OP1P2.area()) - D.t[i].area());

		Vector v0(D.t[i].p0);
		Vector v1(D.t[i].p1);
		Vector v2(D.t[i].p2);

		Vector vp(xpp, ypp, zpp);

		Vector e0 = v2 - v1;
		Vector e1 = v0 - v2;
		Vector e2 = v1 - v0;

		Vector d0 = vp - v0;
		Vector d1 = vp - v1;
		Vector d2 = vp - v2;

		Vector np(D.t[i].nx, D.t[i].ny, D.t[i].nz);

		double AreaTot = (e0.crossProduct(e1)).dotProduct(np)*0.5;
		double Area0 = (e0.crossProduct(d2)).dotProduct(np)*0.5;
		double Area1 = (e1.crossProduct(d0)).dotProduct(np)*0.5;
		double Area2 = (e2.crossProduct(d1)).dotProduct(np)*0.5;

		double b0 = Area0 / AreaTot;
		double b1 = Area1 / AreaTot;
		double b2 = Area2 / AreaTot;



		//double alpha, beta, gamma;
		//alpha = OP0P1.area() / D.t[i].area();
		//beta = OP0P2.area() / D.t[i].area();
		//gamma = OP1P2.area() / D.t[i].area();
		////gamma = 1.0 - alpha - beta;


		//cout << "***** **** : id= " << i << "   " << b0 << "   " <<  b1 << "   " << b2  << endl;

		//if ( (alpha>=0.0 && alpha<=1.0) && (beta >= 0.0 && beta <= 1.0) && (gamma >= 0.0 && gamma <= 1.0) ){
		if ((b0 > 0.0 && b0 < 1.0) && (b1 > 0.0 && b1 < 1.0) && (b2 > 0.0 && b2 < 1.0)) {

		double distance = sqrt((xpp - r.p.x)*(xpp - r.p.x) + (ypp - r.p.y)*(ypp - r.p.y) + (zpp - r.p.z)*(zpp - r.p.z));

		CandidTri.push_back(pair<double, unsigned int>(distance, i));

		//cout << "found one! : id= " << i << "   " << d << endl;
		//cout << "found one! : id= " << i << "   " << b0+b1+b2-1.0 << endl;
		}

	}

		double minDistance = DBL_MAX;
		unsigned int returnId = UINT_MAX;

		for (unsigned int i = 0; i < CandidTri.size(); i++) {
			if (minDistance > CandidTri[i].first) {
				minDistance = CandidTri[i].first;
				returnId = CandidTri[i].second;
			}
		}

		return imp(returnId,pimpact);
	

	
}


void save(char* fname, pt_data &D){

	ofstream file(fname);

	file << "***** NUM POINTS & NUM_TRIANGLES *****" << endl;
	file << NUM_POINTS << endl;
	file << NUM_TRIANGLES << endl;

	file << "***** POINTS *****" << endl;
	for (unsigned int i = 0; i < NUM_POINTS; i++){
		file << D.p[i].first.x << "   " << D.p[i].first.y << "   " << D.p[i].first.z << "   " << D.p[i].first.h << endl;
	}

	file << "***** TRIANGLES *****" << endl;
	for (unsigned int i = 0; i < NUM_TRIANGLES; i++){
		
		file <<  D.t[i].nx << "   " << D.t[i].ny << "   " << D.t[i].nz<<endl;
		//file << D.t[i].p0.x << "   " << D.t[i].p0.y << "   " << D.t[i].p0.z << endl;
		//file << D.t[i].p1.x << "   " << D.t[i].p1.y << "   " << D.t[i].p1.z << endl;
		//file << D.t[i].p2.x << "   " << D.t[i].p2.y << "   " << D.t[i].p2.z << endl;
		file << D.t[i].h << endl;
		file << D.t[i].pid[0] << "   " << D.t[i].pid[1] << "   " << D.t[i].pid[2] << endl;
		file << i<< "   " << D.t[i].neighbor[0] << "   " << D.t[i].neighbor[1] << "   " << D.t[i].neighbor[2] << endl;
		
	}

	file << "***** LIST OF TRINAGLES AT EACH POINT  *****" << endl;
	for (unsigned int i = 0; i < NUM_POINTS; i++) {
	
		for (unsigned int ii : D.v[i]) {
			file << ii << "   ";
		}
		file <<  endl;
	}

	file.close();
	




}

void load(char* fname, pt_data &D){

	ifstream file(fname);

	stringstream ss;
	string s;


	double xmin = 0.0;
	double xmax = 0.0;
	double ymin = 0.0;
	double ymax = 0.0;
	double zmin = 0.0;
	double zmax = 0.0;

	unsigned int temp;
	
	getline(file, s); // line "***** NUM POINTS & NUM_TRIANGLES *****"
	s.clear();

	ss.clear();
	getline(file, s);
	ss << s;
	ss >> NUM_POINTS;

	ss.clear();
	getline(file, s);
	ss << s;
	ss >> NUM_TRIANGLES;


	getline(file, s);  // line "***** POINTS *****"
	s.clear();

	D.p.clear();
	D.p.reserve(NUM_POINTS);
	for (unsigned int i = 0; i < NUM_POINTS; i++){
		
		ss.clear();
		getline(file, s);
		ss << s;
		double x, y, z, h;
		ss >> x >> y >> z>>h;

		if (xmin > x) xmin = x;
		if (xmax < x) xmax = x;
		if (ymin > y) ymin = y;
		if (ymax < y) ymax = y;
		if (zmin > z) zmin = z;
		if (zmax < z) zmax = z;

		D.p.push_back(pair<Point, unsigned int>(Point(x, y, z), i));
		D.p[i].first.h = h;
		//ss >> D.p[i].first.x >> D.p[i].first.y >> D.p[i].first.z;
		//D.p[i].second = i;

	}

	getline(file, s);  // line "***** TRIANGLES *****"
	s.clear();

	D.t.clear();
	D.t.reserve(NUM_TRIANGLES);
	for (unsigned int i = 0; i < NUM_TRIANGLES; i++){

		Triangle tx;
		D.t.push_back(tx);

		double nx, ny, nz,h;
		unsigned int pid0, pid1, pid2;
		unsigned int neighbor0, neighbor1, neighbor2;


		ss.clear();
		getline(file, s);
		ss << s;
		ss >> nx >> ny >> nz;
		//ss >> D.t[i].nx >> D.t[i].ny >> D.t[i].nz;
		D.t[i].nx = nx;
		D.t[i].ny = ny;
		D.t[i].nz = nz;


		ss.clear();
		getline(file, s);
		ss << s;
		ss >> h;
		//ss >> D.t[i].h ;
		D.t[i].h = h;


		ss.clear();
		getline(file, s);
		ss << s;
		ss >> pid0>> pid1>> pid2;
		//ss >> D.t[i].pid[0] >> D.t[i].pid[1] >> D.t[i].pid[2];
		D.t[i].pid[0] = pid0;
		D.t[i].pid[1] = pid1;
		D.t[i].pid[2] = pid2;

		ss.clear();
		getline(file, s);
		ss << s;
		ss >> temp >> neighbor0 >> neighbor1 >> neighbor2;
		D.t[i].neighbor[0] = neighbor0;
		D.t[i].neighbor[1] = neighbor1;
		D.t[i].neighbor[2] = neighbor2;
		//ss >> D.t[i].neighbor[0] >> D.t[i].neighbor[1] >> D.t[i].neighbor[2];
					

		D.t[i].p0 = D.p[D.t[i].pid[0]].first;
		D.t[i].p1 = D.p[D.t[i].pid[1]].first;
		D.t[i].p2 = D.p[D.t[i].pid[2]].first;


	}


	getline(file, s);  // line "***** LIST OF TRINAGLES AT EACH POINT  *****"
	s.clear();

	D.v.clear();
	D.v.reserve(NUM_POINTS);
	for (unsigned int i = 0; i < NUM_POINTS; i++){

		vector <unsigned int> tv;
		

		ss.clear();
		getline(file, s);
		s = s + "-1";
		ss << s;
		int ii;
		while (ss){
			ss >> ii;
			if (ii!=-1) tv.push_back(ii);
						
		}
		D.v.push_back(tv);
		tv.clear();

		


	}


	double offsetx, offsety, offsetz;


	offsetx = abs(xmax - xmin)*0.05;
	offsety = abs(ymax - ymin)*0.05;
	offsetz = abs(zmax - zmin)*0.05;

	if (xmax == xmin) offsetx = (abs(ymax - ymin) < abs(zmax - zmin)) ? 0.05*abs(ymax - ymin) : 0.05*abs(zmax - zmin);
	if (ymax == ymin) offsety = (abs(xmax - xmin) < abs(zmax - zmin)) ? 0.05*abs(xmax - xmin) : 0.05*abs(zmax - zmin);
	if (zmax == zmin) offsetz = (abs(xmax - xmin) < abs(ymax - ymin)) ? 0.05*abs(xmax - xmin) : 0.05*abs(ymax - ymin);


	D.b = BoundingBox(xmin - offsetx, xmax + offsetx, ymin - offsety, ymax + offsety, zmin - offsetz, zmax + offsetz);

	file.close();

}

bool isTriWet(imp& IMP, unsigned int triId, double dimpact, pt_data &D){

	// store location of impact 
	Point p0 = D.t[triId].p0;
	Point p1 = D.t[triId].p1;
	Point p2 = D.t[triId].p2;

	Point pp = IMP.im; // impact point


	if (IMP.id == triId) return true; // if id of triangle where the particle impacts in the same as the triangle where checking then return true;
	if ((p0.x - pp.x)*(p0.x - pp.x) + (p0.y - pp.y)*(p0.y - pp.y) + (p0.z - pp.z)*(p0.z - pp.z) < dimpact*dimpact) return true;
	if ((p1.x - pp.x)*(p1.x - pp.x) + (p1.y - pp.y)*(p1.y - pp.y) + (p1.z - pp.z)*(p1.z - pp.z) < dimpact*dimpact) return true;
	if ((p2.x - pp.x)*(p2.x - pp.x) + (p2.y - pp.y)*(p2.y - pp.y) + (p2.z - pp.z)*(p2.z - pp.z) < dimpact*dimpact) return true;
	


	return false;
}



int main() {

	
	

	//char* fname = "Cyl.stl";

	//find_avglength(fname);
	//DATA = from_stl(fname);
	
	//Triangle* T = DATA.t;
	//Point* P = DATA.p;
	//vector<unsigned int>* V = DATA.v;
	//BoundingBox B = DATA.b;

	//find_neighbors(DATA.t, DATA.p, NUM_TRIANGLES);
	//triangle_report(DATA.t, DATA.p,DATA.v,  NUM_TRIANGLES);
	//write_file_triangles(DATA.t, DATA.p, 50000);
	//save("dump.dat",DATA);
	load("dump.dat", DATA);
	
	write_file_triangles(DATA, DATA.t, DATA.p, 50000);

	//kdtree KDTree(DATA.p, DATA.t, DATA.v, DATA.b);
	//for (unsigned int i=0; i<NUM_POINTS; i++)

	//octree(DATA.p, &DATA.v, DATA.b, 1);
	
	//octree O(DATA,1000);


	//cout << endl << "Press Enter ...";
	//getchar();


	// settings parameters
	coating::settings::Dmean = 58.0;
	coating::settings::Dstd = 18.0;

	coating::settings::Vmean = 60.0;
	coating::settings::Vstd = 5.0;

	coating::settings::Tmean = 1609.0;
	coating::settings::Tstd = 219.0;

	coating::settings::Omean = 0.0;
	coating::settings::Ostd = 50.0;
		

	// gun parameters

	coating::gun::x_g = 0.0;
	coating::gun::y_g = 0.0;
	coating::gun::z_g = 2.0;

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

	particle Particle;
	Point p;

	p.x = Particle.x_p=0.0;
	p.y = Particle.y_p=0.0;
	p.z = Particle.z_p=2.0;

	

	//tStart = cpuSecond();
	//cout << "*** Octree method ***" << endl;
	////r.find_CandidPoints(O);
	//r.find_CandidTriangles(O);
	////id= r.find_Intersection(O);
	////id = r.find_Intersection(O);
	//id = r.find_Tri_Intersection(O);
	//cout << "Ray hits triangle[" << id << "]" << endl ;
	//tStop = cpuSecond();
	//cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
	//cout << endl;

	cout << endl<<endl;

	//tStart = cpuSecond();
	//cout << "*** Brute Force ***" << endl;
	//id = brute_find(r, DATA);
	//cout << "Ray hits triangle[" << id << "]" << endl;
	//tStop = cpuSecond();
	//cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
	//cout << endl;

	tStart = cpuSecond();
	unsigned int numtime = 1000;
	double dt = coating::gun::injtime / numtime;
	double time = 0.0;
	double mass_total = coating::gun::injtime*coating::gun::mfrate;

	
	
	double tStart1 = cpuSecond();
	while (gi<numtime) {

		double mass_dt_init = coating::gun::mfrate*dt;
		double mass_dt = coating::gun::mfrate*dt;
		double mass_dt_old = mass_dt;

		gii = 0;
		

		while (mass_dt >(PI / 6.0)*(1E-6)*(1E-6)*(1.E-6)) { // while remaining mass is larger than mass of 1 micron particle

			Particle.reset();

			Vector v;
			v.vx = Particle.V_x;
			v.vy = Particle.V_y;
			v.vz = Particle.V_z;

			Vector vn;
			vn = v.normal();

			Ray r(p, vn);

			//int id = -1;
			imp im;


			im = brute_find(r, DATA);

			double eff_area ; // total effective wetted area
			vector<unsigned int> eff_tris; // this vector holds a list of all wet triangles
			vector<unsigned int> candid_tris; // this vector holds a list of all candidate wet triangles
			
			eff_area = 0.0;
			eff_tris.clear();
			candid_tris.clear();
						
			if ( im.id != UINT_MAX){

				eff_tris.push_back(im.id); // add impact triangle to the final list
				eff_area = eff_area + sqrt(DATA.t[im.id].area2());

				// add impact triangle neighbors to the candid list
				candid_tris.push_back(DATA.t[im.id].neighbor[0]); 
				candid_tris.push_back(DATA.t[im.id].neighbor[1]);
				candid_tris.push_back(DATA.t[im.id].neighbor[2]);
				
				Vector SurfaceNormal(DATA.t[im.id].nx, DATA.t[im.id].ny, DATA.t[im.id].nz);
				Particle.V_normal = v.dotProduct(SurfaceNormal);

				Particle.dmax_p = Particle.spreadfactor()*Particle.D_p;
				Particle.h_p = 2 * Particle.D_p / (3.0*Particle.spreadfactor()*Particle.spreadfactor());
				
				while (candid_tris.size() != 0){
					unsigned int id = candid_tris[0];

					if (!isTriWet(im, id, Particle.dmax_p, DATA)){
						// triangle is not wet therefore remove it from the candid list
						candid_tris.erase(remove(candid_tris.begin(), candid_tris.end(), id), candid_tris.end());

					}

					else {
						// triangle is wet, so first remove it from the candid list
						
						// 1 -removing from the candid list
						candid_tris.erase(remove(candid_tris.begin(), candid_tris.end(), id), candid_tris.end());
						
						// 2- if it's not already in the final list, then add it to the final list
						if (!isTri_inVector(id, eff_tris)) {
							
							eff_tris.push_back(id);
							eff_area = eff_area + sqrt(DATA.t[id].area2());

							// 3- if it's neghibors are not in the candid list or final list then add them to the candid list

							if (!isTri_inVector(DATA.t[im.id].neighbor[0], eff_tris) && !isTri_inVector(DATA.t[im.id].neighbor[0], candid_tris)) {
								candid_tris.push_back(DATA.t[im.id].neighbor[0]);
							}
							if (!isTri_inVector(DATA.t[im.id].neighbor[1], eff_tris) && !isTri_inVector(DATA.t[im.id].neighbor[1], candid_tris)) {
								candid_tris.push_back(DATA.t[im.id].neighbor[1]);
							}
							if (!isTri_inVector(DATA.t[im.id].neighbor[2], eff_tris) && !isTri_inVector(DATA.t[im.id].neighbor[2], candid_tris)) {
								candid_tris.push_back(DATA.t[im.id].neighbor[2]);
							}

						}

					}

				}

				

				for (unsigned int i = 0; i < eff_tris.size(); i++){
					DATA.t[eff_tris[i]].h = DATA.t[eff_tris[i]].h+ (Particle.h_p*(0.25*PI*Particle.D_p*Particle.D_p) / eff_area);
				}


				mass_dt -= (PI / 6.0)*Particle.D_p*Particle.D_p*Particle.D_p;
			}


			gii++;

			//if (!(gii % 10000)) {
			if (!(gii % 100)) {
				//cout << "*-----------------------------------*" << endl;
				cout << "\rparticle#: " << gii << "   mass/mass0%= " << (mass_dt_old / mass_dt_init)*100.0 << "   Dp= " << Particle.D_p*1E6 << " [micron]   Dmax= " << Particle.dmax_p*1E6 << " [micron]";
				//cout << "*-----------------------------------*" << endl;
			}

			if (mass_dt> (PI / 6.0)*(1E-6)*(1E-6)*(1.E-6)) {
				mass_dt_old = mass_dt;
			}


		}


		cout << "\rparticle#: " << gii << "   mass/mass0%= " << (mass_dt_old / mass_dt_init)*100.0 << "   Dp= " << Particle.D_p*1E6 << " [micron]   Dmax= " << Particle.dmax_p*1E6 << " [micron]";
		gi++;
		cout << endl << "Ntime= " << gi << " time= " << time << "  / endtime= " << coating::gun::injtime << endl;
		tStop = cpuSecond();
		
		cout << "elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
		cout << "*-------------------------------------------*" << endl;
		tStart = cpuSecond();



		//if (!(gi % 100) || gi == 1) {
		if (!(gi % 10) || gi == 1) {

			//write_file_triangles(t, P, gi);
			write_file_triangles(DATA, DATA.t, DATA.p, gi);
		}


		time += dt;
	}

	


	//cout << endl<<"Press Enter to exit ..."<<endl<<endl;
	//getchar();

	double tStop1 = cpuSecond();
	cout << "----->   elapsed time: " << double(tStop1 - tStart1) << " [milisec]" << endl;
	cout << endl;
	return 0;
}






