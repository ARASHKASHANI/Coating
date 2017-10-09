#pragma once


#include <iostream>
#include <iomanip>



#include <vector>
#include<set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include <climits>


#include <string>
#include <sstream>
#include <fstream>



#include<time.h>
#include <sys/timeb.h>
#include <iterator>



using namespace std;

double AVG_L = 0.0; // average edge length of triangles in the STL file
size_t TBL_SIZE; // table size for unordered_multimap
size_t TBL_SIZE_POINTS; // table size for unordered_multimap
size_t NUM_TRIANGLES; // total number of triangles
size_t NUM_POINTS_RAW; // total number of raw points
size_t NUM_POINTS; // total number of points --> duplicate points removed <--

double tStart, tStop;



double cpuSecond() {

	// returns time for CPU in miliseconds

	struct timeb myTime;

	ftime(&myTime);

	return double(1000 * myTime.time + myTime.millitm);
}


size_t find_nextprime(size_t num) {


	size_t v0, v1, v2, v3, v4, v5, v6, v7;

	ifstream f("primes.txt");

	string temp;
	stringstream stemp;





	while (!f.eof()) {

		getline(f, temp);
		stemp << temp;
		stemp >> v0 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7;
		stemp.clear();

		if (num <= v0 && v0<UINT_MAX) { f.close(); return v0; }
		if (num <= v1 && v1<UINT_MAX) { f.close(); return v1; }
		if (num <= v2 && v2<UINT_MAX) { f.close(); return v2; }
		if (num <= v3 && v3<UINT_MAX) { f.close(); return v3; }
		if (num <= v4 && v4<UINT_MAX) { f.close(); return v4; }
		if (num <= v5 && v5<UINT_MAX) { f.close(); return v5; }
		if (num <= v6 && v6<UINT_MAX) { f.close(); return v6; }
		if (num <= v7 && v7<UINT_MAX) { f.close(); return v7; }

	}


	f.close();
	return 49999; //default value



}



struct BoundingBox {
	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;

	BoundingBox() {
		xmin, xmax = 0.0;
		ymin, ymax = 0.0;
		zmin, zmax = 0.0;
	}

	BoundingBox(double xminp, double xmaxp, double yminp, double ymaxp, double zminp, double zmaxp) :
		xmin(xminp), xmax(xmaxp), ymin(yminp), ymax(ymaxp), zmin(zminp), zmax(zmaxp) {}
};



///////////////////////////////////////////////////////////////////////////////////////////////////////
class Point;
namespace std {
	template<>
	class hash<Point> {
	public:
		 size_t operator()(const Point &) const;  // declaration of operator() to allow befriending by Point
	};
}




class Point {
public:
	double x, y, z;
	double h;

	Point() { h = 0.0; }
	Point(const double& x0, const double& y0, const double & z0) : x(x0), y(y0), z(z0) { h = 0; }


	 friend bool operator== (const Point&, const Point&);
	 friend bool operator!= (const Point&, const Point&);
	 friend std::size_t std::hash<Point>::operator()(const Point &) const;


};


unordered_multimap < Point, int> points_map; // int here refers to  triangle id
unordered_multimap < Point, int> rawpoints_map; // int here refers to raw point id




class Vector {
public:
	double vx, vy, vz;

	 Vector() {}

	  Vector(const double& x0, const double&y0, const double& z0, const double& x1, const double&  y1, const double& z1)
	{
		vx = x1 - x0;
		vy = y1 - y0;
		vz = z1 - z0;


	}

	  Vector(const Point& p) {
		
		  vx = p.x ;
		  vy = p.y ;
		  vz = p.z ;


	  }

	  Vector(const Point& p0, const Point& p1) {

		vx = p1.x - p0.x;
		vy = p1.y - p0.y;
		vz = p1.z - p0.z;

	}

	  Vector(const double& vxx, const double& vyy, const double& vzz) : vx(vxx), vy(vyy), vz(vzz) {}

	  double magn() { return sqrt(vx*vx + vy*vy + vz*vz); }

	  Vector normal() {
		return Vector(vx / magn(), vy / magn(), vz / magn());
	}

	  double dotProduct(Vector &second) {
		  return vx*second.vx + vy*second.vy + vz*second.vz;
	}

	  Vector crossProduct(Vector &second) {

		  double nvx = 0, nvy = 0, nvz = 0;
		  nvx = vy*second.vz - vz*second.vy;
		  nvy = -vx*second.vz + vz*second.vx;
		  nvz = vx*second.vy - vy*second.vx;

		  return Vector(nvx, nvy, nvz);
	}

	  Vector operator+ (const Vector& second) {
		  Vector ret;
		  ret.vx = vx + second.vx;
		  ret.vy = vy + second.vy;
		  ret.vz = vz + second.vz;

		  return ret;
	  }

	  Vector operator- (const Vector& second) {
		  Vector ret;
		  ret.vx = vx - second.vx;
		  ret.vy = vy - second.vy;
		  ret.vz = vz - second.vz;

		  return ret;
	  }


};

class Line {
public:

	double x, y, z;
	double nx, ny, nz;

	  Line() {}

	  Line(const Point& p0, const Point& p1) {

		Vector v(p0, p1);

		x = p0.x;
		y = p0.y;
		z = p0.z;

		nx = v.normal().vx;
		ny = v.normal().vy;
		nz = v.normal().vz;


	}

	  Line(const Point&  p0, Vector& v) {


		x = p0.x;
		y = p0.y;
		z = p0.z;

		nx = v.normal().vx;
		ny = v.normal().vy;
		nz = v.normal().vz;


	}

	  Line(Vector& v, const Point& p0) {


		x = p0.x;
		y = p0.y;
		z = p0.z;

		nx = v.normal().vx;
		ny = v.normal().vy;
		nz = v.normal().vz;


	}

	  Line(const double& x0, const double& y0, const double& z0, const double& x1, const double& y1, const double& z1)
	{

		x = x0;
		y = y0;
		z = z0;

		Vector v(x0, y0, z0, x1, y1, z1);
		nx = v.normal().vx;
		ny = v.normal().vy;
		nz = v.normal().vz;


	}



};


class Triangle {
public:
	Point p0, p1, p2;
	double nx, ny, nz; // these are loaded from the stl file, since it has +/- values
	unsigned int neighbor[3]; // holds id of neighboring traingles
	unsigned int pid[3]; // holds id of points in the point array
	BoundingBox b; // bounding box of the triangle

	double h; // holds thickness of coating in each triangle
	
	//Triangle* parent;
	//Triangle* next;

	Triangle() { neighbor[0] = UINT_MAX; neighbor[1] = UINT_MAX; neighbor[2] = UINT_MAX; h = 0.0; /*parent = NULL; next = NULL;*/ }
	  Triangle(const Point& P0, const Point& P1, const Point& P2) : p0(P0), p1(P1), p2(P2)
	{
		  neighbor[0] = UINT_MAX; neighbor[1] = UINT_MAX; neighbor[2] = UINT_MAX; h = 0.0; /*parent = NULL; next = NULL;*/
	}


	  double perimeter() {
		double size0 = sqrt((p1.x - p0.x)*(p1.x - p0.x) + (p1.y - p0.y)*(p1.y - p0.y) + (p1.z - p0.z)*(p1.z - p0.z));
		double size1 = sqrt((p2.x - p0.x)*(p2.x - p0.x) + (p2.y - p0.y)*(p2.y - p0.y) + (p2.z - p0.z)*(p2.z - p0.z));
		double size2 = sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y) + (p2.z - p1.z)*(p2.z - p1.z));

		return size0 + size1 + size2;
	}

	  double area() {

		double size0 = sqrt((p1.x - p0.x)*(p1.x - p0.x) + (p1.y - p0.y)*(p1.y - p0.y) + (p1.z - p0.z)*(p1.z - p0.z));
		double size1 = sqrt((p2.x - p0.x)*(p2.x - p0.x) + (p2.y - p0.y)*(p2.y - p0.y) + (p2.z - p0.z)*(p2.z - p0.z));
		double size2 = sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y) + (p2.z - p1.z)*(p2.z - p1.z));

		double p = (size0 + size1 + size2) / 2.0;
		return sqrt(p*(p - size0)*(p - size1)*(p - size2));
	}

	  double area2() {

		  double size0 = sqrt((p1.x - p0.x)*(p1.x - p0.x) + (p1.y - p0.y)*(p1.y - p0.y) + (p1.z - p0.z)*(p1.z - p0.z));
		  double size1 = sqrt((p2.x - p0.x)*(p2.x - p0.x) + (p2.y - p0.y)*(p2.y - p0.y) + (p2.z - p0.z)*(p2.z - p0.z));
		  double size2 = sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y) + (p2.z - p1.z)*(p2.z - p1.z));

		  double p = (size0 + size1 + size2) / 2.0;
		  return (p*(p - size0)*(p - size1)*(p - size2));
	  }

	  void setBoundingBox() {
		  double xmin=p0.x, xmax=p0.x, ymin=p0.y, ymax=p0.y, zmin=p0.z, zmax=p0.z;

		  if (xmin > p1.x) xmin = p1.x;
		  if (xmin > p2.x) xmin = p2.x;

		  if (xmax < p1.x) xmax = p1.x;
		  if (xmax < p2.x) xmax = p2.x;

		  if (ymin > p1.y) ymin = p1.y;
		  if (ymin > p2.y) ymin = p2.y;

		  if (ymax < p1.y) ymax = p1.y;
		  if (ymax < p2.y) ymax = p2.y;

		  if (zmin > p1.z) zmin = p1.z;
		  if (zmin > p2.z) zmin = p2.z;

		  if (zmax < p1.z) zmax = p1.z;
		  if (zmax < p2.z) zmax = p2.z;


	  }

};


  inline bool operator== (const Point &pOne, const Point &pTwo) {

	if ((pOne.x == pTwo.x) && (pOne.y == pTwo.y) && (pOne.z == pTwo.z)) return true;
	else return false;

}
  inline bool operator!= (const Point &pOne, const Point &pTwo) {

	return !operator== (pOne, pTwo);


}
size_t hash<Point>::operator()(const Point &p) const {

	size_t id = ((int)(abs(p.x / AVG_L) * 73856093) ^ (int)(abs(p.y / AVG_L) * 19349663) ^ (int)(abs(p.z / AVG_L) * 83492791)) % TBL_SIZE;
	//if (id<0) cout << "id =" << id << endl;
	return id;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////

struct pt_data {
	vector<pair<Point, unsigned int>> p; // vector of pair (point, and their original id)
	vector<Triangle> t; // vector of triangles
	vector<vector<unsigned int>> v; // vector of vectors that hold list of triangles that are connected to a point 
	BoundingBox b; // counding box
		
};


void find_avglength(const char* fname) {

	tStart = cpuSecond();

	vector<Point> pv;
	//unordered_set<Point> s;
	//set<Point> s;



	cout << "readig file << " << fname << " >>" << endl;
	ifstream file;
	file.open(fname);


	int i;

	double nx, ny, nz;
	double x, y, z;
	string temp;
	stringstream stemp;

	// find out number of lines in the ascii file and also the total number of triangles	
	int count = 0;
	int count_point = 0;


	double AREA = 0.0;
	double l = 0.0; // characteristic length
	Triangle* tri_temp = new Triangle;



	double l0 = 0.0;
	double l1 = 0.0;
	double l2 = 0.0;


	i = 0;
	getline(file, temp); // reads "solid"
	while (!file.eof()) {


		getline(file, temp); // reads "facet normal nx ny nz"

		if (temp.substr(0,8) == "endsolid") break;
		stemp << temp;
		stemp >> temp >> temp >> nx >> ny >> nz;

		tri_temp->nx = nx;
		tri_temp->ny = ny;
		tri_temp->nz = nz;



		getline(file, temp); // reads "outer loop"

		getline(file, temp); // reads "vertex x y z"
		stemp.clear();
		stemp << temp;
		stemp >> temp >> x >> y >> z;

		tri_temp->p0.x = x;
		tri_temp->p0.y = y;
		tri_temp->p0.z = z;

		//pv.push_back(Point(x, y, z));

		count_point++;

		getline(file, temp); // reads "vertex x y z"
		stemp.clear();
		stemp << temp;
		stemp >> temp >> x >> y >> z;

		tri_temp->p1.x = x;
		tri_temp->p1.y = y;
		tri_temp->p1.z = z;

		//pv.push_back(Point(x, y, z));

		count_point++;


		getline(file, temp); // reads "vertex x y z"
		stemp.clear();
		stemp << temp;
		stemp >> temp >> x >> y >> z;

		tri_temp->p2.x = x;
		tri_temp->p2.y = y;
		tri_temp->p2.z = z;

		//pv.push_back(Point(x, y, z));

		count_point++;

		getline(file, temp); // reads "endloop"

		getline(file, temp); // reads "endfacet"



		count++;


		AREA += tri_temp->area();


		l0 += sqrt((tri_temp->p1.x - tri_temp->p0.x)*(tri_temp->p1.x - tri_temp->p0.x) + (tri_temp->p1.y - tri_temp->p0.y)*(tri_temp->p1.y - tri_temp->p0.y) + (tri_temp->p1.z - tri_temp->p0.z)*(tri_temp->p1.z - tri_temp->p0.z));
		l1 += sqrt((tri_temp->p2.x - tri_temp->p0.x)*(tri_temp->p2.x - tri_temp->p0.x) + (tri_temp->p2.y - tri_temp->p0.y)*(tri_temp->p2.y - tri_temp->p0.y) + (tri_temp->p2.z - tri_temp->p0.z)*(tri_temp->p2.z - tri_temp->p0.z));
		l2 += sqrt((tri_temp->p2.x - tri_temp->p1.x)*(tri_temp->p2.x - tri_temp->p1.x) + (tri_temp->p2.y - tri_temp->p1.y)*(tri_temp->p2.y - tri_temp->p1.y) + (tri_temp->p2.z - tri_temp->p1.z)*(tri_temp->p2.z - tri_temp->p1.z));




		i++;
	}

	l = sqrt(AREA / count);
	l = (l0 + l1 + l2) / (3 * count);

	NUM_TRIANGLES = count;
	NUM_POINTS_RAW = count_point;
	AVG_L = l;


	delete tri_temp;



	tStop = cpuSecond();



	//cout << endl;	cout << "________________________________________" << endl;

	cout << endl;
	cout << "Mini report:" << endl << endl;
	cout << "Number of Triangles =  " << NUM_TRIANGLES << endl;
	cout << "Number of raw points =  " << NUM_POINTS_RAW << endl;
	cout << "Avg Edge Length     =  " << AVG_L << endl;
	cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;




	file.close();


}


bool PointComparison(Point p0, Point p1) {
	if ((p0.x == p1.x) && (p0.y == p1.y) && (p0.z == p1.z)) return true;
	else { return false; }
}

struct sort_by_x {
	bool operator() (Point& p0, Point& p1) {
		return (p0.x < p1.x);
	}
};

struct sort_by_y {
	bool operator() (Point& p0, Point& p1) {
		return (p0.y < p1.y);
	}
};

struct sort_by_z {
	bool operator() (Point& p0, Point& p1) {
		return (p0.z < p1.z);
	}
};

struct sorter {
	bool operator() (Point& p0, Point& p1) {
		return (p0.x<p1.x) || (p0.x == p1.x && p0.y<p1.y) || (p0.x == p1.x && p0.y == p1.y && p0.z<p1.z);
	}
};

pt_data from_stl(const char* fname) {

	// this function returns a pointer to an array of triangles. The dimension of the array is determined from the facet entries in the stl file.


	//--------------------------------------------------------------------

	ifstream file;
	file.open(fname);

	double nx, ny, nz;
	double x, y, z;

	double xmin = 0.0;
	double xmax = 0.0;
	double ymin = 0.0;
	double ymax = 0.0;
	double zmin = 0.0;
	double zmax = 0.0;

	string temp;
	stringstream stemp;

	int count = 0;
	int count_point = 0;

	unsigned int i, j;
	vector<Point> pv;
	pv.reserve(NUM_POINTS_RAW);

	TBL_SIZE = find_nextprime(2 * NUM_TRIANGLES); // looking for the next prime number after 2*number of triangles
	points_map.reserve(TBL_SIZE); // set the hash table size



	//Triangle* tri = new Triangle[NUM_TRIANGLES];
	vector<Triangle> tri(NUM_TRIANGLES);
	
	//cout << "RawPoint Hash Table Size     =  " << TBL_SIZE_POINTS << endl;


	

	//cout << "________________________________________" << endl;
	cout << endl;
	//--------------------------------------------------------------------


	tStart = cpuSecond();

	cout << endl << "Step 1:";
	cout << endl << "Hashing triangle points:" << endl;
	cout << "Triangle Hash Table Size     =  " << TBL_SIZE << endl;
	// going back to the beginning of the file
	/*file.clear();
	file.seekg(0, ios::beg);*/


	int id = 0;
	getline(file, temp); // reads "solid"
	Point p;

	i = 0;
	j = 0;

	while (!file.eof()) {

		stemp.clear();

		getline(file, temp); // reads "facet normal nx ny nz"

		if (temp.substr(0, 8) == "endsolid") break;
		stemp << temp;
		//cout << "i= "<<i<<"   "<< temp<<endl;

		stemp >> temp  >> temp >> nx >> ny >> nz;


		tri[i].nx = nx;
		tri[i].ny = ny;
		tri[i].nz = nz;



		getline(file, temp); // reads "outer loop"

		getline(file, temp); // reads "vertex x y z"
		stemp.clear();
		stemp << temp;
		stemp >> temp >> x >> y >> z;

		tri[i].p0.x = x;
		tri[i].p0.y = y;
		tri[i].p0.z = z;

		if (xmin > tri[i].p0.x) xmin = tri[i].p0.x;
		if (xmax < tri[i].p0.x) xmax = tri[i].p0.x;

		if (ymin > tri[i].p0.y) ymin = tri[i].p0.y;
		if (ymax < tri[i].p0.y) ymax = tri[i].p0.y;

		if (zmin > tri[i].p0.z) zmin = tri[i].p0.z;
		if (zmax < tri[i].p0.z) zmax = tri[i].p0.z;
		

		p.x = x; p.y = y; p.z = z;
		pv.push_back(p);
		j++;


		getline(file, temp); // reads "vertex x y z"
		stemp.clear();
		stemp << temp;
		stemp >> temp >> x >> y >> z;

		tri[i].p1.x = x;
		tri[i].p1.y = y;
		tri[i].p1.z = z;

		if (xmin > tri[i].p1.x) xmin = tri[i].p1.x;
		if (xmax < tri[i].p1.x) xmax = tri[i].p1.x;

		if (ymin > tri[i].p1.y) ymin = tri[i].p1.y;
		if (ymax < tri[i].p1.y) ymax = tri[i].p1.y;

		if (zmin > tri[i].p1.z) zmin = tri[i].p1.z;
		if (zmax < tri[i].p1.z) zmax = tri[i].p1.z;


		p.x = x; p.y = y; p.z = z;
		pv.push_back(p);
		j++;

		getline(file, temp); // reads "vertex x y z"
		stemp.clear();
		stemp << temp;
		stemp >> temp >> x >> y >> z;

		tri[i].p2.x = x;
		tri[i].p2.y = y;
		tri[i].p2.z = z;

		tri[i].setBoundingBox();

		if (xmin > tri[i].p2.x) xmin = tri[i].p2.x;
		if (xmax < tri[i].p2.x) xmax = tri[i].p2.x;

		if (ymin > tri[i].p2.y) ymin = tri[i].p2.y;
		if (ymax < tri[i].p2.y) ymax = tri[i].p2.y;

		if (zmin > tri[i].p2.z) zmin = tri[i].p2.z;
		if (zmax < tri[i].p2.z) zmax = tri[i].p2.z;

		p.x = x; p.y = y; p.z = z;
		pv.push_back(p);
		j++;

		getline(file, temp); // reads "endloop"

		getline(file, temp); // reads "endfacet"


							 //cout << "\r   " << i;

		if (int(i*100.0 / NUM_TRIANGLES) % 10 == 0) {
			cout << "\r " << int(i*100.0/NUM_TRIANGLES) << "% complete ...     ";
		}

		i++;

	}

	//getline(file, temp); // reads "end solid"

	file.close();

	cout << "\r 100% complete ...     ";
	cout << endl;


	tStop = cpuSecond();

	cout << "original size = " << pv.size() << endl;

	//ofstream fp("pv-org.txt");
	//for (int i = 0; i < pv.size(); i++) {
	//	fp << i << "   " << pv[i].x << "   " << pv[i].y << "   " << pv[i].z << endl;
	//}
	//fp.close();

	sort(pv.begin(), pv.end(), sorter());


	vector<Point>::iterator it = unique(pv.begin(), pv.end(), PointComparison);
	pv.resize(distance(pv.begin(), it));

	//fp.open("pv-mod.txt");
	//for (int i = 0; i < pv.size(); i++) {
	//	fp << i << "   " << pv[i].x << "   " << pv[i].y << "   " << pv[i].z << endl;
	//}
	//fp.close();



	cout << "modified size = " << pv.size() << endl;



	cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;

	tStart = cpuSecond();

	cout << endl << "Step 2:";
	cout << endl << "Making map of Points and Id of Triangles containing them ..." ;


	NUM_POINTS = pv.size();
	TBL_SIZE_POINTS = find_nextprime(2 * NUM_POINTS); // looking for the next prime number after 2*number of raw points
	rawpoints_map.reserve(TBL_SIZE_POINTS);

	i = 0;
	for (Point pp : pv) {
		rawpoints_map.insert(pair<Point, int>(pp, i));
		//cout << "\r" << i;
		i++;
	}
	cout << endl;

	
	for (i = 0; i < NUM_TRIANGLES; i++) {
		auto range = rawpoints_map.equal_range(tri[i].p0);
		for (auto it = range.first; it != range.second; ++it) {

			Point c = it->first;
			int val = it->second;

			//cout << "tri[i].p0.x= " << tri[i].p0.x << "  tri[i].p0.y= " << tri[i].p0.y<<"   tri[i].p0.z= "<< tri[i].p0.z<<endl;
			//cout << "        c.x= " <<         c.x << "          c.y= " <<         c.y<<"           c.z= " << c.z << endl;

			if (tri[i].p0 == c) {
				tri[i].pid[0] = val;
				points_map.insert(pair<Point, int>(tri[i].p0, i));
				break;
			}
		}

		range = rawpoints_map.equal_range(tri[i].p1);
		for (auto it = range.first; it != range.second; ++it) {

			Point c = it->first;
			int val = it->second;

			//cout << "tri[i].p1.x= " << tri[i].p1.x << "  tri[i].p1.y= " << tri[i].p1.y << "   tri[i].p1.z= " << tri[i].p1.z << endl;
			//cout << "        c.x= " << c.x << "          c.y= " << c.y << "           c.z= " << c.z << endl;

			if (tri[i].p1 == c) {
				tri[i].pid[1] = val;
				points_map.insert(pair<Point, int>(tri[i].p1, i));
				break;
			}

		}

		range = rawpoints_map.equal_range(tri[i].p2);
		for (auto it = range.first; it != range.second; ++it) {

			Point c = it->first;
			int val = it->second;

			//cout << "tri[i].p2.x= " << tri[i].p2.x << "  tri[i].p2.y= " << tri[i].p2.y << "   tri[i].p2.z= " << tri[i].p2.z << endl;
			//cout << "        c.x= " << c.x << "          c.y= " << c.y << "           c.z= " << c.z << endl;


			if (tri[i].p2 == c) {
				tri[i].pid[2] = val;
				points_map.insert(pair<Point, int>(tri[i].p2, i));
				break;
			}

		}


		if (int(i*100.0 / NUM_TRIANGLES) % 10 == 0) {
			cout << "\r " << int(i*100.0 / NUM_TRIANGLES) << "% complete ...     ";
		}

		//cout << "\r" << i;
	}

	cout << "\r 100% complete ...     ";
	cout << endl;


	tStop = cpuSecond();
	cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl << endl;
	
	tStart = cpuSecond();

	cout << endl << "Step 3:";
	cout << endl  << "Making a list of Traingles that share a each single Point ..." << endl;
	//for (i = 0; i < NUM_TRIANGLES; i++) {
	//	vector<Point>::iterator it= find(pv.begin(), pv.end(),pv[0]);
	//	auto pos = distance(pv.begin(), it);
	//	tri[i].pid[0] = pos;

	//	it = find(pv.begin(), pv.end(), pv[1]);
	//	pos = distance(pv.begin(), it);
	//	tri[i].pid[1] = pos;

	//	it = find(pv.begin(), pv.end(), pv[2]);
	//	pos = distance(pv.begin(), it);
	//	tri[i].pid[2] = pos;


	//}
	//cout << endl << "here!";

	vector<vector<unsigned int>> points_tris (NUM_POINTS); // an array of vectors that holds list of triangles that shares one point pv[i]
	//points_tris = new vector<unsigned int>[NUM_POINTS]; // the size of this array is NUM_POINTS

	for (size_t i = 0; i < NUM_POINTS; i++) {
		auto range = points_map.equal_range(pv[i]);
		for (auto it = range.first; it != range.second; ++it) {

			Point c = it->first;
			int val = it->second;

			if (c == pv[i]) {
				points_tris[i].push_back(val);
			}


		}


		if (int(i*100.0 / NUM_POINTS) % 10 == 0) {
			cout << "\r " << int(i*100.0 / NUM_POINTS) << "% complete ...     ";
		}

	}


	cout << "\r 100% complete ...     ";
	cout << endl;

	tStop = cpuSecond();
	cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl<<endl;


	pt_data mydata;
	//mydata.p = pv.data();

	//Point* outp = new Point[NUM_POINTS];
	//for (unsigned int i = 0; i < pv.size(); i++) {
	//	outp[i] = pv[i];
	//}

	mydata.p.reserve(NUM_POINTS);
	for (unsigned int i = 0; i < NUM_POINTS; i++) {
		
		mydata.p.push_back(pair<Point,unsigned int>(pv[i],i));
					
	}

	//mydata.p = pv;
	mydata.t = tri;
	mydata.v = points_tris;

	double offsetx, offsety, offsetz;
	
	
	offsetx = abs(xmax - xmin)*0.05;
	offsety = abs(ymax - ymin)*0.05;
	offsetz = abs(zmax - zmin)*0.05;
	
	if (xmax == xmin) offsetx = (abs(ymax - ymin) < abs(zmax - zmin)) ? 0.05*abs(ymax - ymin) : 0.05*abs(zmax - zmin);
	if (ymax == ymin) offsety = (abs(xmax - xmin) < abs(zmax - zmin)) ? 0.05*abs(xmax - xmin) : 0.05*abs(zmax - zmin);
	if (zmax == zmin) offsetz = (abs(xmax - xmin) < abs(ymax - ymin)) ? 0.05*abs(xmax - xmin) : 0.05*abs(ymax - ymin);


	mydata.b = BoundingBox(xmin-offsetx, xmax+offsetx, ymin-offsety, ymax+offsety, zmin-offsetz, zmax+offsetz);

	
	
	return mydata;

}




vector<size_t> instersection(vector<size_t> &v1, vector<size_t> &v2)
{

	vector<size_t> v3;


	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));

	if (v3.size() == 0) {
		v3.push_back(UINT_MAX);
	}



	return v3;
}


void find_neighbors(vector<Triangle> &t, vector<pair<Point,unsigned int> >& p, size_t n) {

	tStart = cpuSecond();

	cout << endl << "... now finding neighbors ..." << endl;



	// looping through all triangles and setting id of neighbors
	for (size_t i = 0; i < n; i++) {


		// try with point p0 and p1

		//Point c(0,0,0);
		//int val;

		vector<Point> vp0, vp1, vp2;
		vector<size_t> id0, id1, id2;



		auto range = points_map.equal_range(t[i].p0);
		//auto range = rawpoints_map.equal_range(p[t[i].pid[0]]);
		for (auto it = range.first; it != range.second; ++it) {

			Point c = it->first;
			int val = it->second;

			if (val != i) {
				vp0.push_back(c);
				id0.push_back(val);
			}


		}

		range = points_map.equal_range(t[i].p1);
		//range = rawpoints_map.equal_range(p[t[i].pid[1]]);
		for (auto it = range.first; it != range.second; ++it) {

			Point c = it->first;
			int val = it->second;

			if (val != i) {
				vp1.push_back(c);
				id1.push_back(val);
			}


		}

		range = points_map.equal_range(t[i].p2);
		//range = rawpoints_map.equal_range(p[t[i].pid[2]]);
		for (auto it = range.first; it != range.second; ++it) {

			Point c = it->first;
			int val = it->second;

			if (val != i) {
				vp2.push_back(c);
				id2.push_back(val);
			}


		}

		t[i].neighbor[0] = unsigned int(instersection(id0, id1)[0]);

		t[i].neighbor[1] = unsigned int(instersection(id0, id2)[0]);

		t[i].neighbor[2] = unsigned int(instersection(id1, id2)[0]);



		if (int(i*100.0 / NUM_TRIANGLES) % 10 == 0) {
			cout << "\r " << int(i*100.0 / NUM_TRIANGLES) << "% complete ...     ";
		}






	}

	cout << "\r 100% complete ...     ";
	cout << endl;

	tStop = cpuSecond();

	cout << "Done!" << endl;

	cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
	cout << endl;
}

void triangle_report(vector<Triangle>& t, vector<pair <Point, unsigned int >>& p, vector<vector<unsigned int>>& v, size_t n) {


	tStart = cpuSecond();

	cout << "<<< writing triangle report file >>>" << endl;

	ofstream f("tri-neighbors.txt");
	ofstream ff("tri-coordinates.txt");
	ofstream fff("triangles.txt");
	ofstream ffff("points-tri.txt");
	ofstream fp("points.txt");
	
	fff << NUM_TRIANGLES << endl;
	fff << setprecision(14);


	for (unsigned int i = 0; i < n; i++) {
		f << "Triangle[" << i << "] :  "; // << endl;
										  /*		f << "P0:  " << "x= " << t[i].p0.x << "  y= " << t[i].p0.y << "  z=" << t[i].p0.z << endl;
										  f << "P1:  " << "x= " << t[i].p1.x << "  y= " << t[i].p1.y << "  z=" << t[i].p1.z << endl;
										  f << "P2:  " << "x= " << t[i].p2.x << "  y= " << t[i].p2.y << "  z=" << t[i].p2.z << endl;
										  f << endl;
										  f << "neighbor[0]= " << t[i].neighbor[0] << "  neighbor[1]= " << t[i].neighbor[1] << "  neighbor[2]= " << t[i].neighbor[2] << endl;
										  f << endl;    */

		f << "{ " << t[i].neighbor[0] << ", " << t[i].neighbor[1] << ", " << t[i].neighbor[2] << " }" << endl;

		ff << "[" << i << "] ->  p0: " << "{ x=  " << t[i].p0.x << ", y=  " << t[i].p0.y << ", z=  " << t[i].p0.z << " }" << endl;
		ff << "[" << i << "] ->  p1: " << "{ x=  " << t[i].p1.x << ", y=  " << t[i].p1.y << ", z=  " << t[i].p1.z << " }" << endl;
		ff << "[" << i << "] ->  p2: " << "{ x=  " << t[i].p2.x << ", y=  " << t[i].p2.y << ", z=  " << t[i].p2.z << " }" << endl;

		fff << i << "   " << t[i].nx << "   " << t[i].ny << "   " << t[i].nz << "   " <<
			/*t[i].p0.x << "   " << t[i].p0.y << "   " << t[i].p0.z << "   " <<
			t[i].p1.x << "   " << t[i].p1.y << "   " << t[i].p1.z << "   " <<
			t[i].p2.x << "   " << t[i].p2.y << "   " << t[i].p2.z << "   " <<
			t[i].pid[0] << "   " << t[i].pid[1] << "   " << t[i].pid[2] << "   " <<*/
			p[t[i].pid[0]].first.x << "   " << p[t[i].pid[0]].first.y << "   " << p[t[i].pid[0]].first.z << "   " <<
			p[t[i].pid[1]].first.x << "   " << p[t[i].pid[1]].first.y << "   " << p[t[i].pid[1]].first.z << "   " <<
			p[t[i].pid[2]].first.x << "   " << p[t[i].pid[2]].first.y << "   " << p[t[i].pid[2]].first.z << "   " <<
			t[i].pid[0] << "   " << t[i].pid[1] << "   " << t[i].pid[2] << "   " <<
			t[i].neighbor[0] << "   " << t[i].neighbor[1] << "   " << t[i].neighbor[2] << "   " << endl;

	}

	fff.close();
	ff.close();
	f.close();

	for (unsigned int i = 0; i < NUM_POINTS; i++) {
		ffff << "Point[" << i << "] = {";
		fp << "Point[" << i << "] = { x= " << p[i].first.x << " y= " << p[i].first.y << " z= " << p[i].first.z << " }"<< endl;

		for (unsigned int ii : v[i]) {
			ffff << ii << ",   ";
		}
		ffff << " }" << endl;
	}
	
	ffff.close();
	fp.close();

	tStop = cpuSecond();


	cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
	cout << endl;

}





void write_file_triangles(pt_data& D, vector<Triangle>& t_arg, vector< pair<Point, unsigned int >>& p_arg, unsigned int gi_arg) {
	ofstream f;
	stringstream filename;
	filename << "";
	filename << "output--" << gi_arg << ".dat";


	f.open(filename.str());

	f << "TITLE = \"Thickness plot\"" << endl;
	f << "VARIABLES = \"X\", \"Y\", \"Z\", \"h\" \"Z'\" " << endl;
	f << "ZONE N = " << NUM_POINTS << ", E = " << NUM_TRIANGLES << ", F = FEPOINT,  ET=TRIANGLE " << endl;


	for (unsigned int i = 0; i < NUM_POINTS; i++) {

		//p_arg[i].first.h = 0.0;
		unsigned int trueSize = 0; // 
		for (unsigned int ii = 0; ii < D.v[i].size(); ii++){
			if (D.v[i][ii] != UINT_MAX){
				p_arg[i].first.h = p_arg[i].first.h + D.t[D.v[i][ii]].h;
				trueSize++;
			}
		}
		//p_arg[i].first.h = p_arg[i].first.h / D.v[i].size();
		p_arg[i].first.h = p_arg[i].first.h / trueSize;



		f << p_arg[i].first.x << "   " << p_arg[i].first.y << "   " << p_arg[i].first.z << "   " << p_arg[i].first.h << "   " << p_arg[i].first.z + p_arg[i].first.h << endl;

	}

	for (unsigned int i = 0; i < NUM_TRIANGLES; i++) {
		f << t_arg[i].pid[0] + 1 << "   " << t_arg[i].pid[1] + 1 << "   " << t_arg[i].pid[2] + 1 << endl;

	}

	f.flush();
	f.close();

}





///////////////////////////////////////////////////////////////////////////////////////////////////////



