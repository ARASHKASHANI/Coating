#pragma once
#include "geometry.h"

struct sorter_pair {
	bool operator() (pair<Point, unsigned int>& p0, pair<Point, unsigned int>& p1) {
		return (p0.first.x<p1.first.x) || (p0.first.x == p1.first.x && p0.first.y<p1.first.y) || (p0.first.x == p1.first.x && p0.first.y == p1.first.y && p0.first.z<p1.first.z);
	}
};

bool PointComparison_pair(pair<Point, unsigned int>& p0, pair<Point, unsigned int>& p1) {
	if ((p0.first.x == p1.first.x) && (p0.first.y == p1.first.y) && (p0.first.z == p1.first.z)) return true;
	else { return false; }
}


class Node {

public:

	Node* parent = NULL;
	Node *ijk = NULL, *ipjk = NULL, *ijpk = NULL, *ipjpk = NULL;
	Node *ijkp = NULL, *ipjkp = NULL, *ijpkp = NULL, *ipjpkp = NULL;
	Node* next = NULL;
	bool divided = false;


	size_t numPoints; // number of points enclosed in this region
	vector<pair<Point, unsigned int>> vp; // a vector that holds a vector of Points enclosed in this region
	vector<vector<unsigned int>>* triId; // a vector that holds id of all the Triangles read from the STL file

	vector<unsigned int> triList; // a vector that holds id of Triangle enclosed in the current node

	BoundingBox bbox; // bounding box of enclosing region

	//Node(vector<vector<unsigned int>> &t) : triId(t) {  }
	Node(vector<vector<unsigned int>> *t) : triId(t) {  }



};



class octree {

public:
	unsigned int numcell; // maximum number of allowed points per cell
	Node* Root;
	Node* current;
	vector<Triangle>* tri;
	vector<vector<unsigned int>>* pt;
	vector<pair<Point, unsigned int>>* pp;




	//octree(vector<pair<Point, unsigned int>>& P, vector<vector<unsigned int>>& V, BoundingBox B, unsigned int n) {
	octree(vector<pair<Point, unsigned int>>& P, vector<Triangle>* T, vector<vector<unsigned int>>* V, BoundingBox B, unsigned int n) {

		create(P, T, V, B, n);




	}

	octree(pt_data& Data, unsigned int n) {
		create(Data.p, &Data.t, &Data.v, Data.b, n);

	}

	~octree() {

		tStart = cpuSecond();

		cout << "Deleting leaves of the octree" << endl;

		Node* Ntemp;

		current = Root;

		unsigned int layer = 1;

		ofstream f("octree-delete.txt");

		while (current != NULL) {



			if (current->numPoints > numcell && current->divided == false) {
				//divide(); current->divided = true;
				//if (current != Root) current->vp.clear();
			}


			//if (current->numPoints <= numcell) {
			//write(f, current);
			//}

			if (current->numPoints > numcell) {

				//write(f, current);

				if (current->ijk->numPoints <= numcell) {
					current->ijk->~Node();
					//write(f, current->ijk);
				}
				if (current->ijk->numPoints > numcell) {
					current = current->ijk;
					continue;
				}
				if (current->ipjk->numPoints <= numcell) {
					current->ipjk->~Node();
					//write(f, current->ipjk);
				}
				if (current->ipjk->numPoints > numcell) {
					current = current->ipjk;
					continue;
				}
				if (current->ijpk->numPoints <= numcell) {
					current->ijpk->~Node();
					//write(f, current->ijpk);
				}
				if (current->ijpk->numPoints > numcell) {
					current = current->ijpk;
					continue;
				}
				if (current->ipjpk->numPoints <= numcell) {
					current->ipjpk->~Node();
					//write(f, current->ipjpk);
				}
				if (current->ipjpk->numPoints > numcell) {
					current = current->ipjpk;
					continue;
				}
				if (current->ijkp->numPoints <= numcell) {
					current->ijkp->~Node();
					//write(f, current->ijkp);
				}
				if (current->ijkp->numPoints > numcell) {
					current = current->ijkp;
					continue;
				}
				if (current->ipjkp->numPoints <= numcell) {
					current->ipjkp->~Node();
					//write(f, current->ipjkp);
				}
				if (current->ipjkp->numPoints > numcell) {
					current = current->ipjkp;
					continue;
				}
				if (current->ijpkp->numPoints <= numcell) {
					current->ijpkp->~Node();
					//write(f, current->ijpkp);
				}
				if (current->ijpkp->numPoints > numcell) {
					current = current->ijpkp;
					continue;
				}
				if (current->ipjpkp->numPoints <= numcell) {
					current->ipjpkp->~Node();
					//write(f, current->ipjpkp);
				}
				if (current->ipjpkp->numPoints > numcell) {
					current = current->ipjpkp;
					continue;
				}
			}


			Ntemp = current;
			current = current->next;
			Ntemp->~Node();

			layer++;
			cout << "\r nuber of layer(s)= " << layer;


		}

		f.close();

		tStop = cpuSecond();

		cout << endl << "DONE!" << endl;
		cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
		cout << endl;

		Root = NULL;

	}

	void create(vector<pair<Point, unsigned int>>& P, vector<Triangle>* T, vector<vector<unsigned int>>* V, BoundingBox B, unsigned int n) {

		numcell = n;

		cout << "creating and populating octree" << endl;

		tStart = cpuSecond();

		Root = new Node(V);
		Root->parent = NULL;
		Root->numPoints = NUM_POINTS;
		Root->bbox = B;
		Root->triId = V;
		Root->vp = P;
		Root->next = NULL;
		Root->divided = false;

		Root->triList .reserve((*T).size()); // the Root triList contains all the triangles in the field
		for (unsigned int i = 0; i < (*T).size(); i++) Root->triList.push_back(i);

		tri = T;
		pt = V;
		pp = &P;

		current = Root;

		unsigned int layer = 1;

		ofstream f("octree-points.txt");

		while (current != NULL) {



			if (current->numPoints > numcell && current->divided == false) {
				divide(); current->divided = true;
				if (current != Root) current->vp.clear(); current->triList.clear();
			}


			/*if (current->numPoints <= numcell) {
			write(f, current);
			}*/

			if (current->numPoints > numcell) {

				//write(f, current);

				/*if (current->ijk->numPoints <= numcell) {
				write(f, current->ijk);
				}*/
				if (current->ijk->numPoints > numcell) {
					current = current->ijk;
					continue;
				}
				/*if (current->ipjk->numPoints <= numcell) {
				write(f, current->ipjk);
				}*/
				if (current->ipjk->numPoints > numcell) {
					current = current->ipjk;
					continue;
				}
				/*if (current->ijpk->numPoints <= numcell) {
				write(f, current->ijpk);
				}*/
				if (current->ijpk->numPoints > numcell) {
					current = current->ijpk;
					continue;
				}
				/*if (current->ipjpk->numPoints <= numcell) {
				write(f, current->ipjpk);
				}*/
				if (current->ipjpk->numPoints > numcell) {
					current = current->ipjpk;
					continue;
				}
				/*if (current->ijkp->numPoints <= numcell) {
				write(f, current->ijkp);
				}*/
				if (current->ijkp->numPoints > numcell) {
					current = current->ijkp;
					continue;
				}
				/*if (current->ipjkp->numPoints <= numcell) {
				write(f, current->ipjkp);
				}*/
				if (current->ipjkp->numPoints > numcell) {
					current = current->ipjkp;
					continue;
				}
				/*if (current->ijpkp->numPoints <= numcell) {
				write(f, current->ijpkp);
				}*/
				if (current->ijpkp->numPoints > numcell) {
					current = current->ijpkp;
					continue;
				}
				/*if (current->ipjpkp->numPoints <= numcell) {
				write(f, current->ipjpkp);
				}*/
				if (current->ipjpkp->numPoints > numcell) {
					current = current->ipjpkp;
					continue;
				}
			}



			current = current->next;

			layer++;
			cout << "\r nuber of layer(s)= " << layer;


		}

		f.close();

		tStop = cpuSecond();

		cout << endl << "DONE!" << endl;
		cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
		cout << endl;

	}

	void write(ofstream& f, Node* c) {

		f << c->numPoints << "  : { ";
		if (c->numPoints == 0) f << " } " << endl;
		else {
			for (unsigned int i = 0; i < c->numPoints; i++) {
				f << c->vp[i].second << " ,";
			}
			f << " } " << endl;
		}

	}

	void divide() {

		double xhalf = (current->bbox.xmax + current->bbox.xmin)*0.5;
		double yhalf = (current->bbox.ymax + current->bbox.ymin)*0.5;
		double zhalf = (current->bbox.zmax + current->bbox.zmin)*0.5;


		vector<pair<Point, unsigned int>>::iterator it;


		// ********************** cube ijk ********************** //
		current->ijk = new Node(current->triId);

		current->ijk->parent = current;



		current->ijk->bbox.xmin = current->bbox.xmin;
		current->ijk->bbox.xmax = xhalf;
		current->ijk->bbox.ymin = current->bbox.ymin;
		current->ijk->bbox.ymax = yhalf;
		current->ijk->bbox.zmin = current->bbox.zmin;
		current->ijk->bbox.zmax = zhalf;

		

		// now find the points in cube ijk

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x <= xhalf) && (p.first.y <= yhalf) && (p.first.z <= zhalf)) {
				current->ijk->vp.push_back(p); // update vector of points for cube ijk

			}


		}

		// keep particles sorted ?
		sort(current->ijk->vp.begin(), current->ijk->vp.end(), sorter_pair());
		it = unique(current->ijk->vp.begin(), current->ijk->vp.end(), PointComparison_pair);
		current->ijk->vp.resize(distance(current->ijk->vp.begin(), it));

		current->ijk->numPoints = current->ijk->vp.size();

		current->ijk->divided = false;

		// ********************** cube ijk ********************** //



		// **********************************************************************************************************************************


		// ********************** cube ipjk ********************** //
		current->ipjk = new Node(current->triId);

		current->ipjk->parent = current;



		current->ipjk->bbox.xmin = xhalf;
		current->ipjk->bbox.xmax = current->bbox.xmax;
		current->ipjk->bbox.ymin = current->bbox.ymin;
		current->ipjk->bbox.ymax = yhalf;
		current->ipjk->bbox.zmin = current->bbox.zmin;
		current->ipjk->bbox.zmax = zhalf;

		// now find the points in cube ipjk

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x > xhalf) && (p.first.y <= yhalf) && (p.first.z <= zhalf)) {
				current->ipjk->vp.push_back(p); // update vector of points for cube ipjk

			}


		}

		// keep particles sorted ?
		sort(current->ipjk->vp.begin(), current->ipjk->vp.end(), sorter_pair());
		it = unique(current->ipjk->vp.begin(), current->ipjk->vp.end(), PointComparison_pair);
		current->ipjk->vp.resize(distance(current->ipjk->vp.begin(), it));

		current->ipjk->numPoints = current->ipjk->vp.size();

		current->ipjk->divided = false;

		// ********************** cube ipjk ********************** //


		// **********************************************************************************************************************************


		// ********************** cube ijpk ********************** //
		current->ijpk = new Node(current->triId);

		current->ijpk->parent = current;



		current->ijpk->bbox.xmin = current->bbox.xmin;
		current->ijpk->bbox.xmax = xhalf;
		current->ijpk->bbox.ymin = yhalf;
		current->ijpk->bbox.ymax = current->bbox.ymax;
		current->ijpk->bbox.zmin = current->bbox.zmin;
		current->ijpk->bbox.zmax = zhalf;



		// now find the points in cube ijpk

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x <= xhalf) && (p.first.y > yhalf) && (p.first.z <= zhalf)) {
				current->ijpk->vp.push_back(p); // update vector of points for cube ijpk

			}


		}

		// keep particles sorted ?
		sort(current->ijpk->vp.begin(), current->ijpk->vp.end(), sorter_pair());
		it = unique(current->ijpk->vp.begin(), current->ijpk->vp.end(), PointComparison_pair);
		current->ijpk->vp.resize(distance(current->ijpk->vp.begin(), it));

		current->ijpk->numPoints = current->ijpk->vp.size();

		current->ijpk->divided = false;

		// ********************** cube ijpk ********************** //


		// **********************************************************************************************************************************


		// ********************** cube ipjpk ********************** //
		current->ipjpk = new Node(current->triId);

		current->ipjpk->parent = current;



		current->ipjpk->bbox.xmin = xhalf;
		current->ipjpk->bbox.xmax = current->bbox.xmax;
		current->ipjpk->bbox.ymin = yhalf;
		current->ipjpk->bbox.ymax = current->bbox.ymax;
		current->ipjpk->bbox.zmin = current->bbox.zmin;
		current->ipjpk->bbox.zmax = zhalf;



		// now find the points in cube ipjpk

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x > xhalf) && (p.first.y > yhalf) && (p.first.z <= zhalf)) {
				current->ipjpk->vp.push_back(p); // update vector of points for cube ipjpk

			}


		}

		// keep particles sorted ?
		sort(current->ipjpk->vp.begin(), current->ipjpk->vp.end(), sorter_pair());
		it = unique(current->ipjpk->vp.begin(), current->ipjpk->vp.end(), PointComparison_pair);
		current->ipjpk->vp.resize(distance(current->ipjpk->vp.begin(), it));

		current->ipjpk->numPoints = current->ipjpk->vp.size();

		current->ijpk->divided = false;

		// ********************** cube ipjpk ********************** //








		// ********************** cube ijkp ********************** //
		current->ijkp = new Node(current->triId);

		current->ijkp->parent = current;



		current->ijkp->bbox.xmin = current->bbox.xmin;
		current->ijkp->bbox.xmax = xhalf;
		current->ijkp->bbox.ymin = current->bbox.ymin;
		current->ijkp->bbox.ymax = yhalf;
		current->ijkp->bbox.zmin = zhalf;
		current->ijkp->bbox.zmax = current->bbox.zmax;

		// now find the points in cube ijkp

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x <= xhalf) && (p.first.y <= yhalf) && (p.first.z > zhalf)) {
				current->ijkp->vp.push_back(p); // update vector of points for cube ijkp

			}


		}

		// keep particles sorted ?
		sort(current->ijkp->vp.begin(), current->ijkp->vp.end(), sorter_pair());
		it = unique(current->ijkp->vp.begin(), current->ijkp->vp.end(), PointComparison_pair);
		current->ijkp->vp.resize(distance(current->ijkp->vp.begin(), it));

		current->ijkp->numPoints = current->ijkp->vp.size();

		current->ijkp->divided = false;

		// ********************** cube ijkp ********************** //



		// **********************************************************************************************************************************


		// ********************** cube ipjkp ********************** //
		current->ipjkp = new Node(current->triId);

		current->ipjkp->parent = current;



		current->ipjkp->bbox.xmin = xhalf;
		current->ipjkp->bbox.xmax = current->bbox.xmax;
		current->ipjkp->bbox.ymin = current->bbox.ymin;
		current->ipjkp->bbox.ymax = yhalf;
		current->ipjkp->bbox.zmin = zhalf;
		current->ipjkp->bbox.zmax = current->bbox.zmax;

		// now find the points in cube ipjkp

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x > xhalf) && (p.first.y <= yhalf) && (p.first.z > zhalf)) {
				current->ipjkp->vp.push_back(p); // update vector of points for cube ipjkp

			}


		}

		// keep particles sorted ?
		sort(current->ipjkp->vp.begin(), current->ipjkp->vp.end(), sorter_pair());
		it = unique(current->ipjkp->vp.begin(), current->ipjkp->vp.end(), PointComparison_pair);
		current->ipjkp->vp.resize(distance(current->ipjkp->vp.begin(), it));

		current->ipjkp->numPoints = current->ipjkp->vp.size();

		current->ipjkp->divided = false;

		// ********************** cube ipjkp ********************** //


		// **********************************************************************************************************************************


		// ********************** cube ijpkp ********************** //
		current->ijpkp = new Node(current->triId);

		current->ijpkp->parent = current;



		current->ijpkp->bbox.xmin = current->bbox.xmin;
		current->ijpkp->bbox.xmax = xhalf;
		current->ijpkp->bbox.ymin = yhalf;
		current->ijpkp->bbox.ymax = current->bbox.ymax;
		current->ijpkp->bbox.zmin = zhalf;
		current->ijpkp->bbox.zmax = current->bbox.zmax;



		// now find the points in cube ijpkp

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x <= xhalf) && (p.first.y > yhalf) && (p.first.z > zhalf)) {
				current->ijpkp->vp.push_back(p); // update vector of points for cube ijpkp

			}


		}

		// keep particles sorted ?
		sort(current->ijpkp->vp.begin(), current->ijpkp->vp.end(), sorter_pair());
		it = unique(current->ijpkp->vp.begin(), current->ijpkp->vp.end(), PointComparison_pair);
		current->ijpkp->vp.resize(distance(current->ijpkp->vp.begin(), it));

		current->ijpkp->numPoints = current->ijpkp->vp.size();

		current->ijpkp->divided = false;


		// ********************** cube ijpkp ********************** //


		// **********************************************************************************************************************************


		// ********************** cube ipjpkp ********************** //
		current->ipjpkp = new Node(current->triId);

		current->ipjpkp->parent = current;



		current->ipjpkp->bbox.xmin = xhalf;
		current->ipjpkp->bbox.xmax = current->bbox.xmax;
		current->ipjpkp->bbox.ymin = yhalf;
		current->ipjpkp->bbox.ymax = current->bbox.ymax;
		current->ipjpkp->bbox.zmin = zhalf;
		current->ipjpkp->bbox.zmax = current->bbox.zmax;



		// now find the points in cube ipjpkp

		for (pair<Point, unsigned int> p : current->vp) {

			if ((p.first.x > xhalf) && (p.first.y > yhalf) && (p.first.z > zhalf)) {
				current->ipjpkp->vp.push_back(p); // update vector of points for cube ipjpkp

			}


		}

		// keep particles sorted ?
		sort(current->ipjpkp->vp.begin(), current->ipjpkp->vp.end(), sorter_pair());
		it = unique(current->ipjpkp->vp.begin(), current->ipjpkp->vp.end(), PointComparison_pair);
		current->ipjpkp->vp.resize(distance(current->ipjpkp->vp.begin(), it));

		current->ipjpkp->numPoints = current->ipjpkp->vp.size();

		current->ipjpkp->divided = false;

		// ********************** cube ipjpkp ********************** //

		
		current->ijk->next = current->ipjk;
		current->ipjk->next = current->ijpk;
		current->ijpk->next = current->ipjpk;
		current->ipjpk->next = current->ijkp;
		current->ijkp->next = current->ipjkp;
		current->ipjkp->next = current->ijpkp;
		current->ijpkp->next = current->ipjpkp;
		current->ipjpkp->next = current->next;
		
		
		
		
		// now put triangles in nodes

		for (unsigned int i : current->triList) {


			// first find where the three points are

			//Point p0 = (*tri)[i].p0;;
			//Point p1 = (*tri)[i].p1;
			//Point p2 = (*tri)[i].p2;

			Point p0 = (*pp)[(*tri)[i].pid[0]].first;
			Point p1 = (*pp)[(*tri)[i].pid[1]].first;
			Point p2 = (*pp)[(*tri)[i].pid[2]].first;

			unsigned int iloc0 = 0, iloc1 = 0, iloc2 = 0; // node location of triangle point ijk=1, ijpk=2, ijpk=3, ipjpk=4, ijkp=5, ipjkp=6, ijpkp=7, ipjpkp=8

			if (p0.x <= xhalf && p0.y <= yhalf && p0.z <= zhalf) { iloc0 = 1;	current->ijk->triList.push_back(i); }
			else if (p0.x > xhalf && p0.y <= yhalf && p0.z <= zhalf) { iloc0 = 2;	current->ipjk->triList.push_back(i); }
			else if (p0.x <= xhalf && p0.y > yhalf && p0.z <= zhalf) { iloc0 = 3;	current->ijpk->triList.push_back(i); }
			else if (p0.x > xhalf && p0.y > yhalf && p0.z <= zhalf) { iloc0 = 4;	current->ipjpk->triList.push_back(i); }
			else if (p0.x <= xhalf && p0.y <= yhalf && p0.z > zhalf) { iloc0 = 5;	current->ijkp->triList.push_back(i); }
			else if (p0.x > xhalf && p0.y <= yhalf && p0.z > zhalf) { iloc0 = 6;	current->ipjkp->triList.push_back(i); }
			else if (p0.x <= xhalf && p0.y > yhalf && p0.z > zhalf) { iloc0 = 7;	current->ijpkp->triList.push_back(i); }
			else if (p0.x > xhalf && p0.y > yhalf && p0.z > zhalf) { iloc0 = 8;	current->ipjpkp->triList.push_back(i); }

			if (p1.x <= xhalf && p1.y <= yhalf && p1.z <= zhalf) { iloc1 = 1;	current->ijk->triList.push_back(i); }
			else if (p1.x > xhalf && p1.y <= yhalf && p1.z <= zhalf) { iloc1 = 2;	 current->ipjk->triList.push_back(i); }
			else if (p1.x <= xhalf && p1.y > yhalf && p1.z <= zhalf) { iloc1 = 3;	 current->ijpk->triList.push_back(i); }
			else if (p1.x > xhalf && p1.y > yhalf && p1.z <= zhalf) { iloc1 = 4;	 current->ipjpk->triList.push_back(i); }
			else if (p1.x <= xhalf && p1.y <= yhalf && p1.z > zhalf) { iloc1 = 5;	 current->ijkp->triList.push_back(i); }
			else if (p1.x > xhalf && p1.y <= yhalf && p1.z > zhalf) { iloc1 = 6;	current->ipjkp->triList.push_back(i); }
			else if (p1.x <= xhalf && p1.y > yhalf && p1.z > zhalf) { iloc1 = 7;	current->ijpkp->triList.push_back(i); }
			else if (p1.x > xhalf && p1.y > yhalf && p1.z > zhalf) { iloc1 = 8;	current->ipjpkp->triList.push_back(i); }

			if (p2.x <= xhalf && p2.y <= yhalf && p2.z <= zhalf) { iloc2 = 1;	current->ijk->triList.push_back(i); }
			else if (p2.x > xhalf && p2.y <= yhalf && p2.z <= zhalf) { iloc2 = 2;	current->ipjk->triList.push_back(i); }
			else if (p2.x <= xhalf && p2.y > yhalf && p2.z <= zhalf) { iloc2 = 3;	current->ijpk->triList.push_back(i); }
			else if (p2.x > xhalf && p2.y > yhalf && p2.z <= zhalf) { iloc2 = 4;	current->ipjpk->triList.push_back(i); }
			else if (p2.x <= xhalf && p2.y <= yhalf && p2.z > zhalf) { iloc2 = 5;	current->ijkp->triList.push_back(i); }
			else if (p2.x > xhalf && p2.y <= yhalf && p2.z > zhalf) { iloc2 = 6;	current->ipjkp->triList.push_back(i); }
			else if (p2.x <= xhalf && p2.y > yhalf && p2.z > zhalf) { iloc2 = 7;	current->ijpkp->triList.push_back(i); }
			else if (p2.x > xhalf && p2.y > yhalf && p2.z > zhalf) { iloc2 = 8;	current->ipjpkp->triList.push_back(i); }
				
				
				Vector n;
				
				n= Vector(p0, p1).normal();
					
					if ( (iloc0 == 1 && iloc1 == 4) || (iloc0 == 4 && iloc1 == 1) ){
												
							double t = ((current->ipjk->bbox.xmin)-p0.x)/n.vx;
							double y = p0.y + n.vy*t;
							double z = p0.z + n.vz*t;
							
							if (y <= yhalf && z <= zhalf) current->ipjk->triList.push_back(i);
							else current->ijpk->triList.push_back(i);
						
						}

					else if ((iloc0 == 2 && iloc1 == 3) || (iloc0 == 3 && iloc1 == 2)){

						double t = ((current->ijk->bbox.xmax) - p0.x) / n.vx;
						double y = p0.y + n.vy*t;
						double z = p0.z + n.vz*t;

						if (y <= yhalf && z <= zhalf) current->ijk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);

					}


					else if ((iloc0 == 7 && iloc1 == 6) || (iloc0 == 6 && iloc1 == 7)){

						double t = ((current->ijkp->bbox.xmax) - p0.x) / n.vx;
						double y = p0.y + n.vy*t;
						double z = p0.z + n.vz*t;

						if (y <= yhalf && z > zhalf) current->ijkp->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}


					else if ((iloc0 == 5 && iloc1 == 8) || (iloc0 == 8 && iloc1 == 5)){

						double t = ((current->ipjkp->bbox.xmin) - p0.x) / n.vx;
						double y = p0.y + n.vy*t;
						double z = p0.z + n.vz*t;

						if (y <= yhalf && z > zhalf) current->ipjkp->triList.push_back(i);
						else current->ijpkp->triList.push_back(i);

					}

					else if ((iloc0 == 1 && iloc1 == 7) || (iloc0 == 7 && iloc1 == 1)){

						double t = ((current->ijpk->bbox.ymin) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x <= xhalf && z <= zhalf) current->ijpk->triList.push_back(i);
						else current->ijkp->triList.push_back(i);

					}

					else if ((iloc0 == 3 && iloc1 == 5) || (iloc0 == 5 && iloc1 == 3)){

						double t = ((current->ijk->bbox.ymax) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x <= xhalf && z <= zhalf) current->ijk->triList.push_back(i);
						else current->ijpkp->triList.push_back(i);

					}

					else if ((iloc0 == 2 && iloc1 == 8) || (iloc0 == 8 && iloc1 == 2)){

						double t = ((current->ipjpk->bbox.ymin) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x > xhalf && z <= zhalf) current->ipjpk->triList.push_back(i);
						else current->ipjkp->triList.push_back(i);

					}

					else if ((iloc0 == 4 && iloc1 == 6) || (iloc0 == 6 && iloc1 == 4)){

						double t = ((current->ipjk->bbox.ymax) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x > xhalf && z <= zhalf) current->ipjk->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}

					else if ((iloc0 == 1 && iloc1 == 6) || (iloc0 == 6 && iloc1 == 1)){

						double t = ((current->ijkp->bbox.zmax) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y <= yhalf) current->ijkp->triList.push_back(i);
						else current->ipjk->triList.push_back(i);

					}

					else if ((iloc0 == 5 && iloc1 == 2) || (iloc0 == 2 && iloc1 == 5)){

						double t = ((current->ijk->bbox.zmax) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y <= yhalf) current->ijk->triList.push_back(i);
						else current->ipjkp->triList.push_back(i);

					}

					else if ((iloc0 == 3 && iloc1 == 8) || (iloc0 == 8 && iloc1 == 3)){

						double t = ((current->ijpkp->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpkp->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);

					}

					else if ((iloc0 == 4 && iloc1 == 7) || (iloc0 == 7 && iloc1 == 4)){

						double t = ((current->ijpk->bbox.zmax) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}

					else if ((iloc0 == 2 && iloc1 == 7) || (iloc0 == 7 && iloc1 == 2)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc0 == 3 && iloc1 == 6) || (iloc0 == 6 && iloc1 == 3)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc0 == 1 && iloc1 == 8) || (iloc0 == 8 && iloc1 == 1)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc0 == 4 && iloc1 == 5) || (iloc0 == 5 && iloc1 == 4)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

							
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////

					n = Vector(p0, p2).normal();

					if ((iloc0 == 1 && iloc2 == 4) || (iloc0 == 4 && iloc2 == 1)){

						double t = ((current->ipjk->bbox.xmin) - p0.x) / n.vx;
						double y = p0.y + n.vy*t;
						double z = p0.z + n.vz*t;

						if (y <= yhalf && z <= zhalf) current->ipjk->triList.push_back(i);
						else current->ijpk->triList.push_back(i);

					}

					else if ((iloc0 == 2 && iloc2 == 3) || (iloc0 == 3 && iloc2 == 2)){

						double t = ((current->ijk->bbox.xmax) - p0.x) / n.vx;
						double y = p0.y + n.vy*t;
						double z = p0.z + n.vz*t;

						if (y <= yhalf && z <= zhalf) current->ijk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);

					}


					else if ((iloc0 == 7 && iloc2 == 6) || (iloc0 == 6 && iloc2 == 7)){

						double t = ((current->ijkp->bbox.xmax) - p0.x) / n.vx;
						double y = p0.y + n.vy*t;
						double z = p0.z + n.vz*t;

						if (y <= yhalf && z > zhalf) current->ijkp->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}


					else if ((iloc0 == 5 && iloc2 == 8) || (iloc0 == 8 && iloc2 == 5)){

						double t = ((current->ipjkp->bbox.xmin) - p0.x) / n.vx;
						double y = p0.y + n.vy*t;
						double z = p0.z + n.vz*t;

						if (y <= yhalf && z > zhalf) current->ipjkp->triList.push_back(i);
						else current->ijpkp->triList.push_back(i);

					}

					else if ((iloc0 == 1 && iloc2 == 7) || (iloc0 == 7 && iloc2 == 1)){

						double t = ((current->ijpk->bbox.ymin) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x <= xhalf && z <= zhalf) current->ijpk->triList.push_back(i);
						else current->ijkp->triList.push_back(i);

					}

					else if ((iloc0 == 3 && iloc2 == 5) || (iloc0 == 5 && iloc2 == 3)){

						double t = ((current->ijk->bbox.ymax) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x <= xhalf && z <= zhalf) current->ijk->triList.push_back(i);
						else current->ijpkp->triList.push_back(i);

					}

					else if ((iloc0 == 2 && iloc2 == 8) || (iloc0 == 8 && iloc2 == 2)){

						double t = ((current->ipjpk->bbox.ymin) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x > xhalf && z <= zhalf) current->ipjpk->triList.push_back(i);
						else current->ipjkp->triList.push_back(i);

					}

					else if ((iloc0 == 4 && iloc2 == 6) || (iloc0 == 6 && iloc2 == 4)){

						double t = ((current->ipjk->bbox.ymax) - p0.y) / n.vy;
						double x = p0.x + n.vx*t;
						double z = p0.z + n.vz*t;

						if (x > xhalf && z <= zhalf) current->ipjk->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}

					else if ((iloc0 == 1 && iloc2 == 6) || (iloc0 == 6 && iloc2 == 1)){

						double t = ((current->ijkp->bbox.zmax) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y <= yhalf) current->ijkp->triList.push_back(i);
						else current->ipjk->triList.push_back(i);

					}

					else if ((iloc0 == 5 && iloc2 == 2) || (iloc0 == 2 && iloc2 == 5)){

						double t = ((current->ijk->bbox.zmax) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y <= yhalf) current->ijk->triList.push_back(i);
						else current->ipjkp->triList.push_back(i);

					}

					else if ((iloc0 == 3 && iloc2 == 8) || (iloc0 == 8 && iloc2 == 3)){

						double t = ((current->ijpkp->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpkp->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);

					}

					else if ((iloc0 == 4 && iloc2 == 7) || (iloc0 == 7 && iloc2 == 4)){

						double t = ((current->ijpk->bbox.zmax) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}

					else if ((iloc0 == 2 && iloc2 == 7) || (iloc0 == 7 && iloc2 == 2)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc0 == 3 && iloc2 == 6) || (iloc0 == 6 && iloc2 == 3)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc0 == 1 && iloc2 == 8) || (iloc0 == 8 && iloc2 == 1)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc0 == 4 && iloc2 == 5) || (iloc0 == 5 && iloc2 == 4)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}


					/////////////////////////////////////////////////////////////////////////////////////////////////////////////


					n = Vector(p1, p2).normal();

					if ((iloc1 == 1 && iloc2 == 4) || (iloc1 == 4 && iloc2 == 1)){

						double t = ((current->ipjk->bbox.xmin) - p1.x) / n.vx;
						double y = p1.y + n.vy*t;
						double z = p1.z + n.vz*t;

						if (y <= yhalf && z <= zhalf) current->ipjk->triList.push_back(i);
						else current->ijpk->triList.push_back(i);

					}

					else if ((iloc1 == 2 && iloc2 == 3) || (iloc1 == 3 && iloc2 == 2)){

						double t = ((current->ijk->bbox.xmax) - p1.x) / n.vx;
						double y = p1.y + n.vy*t;
						double z = p1.z + n.vz*t;

						if (y <= yhalf && z <= zhalf) current->ijk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);

					}


					else if ((iloc1 == 7 && iloc2 == 6) || (iloc1 == 6 && iloc2 == 7)){

						double t = ((current->ijkp->bbox.xmax) - p1.x) / n.vx;
						double y = p1.y + n.vy*t;
						double z = p1.z + n.vz*t;

						if (y <= yhalf && z > zhalf) current->ijkp->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}


					else if ((iloc1 == 5 && iloc2 == 8) || (iloc1 == 8 && iloc2 == 5)){

						double t = ((current->ipjkp->bbox.xmin) - p1.x) / n.vx;
						double y = p1.y + n.vy*t;
						double z = p1.z + n.vz*t;

						if (y <= yhalf && z > zhalf) current->ipjkp->triList.push_back(i);
						else current->ijpkp->triList.push_back(i);

					}

					else if ((iloc1 == 1 && iloc2 == 7) || (iloc1 == 7 && iloc2 == 1)){

						double t = ((current->ijpk->bbox.ymin) - p1.y) / n.vy;
						double x = p1.x + n.vx*t;
						double z = p1.z + n.vz*t;

						if (x <= xhalf && z <= zhalf) current->ijpk->triList.push_back(i);
						else current->ijkp->triList.push_back(i);

					}

					else if ((iloc1 == 3 && iloc2 == 5) || (iloc1 == 5 && iloc2 == 3)){

						double t = ((current->ijk->bbox.ymax) - p1.y) / n.vy;
						double x = p1.x + n.vx*t;
						double z = p1.z + n.vz*t;

						if (x <= xhalf && z <= zhalf) current->ijk->triList.push_back(i);
						else current->ijpkp->triList.push_back(i);

					}

					else if ((iloc1 == 2 && iloc2 == 8) || (iloc1 == 8 && iloc2 == 2)){

						double t = ((current->ipjpk->bbox.ymin) - p1.y) / n.vy;
						double x = p1.x + n.vx*t;
						double z = p1.z + n.vz*t;

						if (x > xhalf && z <= zhalf) current->ipjpk->triList.push_back(i);
						else current->ipjkp->triList.push_back(i);

					}

					else if ((iloc1 == 4 && iloc2 == 6) || (iloc1 == 6 && iloc2 == 4)){

						double t = ((current->ipjk->bbox.ymax) - p1.y) / n.vy;
						double x = p1.x + n.vx*t;
						double z = p1.z + n.vz*t;

						if (x > xhalf && z <= zhalf) current->ipjk->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}

					else if ((iloc1 == 1 && iloc2 == 6) || (iloc1 == 6 && iloc2 == 1)){

						double t = ((current->ijkp->bbox.zmax) - p1.z) / n.vz;
						double x = p1.x + n.vx*t;
						double y = p1.y + n.vy*t;

						if (x <= xhalf && y <= yhalf) current->ijkp->triList.push_back(i);
						else current->ipjk->triList.push_back(i);

					}

					else if ((iloc1 == 5 && iloc2 == 2) || (iloc1 == 2 && iloc2 == 5)){

						double t = ((current->ijk->bbox.zmax) - p1.z) / n.vz;
						double x = p1.x + n.vx*t;
						double y = p1.y + n.vy*t;

						if (x <= xhalf && y <= yhalf) current->ijk->triList.push_back(i);
						else current->ipjkp->triList.push_back(i);

					}

					else if ((iloc1 == 3 && iloc2 == 8) || (iloc1 == 8 && iloc2 == 3)){

						double t = ((current->ijpkp->bbox.zmin) - p1.z) / n.vz;
						double x = p1.x + n.vx*t;
						double y = p1.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpkp->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);

					}

					else if ((iloc1 == 4 && iloc2 == 7) || (iloc1 == 7 && iloc2 == 4)){

						double t = ((current->ijpk->bbox.zmax) - p1.z) / n.vz;
						double x = p1.x + n.vx*t;
						double y = p1.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpkp->triList.push_back(i);

					}

					else if ((iloc1 == 2 && iloc2 == 7) || (iloc1 == 7 && iloc2 == 2)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc1 == 3 && iloc2 == 6) || (iloc1 == 6 && iloc2 == 3)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc1 == 1 && iloc2 == 8) || (iloc1 == 8 && iloc2 == 1)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}

					else if ((iloc1 == 4 && iloc2 == 5) || (iloc1 == 5 && iloc2 == 4)){

						/*double t = ((current->ijpk->bbox.zmin) - p0.z) / n.vz;
						double x = p0.x + n.vx*t;
						double y = p0.y + n.vy*t;

						if (x <= xhalf && y > yhalf) current->ijpk->triList.push_back(i);
						else current->ipjpk->triList.push_back(i);*/

					}


					/////////////////////////////////////////////////////////////////////////////////////////////////////////////
				

			}


			



			vector<unsigned int>::iterator itt;

			//1
			sort (current->ijk->triList.begin(), current->ijk->triList.end() );
			itt = unique(current->ijk->triList.begin(), current->ijk->triList.end());
			current->ijk->triList.resize(distance(current->ijk->triList.begin(), itt));

			//2
			sort(current->ipjk->triList.begin(), current->ipjk->triList.end());
			itt = unique(current->ipjk->triList.begin(), current->ipjk->triList.end());
			current->ipjk->triList.resize(distance(current->ipjk->triList.begin(), itt));

			//3
			sort(current->ijpk->triList.begin(), current->ijpk->triList.end());
			itt = unique(current->ijpk->triList.begin(), current->ijpk->triList.end());
			current->ijpk->triList.resize(distance(current->ijpk->triList.begin(), itt));

			//4
			sort(current->ipjpk->triList.begin(), current->ipjpk->triList.end());
			itt = unique(current->ipjpk->triList.begin(), current->ipjpk->triList.end());
			current->ipjpk->triList.resize(distance(current->ipjpk->triList.begin(), itt));

			//5
			sort(current->ijkp->triList.begin(), current->ijkp->triList.end());
			itt = unique(current->ijkp->triList.begin(), current->ijkp->triList.end());
			current->ijkp->triList.resize(distance(current->ijkp->triList.begin(), itt));

			//6
			sort(current->ipjkp->triList.begin(), current->ipjkp->triList.end());
			itt = unique(current->ipjkp->triList.begin(), current->ipjkp->triList.end());
			current->ipjkp->triList.resize(distance(current->ipjkp->triList.begin(), itt));

			//7
			sort(current->ijpkp->triList.begin(), current->ijpkp->triList.end());
			itt = unique(current->ijpkp->triList.begin(), current->ijpkp->triList.end());
			current->ijpkp->triList.resize(distance(current->ijpkp->triList.begin(), itt));

			//8
			sort(current->ipjpkp->triList.begin(), current->ipjpkp->triList.end());
			itt = unique(current->ipjpkp->triList.begin(), current->ipjpkp->triList.end());
			current->ipjpkp->triList.resize(distance(current->ipjpkp->triList.begin(), itt));

			


		}


		

		

		//sort(current->ijk->vp.begin(), current->ijk->vp.end(), sorter_pair());
		//it = unique(current->ijpkp->vp.begin(), current->ijpkp->vp.end(), PointComparison_pair);
		//current->ijpkp->vp.resize(distance(current->ijpkp->vp.begin(), it));





};

///////////////////////////////////////////////////////////////////////////////////////////////

class Ray {
public:

	Point p;
	Vector n;
	vector<pair<Point, unsigned int>> CandidPoint; // list of all Candidate points for ray intersection
	vector<pair<double, unsigned int>> distance_pair; // list of intersection triangle Id's with their distance from the Ray

	vector<unsigned int> CandidTri; // list of all Candidate triangles for ray intersection
	vector<pair<double, unsigned int>> distance_tri; // list of intersection triangle Id's with their distance from the Ray

	Ray(Point p0, Vector n0) : p(p0), n(n0) {}

	Ray(Point p0, double nx, double ny, double nz) : p(p0) {
		n.vx = nx;
		n.vy = ny;
		n.vz = nz;
	}


	//bool ifIntersect(octree &O) {
	bool ifIntersect(Node *N) {
		double x, y, z;
		double t;

		// test 1: Ray intersect with plane xmin of the cell
		if (n.vx) {
			t = (N->bbox.xmin - p.x) / n.vx;
			y = p.y + n.vy*t;
			z = p.z + n.vz*t;
			if (y >= N->bbox.ymin && y < N->bbox.ymax && z >= N->bbox.zmin && z < N->bbox.zmax) return true;
		

		// test 2: Ray intersect with plane xmax of the cell
		
			t = (N->bbox.xmax - p.x) / n.vx;
			y = p.y + n.vy*t;
			z = p.z + n.vz*t;
			if (y >= N->bbox.ymin && y < N->bbox.ymax && z >= N->bbox.zmin && z < N->bbox.zmax) return true;
		}

		// test 3: Ray intersect with plane ymin of the cell
		if (n.vy) {
			t = (N->bbox.ymin - p.y) / n.vy;
			x = p.x + n.vx*t;
			z = p.z + n.vz*t;
			if (x >= N->bbox.xmin && x < N->bbox.xmax && z >= N->bbox.zmin && z < N->bbox.zmax) return true;
		

		// test 4: Ray intersect with plane ymax of the cell
		
			t = (N->bbox.ymax - p.y) / n.vy;
			x = p.x + n.vx*t;
			z = p.z + n.vz*t;
			if (x >= N->bbox.xmin && x < N->bbox.xmax && z >= N->bbox.zmin && z < N->bbox.zmax) return true;
		}

		// test 5: Ray intersect with plane zmin of the cell
		if (n.vz) {
			t = (N->bbox.zmin - p.z) / n.vz;
			x = p.x + n.vx*t;
			y = p.y + n.vy*t;
			if (x >= N->bbox.xmin && x < N->bbox.xmax && y >= N->bbox.ymin && y < N->bbox.ymax) return true;
		

		// test 6: Ray intersect with plane zmax of the cell
		
			t = (N->bbox.zmax - p.z) / n.vz;
			x = p.x + n.vx*t;
			y = p.y + n.vy*t;
			if (x >= N->bbox.xmin && x < N->bbox.xmax && y >= N->bbox.ymin && y < N->bbox.ymax) return true;
		}

		return false;
	}

	



	void find_CandidPoints(octree &O) {

		tStart = cpuSecond();

		cout << "intersecting ray & octree" << endl;

		//Node* Ntemp;

		O.current = O.Root;

		unsigned int layer = 1;

		//if (!ifIntersect(O)) return; // if the ray does not intersect with the octree Root no need for further checking

		ofstream f("Ray-pointlist.txt");


		if (!ifIntersect(O.current)) return;

		while (O.current != NULL) {



			//		if (current->numPoints > numcell && current->divided == false) {
			//			//divide(); current->divided = true;
			//			//if (current != Root) current->vp.clear();
			//		}

			//if (ifIntersect(O.current)) {
			bool result = ifIntersect(O.current);

				if (result && O.current->numPoints <= O.numcell && O.current->numPoints>0) {
					for (auto pp : O.current->vp) {
						this->CandidPoint.push_back(pp);
					}
				}

				if (result && O.current->numPoints > O.numcell) {



					if (ifIntersect(O.current->ijk) && O.current->ijk->numPoints <= O.numcell) {
						//O.current->ijk->~Node();
						for (auto pp : O.current->ijk->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ijk;

					}
					if (ifIntersect(O.current->ijk) && O.current->ijk->numPoints > O.numcell) {
						O.current = O.current->ijk;
						continue;
					}


					if (ifIntersect(O.current->ipjk) && O.current->ipjk->numPoints <= O.numcell) {
						//O.current->ipjk->~Node();
						for (auto pp : O.current->ipjk->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ipjk;

					}
					if (ifIntersect(O.current->ipjk) && O.current->ipjk->numPoints > O.numcell) {
						O.current = O.current->ipjk;
						continue;
					}


					if (ifIntersect(O.current->ijpk) && O.current->ijpk->numPoints <= O.numcell) {
						//O.current->ijpk->~Node();
						for (auto pp : O.current->ijpk->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ijpk;

					}
					if (ifIntersect(O.current->ijpk) && O.current->ijpk->numPoints > O.numcell) {
						O.current = O.current->ijpk;
						continue;
					}


					if (ifIntersect(O.current->ipjpk) && O.current->ipjpk->numPoints <= O.numcell) {
						//O.current->ipjpk->~Node();
						for (auto pp : O.current->ipjpk->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ipjpk;

					}
					if (ifIntersect(O.current->ipjpk) && O.current->ipjpk->numPoints > O.numcell) {
						O.current = O.current->ipjpk;
						continue;
					}


					if (ifIntersect(O.current->ijkp) && O.current->ijkp->numPoints <= O.numcell) {
						//O.current->ijkp->~Node();
						for (auto pp : O.current->ijkp->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ijkp;

					}
					if (ifIntersect(O.current->ijkp) && O.current->ijkp->numPoints > O.numcell) {
						O.current = O.current->ijkp;
						continue;
					}


					if (ifIntersect(O.current->ipjkp) && O.current->ipjkp->numPoints <= O.numcell) {
						//O.current->ipjkp->~Node();
						for (auto pp : O.current->ipjkp->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ipjkp;

					}
					if (ifIntersect(O.current->ipjkp) && O.current->ipjkp->numPoints > O.numcell) {
						O.current = O.current->ipjkp;
						continue;
					}


					if (ifIntersect(O.current->ijpkp) && O.current->ijpkp->numPoints <= O.numcell) {
						//O.current->ijpkp->~Node();
						for (auto pp : O.current->ijpkp->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ijpkp;
					}
					if (ifIntersect(O.current->ijpkp) && O.current->ijpkp->numPoints > O.numcell) {
						O.current = O.current->ijpkp;
						continue;
					}


					if (ifIntersect(O.current->ipjpkp) && O.current->ipjpkp->numPoints <= O.numcell) {
						//O.current->ipjpkp->~Node();
						for (auto pp : O.current->ipjpkp->vp) {
							this->CandidPoint.push_back(pp);
						}
						//O.current = O.current->ipjpkp;

					}
					if (ifIntersect(O.current->ipjpkp) && O.current->ipjpkp->numPoints > O.numcell) {
						O.current = O.current->ipjpkp;
						continue;
					}
				}

			//}

			//		Ntemp = current;
			O.current = O.current->next;
			//		Ntemp->~Node();

			layer++;
			cout << "\r layer(s) scanned= " << layer;


		}

		f.close();

		tStop = cpuSecond();

		cout << endl << "DONE!" << endl;
		cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
		cout << endl;



	}


	//////////////////////////////////////////////////////////////////////////////////////////////

	void find_CandidTriangles(octree &O) {

		//tStart = cpuSecond();

		cout << "intersecting ray & octree" << endl;

		//Node* Ntemp;

		O.current = O.Root;

		unsigned int layer = 1;

		//if (!ifIntersect(O)) return; // if the ray does not intersect with the octree Root no need for further checking

		ofstream f("Ray-pointlist.txt");


		if (!ifIntersect(O.current)) return;

		while (O.current != NULL) {



			//		if (current->numPoints > numcell && current->divided == false) {
			//			//divide(); current->divided = true;
			//			//if (current != Root) current->vp.clear();
			//		}

			//if (ifIntersect(O.current)) {
			bool result = ifIntersect(O.current);

			if (result && O.current->numPoints <= O.numcell && O.current->numPoints>0) {
				for (auto t : O.current->triList) {
					this->CandidTri.push_back(t);
				}
			}

			else if (result && O.current->numPoints > O.numcell) {



				if (ifIntersect(O.current->ijk) && O.current->ijk->numPoints <= O.numcell) {
					//O.current->ijk->~Node();
					for (auto t : O.current->ijk->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ijk;

				}
				else if (ifIntersect(O.current->ijk) && O.current->ijk->numPoints > O.numcell) {
					O.current = O.current->ijk;
					continue;
				}


				if (ifIntersect(O.current->ipjk) && O.current->ipjk->numPoints <= O.numcell) {
					//O.current->ipjk->~Node();
					for (auto t : O.current->ipjk->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ipjk;

				}
				else if (ifIntersect(O.current->ipjk) && O.current->ipjk->numPoints > O.numcell) {
					O.current = O.current->ipjk;
					continue;
				}


				if (ifIntersect(O.current->ijpk) && O.current->ijpk->numPoints <= O.numcell) {
					//O.current->ijpk->~Node();
					for (auto t : O.current->ijpk->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ijpk;

				}
				else if (ifIntersect(O.current->ijpk) && O.current->ijpk->numPoints > O.numcell) {
					O.current = O.current->ijpk;
					continue;
				}


				if (ifIntersect(O.current->ipjpk) && O.current->ipjpk->numPoints <= O.numcell) {
					//O.current->ipjpk->~Node();
					for (auto t : O.current->ipjpk->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ipjpk;

				}
				else if (ifIntersect(O.current->ipjpk) && O.current->ipjpk->numPoints > O.numcell) {
					O.current = O.current->ipjpk;
					continue;
				}


				if (ifIntersect(O.current->ijkp) && O.current->ijkp->numPoints <= O.numcell) {
					//O.current->ijkp->~Node();
					for (auto t : O.current->ijkp->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ijkp;

				}
				else if (ifIntersect(O.current->ijkp) && O.current->ijkp->numPoints > O.numcell) {
					O.current = O.current->ijkp;
					continue;
				}


				if (ifIntersect(O.current->ipjkp) && O.current->ipjkp->numPoints <= O.numcell) {
					//O.current->ipjkp->~Node();
					for (auto t : O.current->ipjkp->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ipjkp;

				}
				else if (ifIntersect(O.current->ipjkp) && O.current->ipjkp->numPoints > O.numcell) {
					O.current = O.current->ipjkp;
					continue;
				}


				if (ifIntersect(O.current->ijpkp) && O.current->ijpkp->numPoints <= O.numcell) {
					//O.current->ijpkp->~Node();
					for (auto t : O.current->ijpkp->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ijpkp;
				}
				else if (ifIntersect(O.current->ijpkp) && O.current->ijpkp->numPoints > O.numcell) {
					O.current = O.current->ijpkp;
					continue;
				}


				if (ifIntersect(O.current->ipjpkp) && O.current->ipjpkp->numPoints <= O.numcell) {
					//O.current->ipjpkp->~Node();
					for (auto t : O.current->ipjpkp->triList) {
						this->CandidTri.push_back(t);
					}
					//O.current = O.current->ipjpkp;

				}
				else if (ifIntersect(O.current->ipjpkp) && O.current->ipjpkp->numPoints > O.numcell) {
					O.current = O.current->ipjpkp;
					continue;
				}
			}

			//}

			//		Ntemp = current;
			O.current = O.current->next;
			//		Ntemp->~Node();

			layer++;
			cout << "\r layer(s) scanned= " << layer;


		}

		f.close();

		//tStop = cpuSecond();

		cout << endl << "DONE!" << endl;
		//cout << "----->   elapsed time: " << double(tStop - tStart) << " [milisec]" << endl;
		//cout << endl;



	}

	/////////////////////////////////////////////////////////////////////////////////////////////

	unsigned int find_Intersection(octree &O) {

		// this function returns the Id of the triangle that the ray intersects with

		for (auto pp : this->CandidPoint) {

			// for each point in the CandidPoint list, find each Triangle its connecting to
			// and perform intersection test
			for (unsigned int numTri = 0; numTri < (*O.pt)[pp.second].size(); numTri++) {


				// find intersection point of the Ray with the plane defined by Triangle
				Triangle TRI((*O.tri)[numTri]);

				double nx = this->n.vx;
				double ny = this->n.vy;
				double nz = this->n.vz;

				double num = TRI.nx*(TRI.p0.x - this->p.x) + TRI.ny*(TRI.p0.y - this->p.y) + TRI.nz*(TRI.p0.z - this->p.z);
				double denum = TRI.nx*nx + TRI.ny*ny + TRI.nz*nz;

				double t = 0;
				if (denum != 0) {
					t = num / denum;
				}
				// if (denum==0.0) then what



				double xpp = nx*t + this->p.x;
				double ypp = ny*t + this->p.y;
				double zpp = nz*t + this->p.z;

				Point OO(xpp, ypp, zpp);

				// now check if this point is inside the triangle

				Triangle OP0P1(OO, TRI.p0, TRI.p1);
				Triangle OP0P2(OO, TRI.p0, TRI.p2);
				Triangle OP1P2(OO, TRI.p1, TRI.p2);

				if (double d = abs((OP0P1.area2() + OP0P2.area2() + OP1P2.area2()) - TRI.area2()) < 1e-5) {
					distance_pair.push_back(pair<double, unsigned int>(d, numTri));
				}


			}

		}

		// now find the minimum distance to find the final intersection

		double min_distance = 0.0;
		unsigned int return_id;

		for (unsigned int i = 0; i < distance_pair.size(); i++) {
			if (min_distance > distance_pair[i].first) {
				min_distance = distance_pair[i].first;
				return_id = distance_pair[i].second;
			}

			i++;
		}

		return 0;
	}









	/////////////////////////////////////////////////////////////////////////////////////////////

	

	/////////////////////////////////////////////////////////////////////////////////////////////

	unsigned int find_Tri_Intersection(octree &O) {

		// this function returns the Id of the triangle that the ray intersects with

		for (auto tt : this->CandidTri) {
						

				// find intersection point of the Ray with the plane defined by Triangle
				Triangle TRI((*O.tri)[tt]);

				double nx = this->n.vx;
				double ny = this->n.vy;
				double nz = this->n.vz;

				double num = TRI.nx*(TRI.p0.x - this->p.x) + TRI.ny*(TRI.p0.y - this->p.y) + TRI.nz*(TRI.p0.z - this->p.z);
				double denum = TRI.nx*nx + TRI.ny*ny + TRI.nz*nz;

				double t = 0;
				if (denum != 0) {
					t = num / denum;
				}
				// if (denum==0.0) then what



				double xpp = nx*t + this->p.x;
				double ypp = ny*t + this->p.y;
				double zpp = nz*t + this->p.z;

				Point OO(xpp, ypp, zpp);

				// now check if this point is inside the triangle

				Triangle OP0P1(OO, TRI.p0, TRI.p1);
				Triangle OP0P2(OO, TRI.p0, TRI.p2);
				Triangle OP1P2(OO, TRI.p1, TRI.p2);

				// original
				//if (double d = abs((OP0P1.area2() + OP0P2.area2() + OP1P2.area2()) - TRI.area2()) < 1e-5) {
				//	distance_tri.push_back(pair<double, unsigned int>(d, tt));
				//}


				// now check if this point is inside the triangle

				//Triangle OP0P1(O, D.t[i].p0, D.t[i].p1);
				//Triangle OP0P2(O, D.t[i].p0, D.t[i].p2);
				//Triangle OP1P2(O, D.t[i].p1, D.t[i].p2);

				//double d;
				//if (d = abs((OP0P1.area2() + OP0P2.area2() + OP1P2.area2()) - D.t[i].area2()) < 1e-5) {
				//distance_pair.push_back(pair<double, unsigned int>(d, numTri));
				//d = ((OP0P1.area() + OP0P2.area() + OP1P2.area()) - D.t[i].area());

				Vector v0(TRI.p0);
				Vector v1(TRI.p1);
				Vector v2(TRI.p2);

				Vector vp(xpp, ypp, zpp);

				Vector e0 = v2 - v1;
				Vector e1 = v0 - v2;
				Vector e2 = v1 - v0;

				Vector d0 = vp - v0;
				Vector d1 = vp - v1;
				Vector d2 = vp - v2;

				Vector np(TRI.nx, TRI.ny, TRI.nz);

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
				if ((b0 >= 0.0 && b0 <= 1.0) && (b1 >= 0.0 && b1 <= 1.0) && (b2 >= 0.0 && b2 <= 1.0)) {

					double distance = sqrt((xpp - p.x)*(xpp - p.x) + (ypp - p.y)*(ypp - p.y) + (zpp - p.z)*(zpp - p.z));

					//CandidTri.push_back(pair<double, unsigned int>(distance, i));
					distance_tri.push_back(pair<double, unsigned int>(distance, tt));

					//cout << "found one! : id= " << i << "   " << d << endl;
					cout << "found one! : id= " << tt << "   " << b0 + b1 + b2 - 1.0 << endl;
				}

			

		}

		// now find the minimum distance to find the final intersection

		double min_distance = DBL_MAX;
		unsigned int return_id = UINT_MAX;;

		for (unsigned int i = 0; i < distance_tri.size(); i++) {
			if (min_distance > distance_tri[i].first) {
				min_distance = distance_tri[i].first;
				return_id = distance_tri[i].second;
			}

			
		}

		return return_id;
	}




};


class TriangleNode : public Triangle{
//
//public:
//
//	TriangleNode* parent = NULL;
//	TriangleNode *first = NULL;
//	TriangleNode *second = NULL;
//	TriangleNode *third = NULL;
//	TriangleNode* next = NULL;
//	Triangle t;
//
//	TriangleNode(Triangle tt) : t(tt) {}
//	
//
//
//};
//
//class TriangleTree {
//
//public:
//	
//	TriangleNode* Root;
//	TriangleNode* current;
//	
//	TriangleTree(Triangle t, double landingpoint_radius) {
//		Root = new TriangleNode(t);
//		Root->parent = NULL;
//		Root->next = NULL;
//		//Root->first = &(t.neighbor[0]);
//
//	}


};