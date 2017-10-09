#pragma once
#include "geometry.h"


struct Median{
	double val; // Median value in the original unsorted vector<Point>
	unsigned int id; // id in the original unsorted vector<Point>

	Median(double v, unsigned int i) : val(v), id(i) {}
	Median() {}
};



class Node {

public:

	Node* parent;
	Node* rhs;
	Node* lhs;

	size_t numPoints; // number of points enclosed in this region
	vector<Point> vp; // a vector that holds a vector of Points enclosed in this region

	vector<vector<unsigned int>> triId; // a vector that holds id of Triangle enclosed in this region

	BoundingBox bbox; // bounding box of enclosing region
	
	char splitplane; // =0, 1, 2 (split by x-plane, split by y-plane, split by z-plane)

};

class kdtree {

public:
	Node* Root;
	Node* current;

	kdtree(vector<Point>& P, vector<Triangle>& T, vector<vector<unsigned int>>& V, BoundingBox B) {
		
		Root = new Node;
		Root->parent = NULL;
		Root->numPoints = NUM_POINTS;
		Root->bbox = B;
		Root->splitplane = 0;
		Root->triId = V;
		Root->vp = P;
		
		current = Root;

		double xMedian, yMedian, zMedian;
		unsigned int id;
		Median m;
		

		if (current->numPoints > 1) { // must be changed to while(...) ?

			if (current->splitplane == 0) m=findMedianX();
			//if (current->splitplane == 0) yMedian=findMedianX();
			//if (current->splitplane == 0) zMedian=findMedianX();


		}
	}

	kdtree(pt_data Data) {
		kdtree(Data.p, Data.t, Data.v, Data.b);
	}

	Median findMedianX() {
		
		vector<Point> vp_xsorted = current->vp;
		sort(vp_xsorted.begin(), vp_xsorted.end(), sort_by_x());
				
		unsigned int i = 0; // holds number of elements with the Median value
		unsigned int j = 0; // finds Id first element with the Median value
		unsigned int k = 0;
		for (Point p : vp_xsorted) {
			if (p.x == vp_xsorted[vp_xsorted.size() / 2].x) {
				i++;
				if (i == 1) j = k;
			}
			k++;
		}
			
		return Median(vp_xsorted[vp_xsorted.size()/2].x, j);
	}
	//Median findMedianY() {
	//	vector<Point> vp_ysorted = current->vp;
	//	sort(vp_ysorted.begin(), vp_ysorted.end(), sort_by_y());

	//	return vp_ysorted[vp_ysorted.size() / 2].y;
	//}
	//Median findMedianZ() {
	//	vector<Point> vp_zsorted = current->vp;
	//	sort(vp_zsorted.begin(), vp_zsorted.end(), sort_by_z());

	//	return vp_zsorted[vp_zsorted.size() / 2].z;
	//}

	BoundingBox findBoundingBox() {
		return BoundingBox();
	}
};
