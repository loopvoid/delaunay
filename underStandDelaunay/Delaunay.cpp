#include "Delaunay.h"

using namespace std;

//判断点是否在三角形外接圆内
//p(px,py) 在 (x1,y1),(x2,y2),(x3,y3)组成的三角形的外接圆内 返回true
//外接圆心为(xc,yc)，半径为r
//在外接圆边界上算圆内
bool CircumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3,
	double &xc, double &yc, double &r) {
	double m1, m2, mx1, mx2, my1, my2;
	double dx, dy, rsqr, drsqr;

	//判断输入点
	if (abs(y1 - y2) < EPSILON && abs(y2 - y3) < EPSILON)
		return false;

	if (abs(y2 - y1) < EPSILON) {
		m2 = -(x3 - x2) / (y3 - y2);
		mx2 = (x2 + x3) / 2.0;
		my2 = (y2 + y3) / 2.0;
		xc = (x2 + x1) / 2.0;
		yc = m2 * (xc - mx2) + my2;
	}
	else if (abs(y3 - y2) < EPSILON) {
		m1 = -(x2 - x1) / (y2 - y1);
		mx1 = (x1 + x2) / 2.0;
		my1 = (y1 + y2) / 2.0;
		xc = (x3 + x2) / 2.0;
		yc = m1 * (xc - mx1) + my1;
	}
	else {
		m1 = -(x2 - x1) / (y2 - y1);
		m2 = -(x3 - x2) / (y3 - y2);
		mx1 = (x1 + x2) / 2.0;
		mx2 = (x2 + x3) / 2.0;
		my1 = (y1 + y2) / 2.0;
		my2 = (y2 + y3) / 2.0;
		xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		yc = m1 * (xc - mx1) + my1;
	}

	dx = x2 - xc;
	dy = y2 - yc;
	rsqr = dx * dx + dy * dy;
	r = sqrt(rsqr);
	dx = xp - xc;
	dy = yp - yc;
	drsqr = dx * dx + dy * dy;

	return((drsqr <= rsqr) ? true : false);
}


//将数组pxyz中的顶点nv作为输入
//返回数组v中的ntri三角形的列表
//这些三角形以一致的顺时针顺序排列
//三角形数组v分配3 * nv
//顶点数组pxyz必须足够大，可以容纳3个以上的点
//顶点数组按x值递增排序   qsort(p,nv,sizeof(XYZ),XYZCompare);
int Triangulate(int nv, XYZ pxyz[], ITRIANGLE v[], int &ntri) {
	IEDGE *p_EdgeTemp;
	int nedge = 0;
	int trimax, emax = 200;
	int status = 0;

	bool inside;
	//int i, j, k;
	double xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;

	//分配内存，每个三角形一个标志位
	trimax = 4 * nv;
	bool *complete = new bool[trimax];
	//分配内存给三角形的边
	IEDGE *edges = new IEDGE[emax];

	//计算最大最小顶点边界,计算超级三角形
	double xmin = pxyz[0].x;
	double ymin = pxyz[0].y;
	double xmax = xmin;
	double ymax = ymin;
	for (int i = 0; i < nv; ++i) {
		if (pxyz[i].x < xmin) xmin = pxyz[i].x;
		if (pxyz[i].x > xmax) xmax = pxyz[i].x;
		if (pxyz[i].y < ymin) ymin = pxyz[i].y;
		if (pxyz[i].y > ymax) ymax = pxyz[i].y;
	}
	double dx = xmax - xmin;
	double dy = ymax - ymin;
	double dmax = (dx > dy) ? dx : dy;

	double xmid = (xmax + xmin) / 2.0;
	double ymid = (ymax + ymin) / 2.0;


	//初始化超级三角形，包含所有的点
	//超级三角形坐标添加到顶点列表vertex list末尾
	//超级三角形是triangle list的第一个三角形
	pxyz[nv + 0].x = xmid - 20 * dmax;
	pxyz[nv + 0].y = ymid - dmax;
	pxyz[nv + 0].z = 0.0;
	pxyz[nv + 1].x = xmid;
	pxyz[nv + 1].y = ymid + 20 * dmax;
	pxyz[nv + 1].z = 0.0;
	pxyz[nv + 2].x = xmid + 20 * dmax;
	pxyz[nv + 2].y = ymid - dmax;
	pxyz[nv + 0].z = 0.0;
	v[0].pt1 = nv;
	v[0].pt2 = nv + 1;
	v[0].pt3 = nv + 2;
	complete[0] = false;
	ntri = 1;

	//将在三角形网格内部的点包括进来
	for (int i = 0; i < nv; ++i) {
		xp = pxyz[i].x;
		yp = pxyz[i].y;
		nedge = 0;


		//初始化边缓存
		//如果点(xp,yp)在外接圆内部，则该三角形的三条边存入edge buffer，并删除该三角形
		for (int j = 0; j < ntri; ++j) {
			if(complete[j]) continue;
			x1 = pxyz[v[j].pt1].x;
			y1 = pxyz[v[j].pt1].y;
			x2 = pxyz[v[j].pt2].x;
			y2 = pxyz[v[j].pt2].y;
			x3 = pxyz[v[j].pt3].x;
			y3 = pxyz[v[j].pt3].y;
			inside = CircumCircle(xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r);

			//如果 xc+r+EPSILON<xp
			if (xc < xp && ((xp - xc)*(xp - xc)) > r)
				complete[j] = true;

			if (inside) {
				if (nedge + 3 >= emax) {
					emax += 100;
					p_EdgeTemp = new IEDGE[emax];
					for (int i = 0; i < nedge; i++) { // Fix by John Bowman
						p_EdgeTemp[i] = edges[i];
					}
					delete[]edges;
					edges = p_EdgeTemp;
				}

				edges[nedge + 0].pt1 = v[j].pt1;
				edges[nedge + 0].pt2 = v[j].pt2;
				edges[nedge + 1].pt1 = v[j].pt2;
				edges[nedge + 1].pt2 = v[j].pt3;
				edges[nedge + 2].pt1 = v[j].pt3;
				edges[nedge + 2].pt2 = v[j].pt1;
				nedge += 3;
				v[j] = v[ntri - 1];
				complete[j] = complete[ntri - 1];
				ntri--;
				j--;
			}
		}


		//如果所有三角形都是逆时针，则所有边都标记为反方向
		for (int j = 0; j < nedge - 1; j++) {
			for (int k = j + 1; k < nedge; k++) {
				if ((edges[j].pt1 == edges[k].pt2) && (edges[j].pt2 == edges[k].pt1)) {
					edges[j].pt1 = -1;
					edges[j].pt2 = -1;
					edges[k].pt1 = -1;
					edges[k].pt2 = -1;
				}
				if ((edges[j].pt1 == edges[k].pt1) && (edges[j].pt2 == edges[k].pt2)) {
					edges[j].pt1 = -1;
					edges[j].pt2 = -1;
					edges[k].pt1 = -1;
					edges[k].pt2 = -1;
				}
			}
		}
		//为当前点形成新的三角形跳过所有标记过的边
		//所有边都逆时针排序
		for (int j = 0; j < nedge; j++) {
			if (edges[j].pt1 < 0 || edges[j].pt2 < 0)
				continue;

			if (ntri >= trimax) {
				break;
			}

			v[ntri].pt1 = edges[j].pt1;
			v[ntri].pt2 = edges[j].pt2;
			v[ntri].pt3 = i;
			complete[ntri] = false;
			ntri++;
		}

		if (ntri >= trimax) {
			break;
		}
	}

	//删除有超三角形顶点并且顶点数大于nv的三角形
	for (int i = 0; i < ntri; i++) {
		if (v[i].pt1 >= nv || v[i].pt2 >= nv || v[i].pt3 >= nv) {
			v[i] = v[ntri - 1];
			ntri--;
			i--;
		}
	}

	delete[] edges;
	delete[] complete;
	return 0;
}

int XYZCompare(const void *v1, const void *v2) {
	XYZ *pt1 = (XYZ*)v1;
	XYZ *pt2 = (XYZ*)v2;

	if (pt1->x < pt2->x) return -1;
	else if (pt1->x > pt2->x)
		return 1;
	else
		return 0;
}