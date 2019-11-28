#include "Delaunay.h"

using namespace std;

//�жϵ��Ƿ������������Բ��
//p(px,py) �� (x1,y1),(x2,y2),(x3,y3)��ɵ������ε����Բ�� ����true
//���Բ��Ϊ(xc,yc)���뾶Ϊr
//�����Բ�߽�����Բ��
bool CircumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3,
	double &xc, double &yc, double &r) {
	double m1, m2, mx1, mx2, my1, my2;
	double dx, dy, rsqr, drsqr;

	//�ж������
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


//������pxyz�еĶ���nv��Ϊ����
//��������v�е�ntri�����ε��б�
//��Щ��������һ�µ�˳ʱ��˳������
//����������v����3 * nv
//��������pxyz�����㹻�󣬿�������3�����ϵĵ�
//�������鰴xֵ��������   qsort(p,nv,sizeof(XYZ),XYZCompare);
int Triangulate(int nv, XYZ pxyz[], ITRIANGLE v[], int &ntri) {
	IEDGE *p_EdgeTemp;
	int nedge = 0;
	int trimax, emax = 200;
	int status = 0;

	bool inside;
	//int i, j, k;
	double xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;

	//�����ڴ棬ÿ��������һ����־λ
	trimax = 4 * nv;
	bool *complete = new bool[trimax];
	//�����ڴ�������εı�
	IEDGE *edges = new IEDGE[emax];

	//���������С����߽�,���㳬��������
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


	//��ʼ�����������Σ��������еĵ�
	//����������������ӵ������б�vertex listĩβ
	//������������triangle list�ĵ�һ��������
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

	//���������������ڲ��ĵ��������
	for (int i = 0; i < nv; ++i) {
		xp = pxyz[i].x;
		yp = pxyz[i].y;
		nedge = 0;


		//��ʼ���߻���
		//�����(xp,yp)�����Բ�ڲ�����������ε������ߴ���edge buffer����ɾ����������
		for (int j = 0; j < ntri; ++j) {
			if(complete[j]) continue;
			x1 = pxyz[v[j].pt1].x;
			y1 = pxyz[v[j].pt1].y;
			x2 = pxyz[v[j].pt2].x;
			y2 = pxyz[v[j].pt2].y;
			x3 = pxyz[v[j].pt3].x;
			y3 = pxyz[v[j].pt3].y;
			inside = CircumCircle(xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r);

			//��� xc+r+EPSILON<xp
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


		//������������ζ�����ʱ�룬�����б߶����Ϊ������
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
		//Ϊ��ǰ���γ��µ��������������б�ǹ��ı�
		//���б߶���ʱ������
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

	//ɾ���г������ζ��㲢�Ҷ���������nv��������
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