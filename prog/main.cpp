//---------------------------------------------------------------------------
//#include <vcl.h>
//#pragma hdrstop

//#include <tchar.h>
//#include <conio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>


using namespace std;

#include "matrix3.h"
#include "poly.h"
#include "tensor.h"
#include "tensor3.h"
#include "vector3.h"
#include "util.h"
#include "types.h"
#include "pov_ray.h"

#include <fstream>


//---------------------------------------------------------------------------

void work(){
	Tensor t=make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);
	Matrix3 cr=t.crist(Vector3(1/sqrt(2.),1/sqrt(2.),0));
	cout << cr << endl;
	cout << t << endl;

	Poly pol = cr.characteristic();

	cout << pol << endl;

	cout << pol.root() << endl;

	cout << pol.all_roots() << endl;

	cout << rho2v(pol.all_roots(),5.96E3) << endl;


	//return;


	Matrix3 matx = rotx (0.3);
	cout << matx << endl;

	cout << endl;

	Matrix3 maty = roty (0.3);
	cout << maty << endl;

	cout << endl;

	Matrix3 matz = rotz (0.3);
	cout << matz << endl;

	cout << endl;

	//cout << matx*maty*matz << endl;
	cout << euler(0.3, 0.3, 0.3) << endl;

}

void test_Poly(){
	Poly a(1, 5, -4, -20);
	cout << a << endl;
	cout << a.all_roots() << endl;
}

Vector3 get_slow(Vector3 rn, const Tensor& t, double rho) {
    Matrix3 cr=t.crist(rn);
    Poly pol = cr.characteristic();
    VCD roots = pol.all_roots();
    VD vels = rho2v(pol.all_roots(), rho);
    sort(vels.begin(),vels.end());
    rn*=1000/vels[2];
    return rn;
}

void test_PR() {
	PovRay pov("tst.pov");
	int n=200;

	Tensor tt = make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);
	Tensor t = tt * euler( M_PI/2, M_PI/3, 0);

	for (int p = 0; p <=n; p++) {
		cout << p << endl;
		for (int q = 0; q < n; q++) {
			Vector3 rn = pq2vec(p,q,n);
			Matrix3 cr=t.crist(rn);
			Poly pol = cr.characteristic();
			VCD roots = pol.all_roots();
			VD vels = rho2v(pol.all_roots(),5.96E3);
			sort(vels.begin(),vels.end());
			rn*=1000/vels[0];
			pov.point(rn);
		}
	}
}

void test_PR_te() {
    cout << "test_PR_te" << endl;
    PovRay pov("te_PR_group_norm_2.pov");
    int n=10;
    int ind=2;
    double rad=0.04, radc=0.005;

    Tensor t = make_trigonal_tensor(3.27E10,0.86E10,2.49E10,7.22E10,1.24E10,3.12E10,1.21E10);
    Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
    Tensor3 pc = make_trigonal_piezo_tensor_32(0.42*0, 0.17*0);
    double rho = 6.25E3;

    for (int p = 0; p <=n; p++) {
        cout << p << endl;
        for (int q = 0; q < 2*n; q++) {
            Vector3 rn = pq2vec(p,q,n);

            Matrix3 cr=t.make_piezo_crist(rn, epsilon, pc);
            Poly pol = cr.characteristic();
            VCD roots = pol.all_roots();
            VD vels = rho2v(pol.all_roots(),rho);
            sort(vels.begin(),vels.end());
            rn*=1000/vels[ind];

            Vector3 drn(rn + slow_normal(rn, ind, epsilon, t, pc)*3*rad);
            Vector3 drnn(rn + rn*3*rad);

            pov.cylinder(rn, drn, radc, Vector3(0,0.9,0.1));
            pov.cylinder(rn, drnn, radc, Vector3(0,0.3,0.7));
            //pov.point(rn, rad);

        }
    }

    int nd=100;

    pov.begin_mesh();

    for (int q = 0; q < 2*nd; q++) {
        Vector3 r1 = pq2vec(0,q,nd);
        Vector3 r2 = pq2vec(1,q,nd);
        Vector3 r3 = pq2vec(1,q+1,nd);
        pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
    }

    for (int p = 1; p < nd; p++) {
        for (int q = 0; q < 2*nd; q++) {
            Vector3 r1 = pq2vec(p,q,nd);
            Vector3 r2 = pq2vec(p+1,q,nd);
            Vector3 r3 = pq2vec(p+1,q+1,nd);
            pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
        }
    }

    for (int p = 1; p < nd; p++) {
        for (int q = 0; q < 2*nd; q++) {
            Vector3 r1 = pq2vec(p,q,nd);
            Vector3 r2 = pq2vec(p,q+1,nd);
            Vector3 r3 = pq2vec(p+1,q+1,nd);
            pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
        }
    }

    for (int q = 0; q < 2*nd; q++) {
        Vector3 r1 = pq2vec(nd-1,q,nd);
        Vector3 r2 = pq2vec(nd,q,nd);
        Vector3 r3 = pq2vec(nd-1,q+1,nd);
        pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
    }

    pov.end_mesh();

}

void group() {

    ofstream dest("group.txt");

    double phi=0*M_PI/180;

    for (double theta = 0; theta < M_PI/3; theta+=0.001) {
        Vector3 n(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

        //int i=40;
        int ind=2;

        Tensor t = make_trigonal_tensor(3.27E10,0.86E10,2.49E10,7.22E10,1.24E10,3.12E10,1.21E10);

        Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
        Tensor3 pc = make_trigonal_piezo_tensor_32(0.42, 0.17);

        //for (int p = 0; p <=i; p++) {
            //cout << p << endl;
            //for (int q = 0; q < 2*i; q++) {
                //Vector3 rn = pq2vec(p,q,i);

                Matrix3 cr=t.make_piezo_crist(n, epsilon, pc);
                Poly pol = cr.characteristic();
                VCD roots = pol.all_roots();
                VD vels = rho2v(pol.all_roots(),6.25E3);
                sort(vels.begin(),vels.end());
                n*=1000/vels[ind];

                Vector3 Vgr = slow_normal(n, ind, epsilon, t, pc);

                double psi = acos(n*Vgr);
                dest << theta*180/M_PI << " " << psi*180/M_PI << endl;

    }
}

void test_beans() {
	Tensor tt = make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);

    int n_phi=1000;
    int n_theta=10000;
	double dtheta=2*M_PI/n_theta;
	double dphi=2*M_PI/n_phi;
	double v=1/sqrt(2.);

	ofstream dest("a.txt");

	PovRay pic("bean.pov");

	for (double theta = 0; theta < 2*M_PI; theta+=dtheta)
	//double theta=M_PI/2;
	{
		Vector3 n1(1, 0, 0);
		Vector3 n2(0, sin(theta), cos(theta));
		cout << "theta=" << theta << endl;
		//dest << "theta=" << theta << endl << "\t";

		bool flag = true;
		double vp=0, vc=0;
		for (double phi = 0; phi < 2*M_PI; phi+=dphi) {
			Vector3 n=n1*cos(phi)+n2*sin(phi);
			Matrix3 cr=tt.crist(n);
			Poly pol = cr.characteristic();
			VCD roots = pol.all_roots();
			VD gamms = arr_re(roots);
			sort(gamms.begin(),gamms.end());
			Vector3 q = get_polarization( cr, gamms[2]);
			//cout << " q= " << q << endl;
			vp=vc;
			vc=abs(q*n);


			if (flag) {
				flag = false;
			} else {
				if ((v-vp)*(v-vc)<0) {
                    //Vector3 nn=n*(1000/sqrt(gamms[2]/5.96e3));
					dest << theta*180/M_PI << " " << phi*180/M_PI << endl;
                    pic.cylinder(Vector3 (0,0,0), n, 0.005);

				}
			}
		}
	}

}

void test_insect() {
	Tensor tt = make_trigonal_tensor(3.27E10,0.86E10,2.49E10,7.22E10,1.24E10,3.12E10,1.21E10);

	int n_phi=1000;
	int n_theta=10000;
	double dtheta=2*M_PI/n_theta;
	double dphi=2*M_PI/n_phi;
	double v=1/sqrt(2.);

    ofstream dest("te_insect.txt");

    //PovRay pic("slow_insect.pov");

	for (double theta = 0; theta < 2*M_PI; theta+=dtheta)
    {
		Vector3 n1(1, 0, 0);
		Vector3 n2(0, sin(theta), cos(theta));
		cout << "theta=" << theta << endl;

		bool flag = true;
		double vp=0, vc=0;
		for (double phi = 0; phi < 2*M_PI; phi+=dphi) {
			Vector3 n=n1*cos(phi)+n2*sin(phi);
			Matrix3 cr=tt.crist(n);
			Poly pol = cr.characteristic();
			VCD roots = pol.all_roots();
			VD gamms = arr_re(roots);
			sort(gamms.begin(),gamms.end());
			Vector3 q = get_polarization( cr, gamms[2]);

			vp=vc;
			vc=abs(q*n);

            if (flag) {
				flag = false;
			} else {
				if ((v-vp)*(v-vc)<0) {
					Vector3 nn=n*(1000/sqrt(gamms[2]/5.96e3));
                    dest << theta*180/M_PI << " " << phi*180/M_PI << endl;
                    //pic.point(nn);
				}
			}
		}
    }

}

void test_insect_piezo() {
    //Tensor tt = make_trigonal_tensor(3.27E10,0.86E10,2.49E10,7.22E10,1.24E10,3.12E10,1.21E10);

    Tensor t = make_trigonal_tensor(3.257E10, 0.845E10, 2.57E10, 7.17E10, 1.238E10, 3.094E10, 1.206E10);

    ofstream dest("te_piezo_insect.txt");

    //double theta=0;

    //for (double phi = 0; phi < 2*M_PI; phi+=0.01) {

    //Vector3 n(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

    //Vector3 n(0,cos(phi),sin(phi));

    Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
    Tensor3 pc = make_trigonal_piezo_tensor_32(0.42, 0.17);

    //Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);

    int n_phi=10000;
    int n_theta=100000;
    double dtheta=2*M_PI/n_theta;
    double dphi=2*M_PI/n_phi;
    double v=1/sqrt(2.);

    //ofstream dest("b.txt");

    PovRay pic("slow_piezo_insect.pov");

    for (double theta = 0; theta < 2*M_PI; theta+=dtheta)

    {
        Vector3 n1(1, 0, 0);
        Vector3 n2(0, sin(theta), cos(theta));
        cout << theta << endl;
        //dest << "theta=" << theta << endl << "\t";

        bool flag = true;
        double vp=0, vc=0;
        for (double phi = 0; phi < 2*M_PI; phi+=dphi) {
            Vector3 n=n1*cos(phi)+n2*sin(phi);

            Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);

            //Matrix3 cr=tt.crist(n);
            Poly pol = tt.characteristic();
            VCD roots = pol.all_roots();
            VD gamms = arr_re(roots);
            sort(gamms.begin(),gamms.end());
            Vector3 q = get_polarization( tt, gamms[2]);
            //cout << " q= " << q << endl;
            vp=vc;
            vc=abs(q*n);


            //dest << vc << " ";

            if (flag) {
                flag = false;
            } else {
                if ((v-vp)*(v-vc)<0) {
                    Vector3 nn=n; //*(1000/sqrt(gamms[2]/5.96e3));
                    //dest << theta*180/M_PI << " " << phi*180/M_PI << endl;
                    pic.point(nn);
                }
            }
        }
        //dest << endl;
    }
}

double deviant(double phi, const Vector3& n1, const Vector3& n2, const Tensor& t, const Matrix3& epsilon, const Tensor3& pc) {

    double v=1/sqrt(2.);

    Vector3 n=n1*cos(phi)+n2*sin(phi);

    Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);

    Poly pol = tt.characteristic();
    VCD roots = pol.all_roots();
    VD gamms = arr_re(roots);
    sort(gamms.begin(),gamms.end());
    Vector3 q = get_polarization( tt, gamms[2]);

    return abs(q*n)-v;

}

void test_insect_piezo_dihotom() {

    Tensor t = make_trigonal_tensor(3.257E10, 0.845E10, 2.57E10, 7.17E10, 1.238E10, 3.094E10, 1.206E10);

    ofstream dest("te_piezo_insect_dihotom.txt");

    //double theta=0;

    //for (double phi = 0; phi < 2*M_PI; phi+=0.01) {

    //Vector3 n(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

    //Vector3 n(0,cos(phi),sin(phi));

    Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
    Tensor3 pc = make_trigonal_piezo_tensor_32(0.42, 0.17);


    int n_phi=1000;
    int n_theta=10000;
    double dtheta=2*M_PI/n_theta;
    double dphi=2*M_PI/n_phi;

    //PovRay pic("slow_piezo_insect.pov");

    for (double theta = 0; theta < 2*M_PI; theta+=dtheta)

    {
        Vector3 n1(1, 0, 0);
        Vector3 n2(0, sin(theta), cos(theta));
        cout << theta << endl;
        //dest << "theta=" << theta << endl << "\t";

        for (double phi = 0; phi < 2*M_PI; phi+=dphi) {

            double a = phi;
            double b = phi+dphi;
            double fin_eps = 1E-7;

            if (deviant(a, n1, n2, t, epsilon, pc)*deviant(b, n1, n2, t, epsilon, pc)<0) {

                double c=0;

                while (abs(b-a)>fin_eps) {
                    c = (a+b)/2;
                    if (deviant(a, n1, n2, t, epsilon, pc)*deviant(c, n1, n2, t, epsilon, pc)<=0)
                        b=c;
                    else
                        a=c;
                }

                //Vector3 nn=n; //*(1000/sqrt(gamms[2]/5.96e3));
                dest /*<< theta*180/M_PI << " "*/ << c*180/M_PI << endl;
                //pic.point(nn);
            }
        }
    }
}

void test_all_roots(){
stringstream sour( "(-1,0) (1.58177e+11,0) (-6.51017e+21,0) (8.00479e+31,0)");
Poly tst(sour);
cout << "final roots = " << tst.all_roots();
}

void test_origin(){
	Tensor t=make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);
	ofstream dest("0.txt");
	for (double phi = 0; phi < 2*M_PI; phi+=0.01) {
		Matrix3 cr=t.crist(Vector3(cos(phi),sin(phi),0));
		Poly pol= cr.characteristic();
		VCD roots = pol.all_roots();
		VD vels = rho2v(pol.all_roots(),5.96E3);
		sort(vels.begin(),vels.end());
		dest << phi << " " << 1/vels[0] << endl;
	}

}

void test_origin_te(){
	Tensor t=make_trigonal_tensor(3.27E10,0.86E10,2.49E10,7.22E10,1.24E10,3.12E10,1.21E10);
    ofstream dest("te_0.txt");

    for (double phi = 0; phi < 2*M_PI; phi+=0.001) {

        double theta = 90.*M_PI/180;

        Vector3 n(cos(phi), sin(theta)*sin(phi), sin(phi)*cos(theta));

        Matrix3 cr=t.crist(n);
		Poly pol= cr.characteristic();
		VCD roots = pol.all_roots();
		VD vels = rho2v(pol.all_roots(),6.25E3);
		sort(vels.begin(),vels.end());
        dest << /*phi << " " <<*/ 1/vels[0] << " " << 1/vels[1] << " " << 1/vels[2] << endl;
	}

}

void test_polarization(){

	Tensor t=make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);
    //double phi=M_PI/8+1e-3;
    //double theta=0.1;

    ofstream dest("pol_teo2.txt");

    for (double phi = 0; phi < M_PI/2; phi+=0.0001) {
    //double phi=M_PI/8+1e-3;
    double theta=90*M_PI/180;
    //Vector3 n(cos(phi),sin(phi),0);

    Vector3 n(cos(phi), sin(theta)*sin(phi), sin(phi)*cos(theta));
    //Vector3 n(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
    //cout << "n= " << n << endl;
	Matrix3 cr=t.crist(n);
    //cout << "matrix= " << endl << cr << endl;

	Poly pol = cr.characteristic();

	VCD roots = pol.all_roots();

    //dest << phi*180/M_PI << " ";

	for (int i = 0; i < 3; i++) {
		double gamma = roots[i].real();

		Vector3 q = get_polarization( cr, gamma);
        //cout << "V= " << V;
        //cout << " q= " << q << endl;

        double vm=abs(q*n);
        double angle=acos(vm);
        dest << angle*180/M_PI << " ";
        }
        dest << endl;
    }
}

void test_polarization_te(){

    Tensor t=make_trigonal_tensor(3.27E10, 0.86E10, 2.49E10, 7.22E10, 1.24E10, 3.12E10, 1.21E10);

	ofstream dest("pol_te.txt");

    for (double phi = 0; phi < M_PI/3; phi+=0.0001) {
	//double phi=M_PI/8+1e-3;
	//double theta=M_PI/2;
    Vector3 n(cos(phi),sin(phi),0);

	Matrix3 cr=t.crist(n);


	Poly pol = cr.characteristic();

	VCD comp_roots = pol.all_roots();
	vector<double> roots(comp_roots.size());
	for(int t = 0; t < int(comp_roots.size()); ++t)
			roots[t] = real(comp_roots[t]);
	sort(roots.begin(),roots.end());

	dest << phi*180/M_PI << " ";
	for (int i = 0; i < 3; i++) {
		double gamma = roots[i];//.real();

		Vector3 q = get_polarization( cr, gamma);

		double vm=abs(q*n);
		double angle=acos(vm);
		//dest << "V= " << V;
		dest << angle*180/M_PI << " ";
		}
		dest << endl;
	}
}


void test_ort(const Vector3& t) {

    Vector3 to=ort(t);

    cout << "t=" << t << " to=" << to << " t*to=" << t*to << endl;
}

void test_orts(){

    test_ort(Vector3(0.718,0.355,0.815));
    test_ort(Vector3(0,-55,0));
    test_ort(Vector3(10,10,10));
}

void test_polarization_te_piezo(){

	Tensor t=make_trigonal_tensor(3.257E10, 0.845E10, 2.57E10, 7.17E10, 1.238E10, 3.094E10, 1.206E10);

	ofstream dest("pol_te_piezo.txt");

    for (double phi = 0; phi < M_PI/3; phi+=0.0001) {
	//double phi=M_PI/8+1e-3;
	//double theta=M_PI/2;
        Vector3 n(cos(phi),sin(phi),0);

	Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
	Tensor3 pc = make_trigonal_piezo_tensor_32(0.42, 0.17);

	Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);

	//Matrix3 cr=t.crist(n);


	Poly pol = tt.characteristic();

	VCD comp_roots = pol.all_roots();
	vector<double> roots(comp_roots.size());
	for(int t = 0; t < int(comp_roots.size()); ++t)
			roots[t] = real(comp_roots[t]);
	sort(roots.begin(),roots.end());

    //dest << phi*180/M_PI << " ";
	for (int i = 0; i < 3; i++) {
		double gamma = roots[i];//.real();

		Vector3 q = get_polarization( tt, gamma);

		double vm=abs(q*n);
		double angle=acos(vm);
		//dest << "V= " << V;
		dest << angle*180/M_PI << " ";
		}
		dest << endl;
	}
}

void test_cinnabar(){
Tensor t=make_tetragonal_tensor(11.7E10,1.7E10,2.9E10,28.9E10,5E10,11.3E10);
	ofstream dest("HgS_2.txt");
	double phi=0;
	for (double theta = 0; theta < 2*M_PI; theta+=0.01) {
		Matrix3 cr=t.crist(Vector3(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)));
		Poly pol= cr.characteristic();
		VCD roots = pol.all_roots();
		VD vels = rho2v(pol.all_roots(),8.1E3);
		sort(vels.begin(),vels.end());
		dest << theta*180/3.1415 << " " << 1/vels[2] << endl;
	}
}

void test_rutil(){
    Tensor t=make_tetragonal_tensor(27.3E10, 17.6E10, 14.9E10, 48.4E10, 19.4E10, 12.5E10);
    ofstream dest("TiO2_z.txt");

    Vector3 n(0,0,1);

    Matrix3 tt =t.crist(n);

    Poly pol= tt.characteristic();
    VCD roots = pol.all_roots();
    VD vels = rho2v(pol.all_roots(),4.25E3);
    sort(vels.begin(),vels.end());
    dest << /*theta*180/3.1415 <<*/ vels[0] << " " << vels[1] << " " << vels[2] << endl;
}

void test_indum(){
    Tensor t=make_tetragonal_tensor(4.53E10, 4E10, 4.15E10, 4.51E10, 1.21E10, 0.651E10);
    ofstream dest("In_z.txt");

    Vector3 n(0,0,1);

    Matrix3 tt = t.crist(n);

    Poly pol= tt.characteristic();
    VCD roots = pol.all_roots();
    VD vels = rho2v(pol.all_roots(),7.28E3);
    sort(vels.begin(),vels.end());
    dest << /*theta*180/3.1415 <<*/ vels[0] << " " << vels[1] << " " << vels[2] << endl;
}

void test_linbo3(){
	Tensor t = make_trigonal_tensor(20.3E10, 5.3E10, 7.5E10, 24.5E10, 0.9E10, 6E10, (20.3E10-5.3E10)/2);

    ofstream dest("linbo3_xy.txt");

    //double theta=0;

    for (double phi = 0; phi < 2*M_PI; phi+=0.01) {

	//Vector3 n(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

    Vector3 n(cos(phi),sin(phi),0);

    Matrix3 epsilon = eps2mat(8.854E-12*44, 8.854E-12*29);
    Tensor3 pc = make_trigonal_piezo_tensor_3m(3.7, 2.5, 0.2, 1.3);

    Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);

        Poly pol = tt.characteristic();

		VCD roots = pol.all_roots();

		VD vels = rho2v(pol.all_roots(),4.76E3);
		sort(vels.begin(),vels.end());
        dest /*<< phi << " "*/ << vels[0] << " " << vels[1] << " " << vels[2] << endl;
    }
}

void test_linbo3_koshibiki(){
	Tensor t = make_trigonal_tensor(19.886E10, 5.467E10, 6.799E10, 23.418E10, 0.783E10, 5.985E10, (19.886E10-5.467E10)/2);

	ofstream dest("linbo3_test.txt");

	//for (double phi = 0; phi < 2*M_PI; phi+=0.01) {

	double phi=30*M_PI/180;
	Vector3 n(cos(phi),sin(phi),0);

	Matrix3 epsilon = eps2mat(8.8542E-12*44.9, 8.8542E-12*26.7);
	Tensor3 pc = make_trigonal_piezo_tensor_3m(3.65, 2.407, 0.328, 1.894);

	Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);


		Poly pol = tt.characteristic();

		VCD roots = pol.all_roots();

		VD vels = rho2v(pol.all_roots(),4.6428E3);
		sort(vels.begin(),vels.end());

		dest << vels[0] << " " << vels[1] << " " << vels[2] << endl;
		dest << n << " " << endl;
		//dest << pol << " " << endl;
		//dest << roots << " " << endl;
	//}
}

void test_te_piezo(){
	Tensor t = make_trigonal_tensor(3.257E10, 0.845E10, 2.57E10, 7.17E10, 1.238E10, 3.094E10, 1.206E10);

    ofstream dest("te_piezo_0.txt");

    for (double phi = 0; phi < 2*M_PI; phi+=0.001) {

        double theta = 90.*M_PI/180;

        Vector3 n(cos(phi), sin(theta)*sin(phi), sin(phi)*cos(theta));

    //Vector3 n(cos(phi),sin(phi),0);

        Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
        Tensor3 pc = make_trigonal_piezo_tensor_32(0.42, 0.17);

        Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);

        Poly pol = tt.characteristic();

		VCD roots = pol.all_roots();

		VD vels = rho2v(pol.all_roots(),6.25E3);
		sort(vels.begin(),vels.end());
        dest <</* phi*180/3.14 << " " <<*/ 1/vels[0] << " " << 1/vels[1] << " " << 1/vels[2] << endl;
	}
}

void test_graph(){
	Tensor t=make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);

	ofstream dest("xz.txt");

	for (double phi=0; phi<2*M_PI; phi+=0.01){
		Vector3 n(cos(phi), 0, sin(phi));

		Matrix3 cr=t.crist(n);

		Poly pol = cr.characteristic();

		VCD roots = pol.all_roots();

		VD gamms = arr_re(roots);
		sort(gamms.begin(),gamms.end());
		Vector3 q = get_polarization( cr, gamms[2]);

		dest << phi*180/M_PI << " " << acos(abs(n*q))*180/M_PI << endl;
	}
}

void test_vectors() {
	Vector3 a(1,2,3);
	Vector3 b(3,2,1);
	Vector3 c=a * 2 + b * 3;
}

void test_slowness() {
	PovRay pov("slow.pov");
	int n=100;
    double rho = 5.96E3;

	Tensor tt = make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);
	//Tensor t = tt * euler( M_PI/2, M_PI/3, 0);

	/*for (int p = 0; p <=n; p++) {
		cout << p << endl;
		for (int q = 0; q < 2*n; q++) {
			Vector3 rn = pq2vec(p,q,n);
			pov.point(rn);
		}
	} */

	pov.begin_mesh();

	for (int q = 0; q < 2*n; q++) {
		Vector3 r1 = pq2vec(0,q,n);
		Vector3 r2 = pq2vec(1,q,n);
		Vector3 r3 = pq2vec(1,q+1,n);
        pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
	}

	for (int p = 1; p < n; p++) {
		for (int q = 0; q < 2*n; q++) {
			Vector3 r1 = pq2vec(p,q,n);
			Vector3 r2 = pq2vec(p+1,q,n);
			Vector3 r3 = pq2vec(p+1,q+1,n);
            pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
		}
	}

	for (int p = 1; p < n; p++) {
		for (int q = 0; q < 2*n; q++) {
			Vector3 r1 = pq2vec(p,q,n);
			Vector3 r2 = pq2vec(p,q+1,n);
			Vector3 r3 = pq2vec(p+1,q+1,n);
            pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
		}
	}

	for (int q = 0; q < 2*n; q++) {
		Vector3 r1 = pq2vec(n-1,q,n);
		Vector3 r2 = pq2vec(n,q,n);
		Vector3 r3 = pq2vec(n-1,q+1,n);
        pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
	}

	pov.end_mesh();
}

void test_slowness_te() {
	PovRay pov("slow_te_0.pov");
    int n=100;
    double rho = 6.25E3;


	Tensor tt = make_trigonal_tensor(3.27E10,0.86E10,2.49E10,7.22E10,1.24E10,3.12E10,1.21E10);

	//Tensor t = tt * euler( M_PI/2, M_PI/3, 0);

	/*for (int p = 0; p <=n; p++) {
		cout << p << endl;
		for (int q = 0; q < 2*n; q++) {
			Vector3 rn = pq2vec(p,q,n);
			pov.point(rn);
		}
	} */

	pov.begin_mesh();

	for (int q = 0; q < 2*n; q++) {
		Vector3 r1 = pq2vec(0,q,n);
		Vector3 r2 = pq2vec(1,q,n);
		Vector3 r3 = pq2vec(1,q+1,n);
        pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
	}

	for (int p = 1; p < n; p++) {
		for (int q = 0; q < 2*n; q++) {
			Vector3 r1 = pq2vec(p,q,n);
			Vector3 r2 = pq2vec(p+1,q,n);
			Vector3 r3 = pq2vec(p+1,q+1,n);
            pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
		}
	}

	for (int p = 1; p < n; p++) {
		for (int q = 0; q < 2*n; q++) {
			Vector3 r1 = pq2vec(p,q,n);
			Vector3 r2 = pq2vec(p,q+1,n);
			Vector3 r3 = pq2vec(p+1,q+1,n);
            pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
		}
	}

	for (int q = 0; q < 2*n; q++) {
		Vector3 r1 = pq2vec(n-1,q,n);
		Vector3 r2 = pq2vec(n,q,n);
		Vector3 r3 = pq2vec(n-1,q+1,n);
        pov.triangle(get_slow(r1,tt,rho),get_slow(r2,tt,rho),get_slow(r3,tt,rho));
    }

	pov.end_mesh();
}

void test_slowness_normal(){

    Vector3 n=Vector3(1,1,1).normalized();

    Tensor t = make_trigonal_tensor(3.257E10, 0.845E10, 2.57E10, 7.17E10, 1.238E10, 3.094E10, 1.206E10);

    Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
    Tensor3 pc = make_trigonal_piezo_tensor_32(0.42, 0.17);

    Vector3 normal=slow_normal(n,0,epsilon,t,pc);

    cout << "n=" << n << endl << "normal=" << normal << endl;
}

void test_piezo_group_vels()
{
    PovRay pov("slow_te_piezo_gr_0.pov");
    int n=100;
    double rho=6.25E3;

    Tensor t = make_trigonal_tensor(3.257E10, 0.845E10, 2.57E10, 7.17E10, 1.238E10, 3.094E10, 1.206E10);

    //ofstream dest("te_piezo_gr.txt");

    //for (double theta = 0; theta < 2*M_PI; theta+=0.01) {

        //for (double phi = 0; phi < 2*M_PI; phi+=0.01) {

            //Vector3 n(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

                //Vector3 n(cos(phi),sin(phi),0);

            Matrix3 epsilon = eps2mat(8.854E-12*23, 8.854E-12*36);
            Tensor3 pc = make_trigonal_piezo_tensor_32(0.42, 0.17);

            //Matrix3 tt = t.make_piezo_crist(n, epsilon, pc);

            //Vector3 normal = slow_normal(n,0,epsilon,t,pc);

            //Poly pol = tt.characteristic();

            //VCD roots = pol.all_roots();

            //VD vels = rho2v(pol.all_roots(),6.25E3);
            //sort(vels.begin(),vels.end());

            //dest /*<< phi*180/3.14 << " "*/ << 1/vels[0] << normal << endl;


            pov.begin_mesh();

            for (int q = 0; q < 2*n; q++) {
                Vector3 r1 = pq2vec(0,q,n);
                Vector3 r2 = pq2vec(1,q,n);
                Vector3 r3 = pq2vec(1,q+1,n);
                pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
            }

            for (int p = 1; p < n; p++) {
                for (int q = 0; q < 2*n; q++) {
                    Vector3 r1 = pq2vec(p,q,n);
                    Vector3 r2 = pq2vec(p+1,q,n);
                    Vector3 r3 = pq2vec(p+1,q+1,n);
                    pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
                }
            }

            for (int p = 1; p < n; p++) {
                for (int q = 0; q < 2*n; q++) {
                    Vector3 r1 = pq2vec(p,q,n);
                    Vector3 r2 = pq2vec(p,q+1,n);
                    Vector3 r3 = pq2vec(p+1,q+1,n);
                    pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
                }
            }

            for (int q = 0; q < 2*n; q++) {
                Vector3 r1 = pq2vec(n-1,q,n);
                Vector3 r2 = pq2vec(n,q,n);
                Vector3 r3 = pq2vec(n-1,q+1,n);
                pov.triangle(get_slow(r1,t,rho),get_slow(r2,t,rho),get_slow(r3,t,rho));
            }

            pov.end_mesh();

        //}
    //}
}

void make_cadr(int j, const Tensor& t, const Vector3& n) {

    stringstream st;
    st << "mult_st\\sphere_travel_";
    st << setfill ('0') << setw (4);
    st << j << ".pov";
    PovRay pov(st.str());

    Matrix3 cr=t.crist(n);
    Poly pol = cr.characteristic();

    VCD roots = pol.all_roots();

    for (int i = 0; i < 3; i++) {
        double gamma = roots[i].real();

        Vector3 q = get_polarization(cr, gamma);

        int r, g, b;
        if (i == 0) {r=1, g=0, b=0;}
        else
            if (i == 1) {r=0, g=1, b=0;}
            else {r=0, g=0, b=1;}
        Vector3 rgb (r,g,b);
        pov.cylinder(n, q*0.25+n, 0.01, rgb);
    }

    pov.cylinder(Vector3(0,0,0), n*1.25, 0.01, Vector3(0.9,0.2,0));
    //pov.cylinder(Vector3(0,0,0), m);
}

void sphere_travel() {

    int i=0;

    Tensor t=make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);

    for (double phi=0; phi<M_PI/2; phi+=0.005){
        ++i;
        Vector3 n(sin(phi), 0, cos(phi));
        make_cadr(i,t,n);
    }

    for (double phi=0; phi<M_PI/4; phi+=0.005){

        ++i;
        Vector3 n(cos(phi), sin(phi),0);
        make_cadr(i,t,n);
    }

    for (double phi=M_PI/2; phi>0; phi-=0.005){

        ++i;
        Vector3 n(cos(M_PI/4)*sin(phi), sin(M_PI/4)*sin(phi), cos(phi));
        make_cadr(i,t,n);
    }
}

void turn_matrix() {

    Tensor t=make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);

    double phi = 10*M_PI/180;
    Matrix3 m = turn_mat(phi);
    Tensor tt=t*m;

    ofstream dest("test_turn.txt");
    for (double phi = 0; phi < 2*M_PI; phi+=0.01) {
        Matrix3 cr=tt.crist(Vector3(cos(phi),sin(phi),0));
        Poly pol= cr.characteristic();
        VCD roots = pol.all_roots();
        VD vels = rho2v(pol.all_roots(),5.96E3);
        sort(vels.begin(),vels.end());
        dest << phi*180/M_PI << " " << 1/vels[0] << endl;
    }
}

void test_make_strain() {

    Tensor t=make_tetragonal_tensor(5.6E10,5.145E10,2.2E10,10.6E10,6.6E10,2.65E10);

    Vector3 n(1,0,0);
    int ind = 0;
    double rho = 5.96E3;

    cout << make_strain(polaris(n,ind,t,rho),slowness_vec(n,ind,t,rho)) << endl;
}

int main(int argc, char* argv[])
{
    test_make_strain();

    //turn_matrix();

    //sphere_travel();

	//test_jum();
	//test_mat();

	//test_vectors();

	//test_Poly();

    //test_polarization();

    //test_polarization_te();

    //test_polarization_te_piezo();

	//test_graph();

    //test_beans();

	//test_slowness();

	//test_slowness_te();

    //test_insect();

    //test_insect_piezo();

    //test_insect_piezo_dihotom();

	//test_cinnabar();

    //test_linbo3();

    //test_rutil();

    //test_indum();

    //test_te_piezo();

    //test_PR_te();

    //group();

	//test_all_roots();

	//test_origin();

    //test_origin_te();

    //test_orts();

    //test_slowness_normal();

    //test_piezo_group_vels();

	//work();

	//getch();

	return 0;
}
//---------------------------------------------------------------------------
