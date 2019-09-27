//////////////////////////////////////////////////////////////////////////////////
//      FUNKIFIT_ROOT v1.0
//      September 2019
//      ROOT version written by J. Bishop
//      email: jackedwardbishop@gmail.com
//      --
//      This header file contains the class definitions needed to call FUNKIFIT
//      The KinFit takes in two arrays, one for the parameters P and one for
//	for the covariance matrix V. Additionally, the mass of the beam and the
//	light products is needed
//////////////////////////////////////////////////////////////////////////////////



#include "TMatrixD.h"
#include <iostream>
#include <TRandom3.h>
using namespace std;
int arrayrow=4;
TRandom3 *rndm;

class KinFit {
public:
	KinFit(double *p,double *V,double Mb,double M1,double M2);//pass parameter and covariance matrix
	~KinFit();
	void kinfit();
	TMatrixD Evaluate(TMatrixD);
	void CalcConstraint(TMatrixD&,TMatrixD);
	void SetConstraint(double *powersin);
	double weighting=.01;
	TMatrixD alpha,alpha0;
	TMatrixD Evaluatedelta(TMatrixD,TMatrixD);
	double firstchisqrd;
	double GetFirstChisqrd() {return firstchisqrd;};
	int parameters;
	int constraints;
	double *powers;
	double Mb,M1,M2;
};

KinFit::~KinFit() {
	return;
}
void KinFit::CalcConstraint(TMatrixD &D,TMatrixD alpha0) {//pb,thetab,phib,p1,theta1,phi1,p2,theta2,phi2
//	H1: p_b - p_1 + p_2 (x) = 0
//	H2: p_b - p_1 + p_2 (y) = 0
//	H3: p_b - p_1 + p_2 (z) = 0
//	H4: p_b^2/(2.*Mb) - p_1^2/(2.*M1) + p_2^2/(2.*M2) = 0
// 	H = 	(sin(thetab)*cos(phib),pb*cos(thetab)*cos(phib),-pb*sin(thetab)*sin(phib),-sin(theta1)*cos(phi1),-p1*cos(theta1)*cos(phi1),p1*sin(theta1)*sin(phi1),-sin(theta2)*cos(phi2),-p2*cos(theta2)*cos(phi2),p2*sin(theta2)*sin(phi2))//px cons
//		(0,1,0,0,-1,0,0,-1,0,0,0,0,0)//py cons
//		(0,0,1,0,0,-1,0,0,-1,0,0,0,0)//pz cons
//		(2*pxb/Mb,2*pyb/Mb,2*pzb/Mb,-2*px1/M1,-2*py1/M1,-2*pz1/M1,-2*px2/M2,-2*py2/M2,-2*pz2/M2,0,0,0,0)//E total
	double pb,p1,p2,thetab,theta1,theta2,phib,phi1,phi2;
	pb=alpha0[0][0];
	thetab=alpha0[1][0];
	phib=alpha0[2][0];
	p1=alpha0[3][0];
	theta1=alpha0[4][0];
	phi1=alpha0[5][0];
	p2=alpha0[6][0];
	theta2=alpha0[7][0];
	phi2=alpha0[8][0];
// H1/dalpha
	D[0][0]=sin(thetab)*cos(phib);//d/dpb
	D[0][1]=pb*cos(thetab)*cos(phib);//d/dthetab
	D[0][2]=-pb*sin(thetab)*sin(phib);//d/dphib
	D[0][3]=-sin(theta1)*cos(phi1);//
	D[0][4]=-p1*cos(theta1)*cos(phi1);//
	D[0][5]=p1*sin(theta1)*sin(phi1);//
	D[0][6]=-sin(theta2)*cos(phi2);
	D[0][7]=-p2*cos(theta2)*cos(phi2);
	D[0][8]=p2*sin(theta2)*sin(phi2);
// H2/dalpha
	D[1][0]=cos(thetab);//d/dpb
	D[1][1]=-pb*sin(thetab);//d/dthetab
	D[1][2]=0;//d/dphib
	D[1][3]=-cos(theta1);//
	D[1][4]=p1*sin(theta1);//
	D[1][5]=0;//
	D[1][6]=-cos(theta2);
	D[1][7]=p2*sin(theta2);
	D[1][8]=0;
// H3/dalpha
	D[2][0]=sin(thetab)*sin(phib);//d/dpb
	D[2][1]=pb*cos(thetab)*sin(phib);//d/dthetab
	D[2][2]=pb*sin(thetab)*cos(phib);//d/dphib
	D[2][3]=-sin(theta1)*sin(phi1);//
	D[2][4]=-p1*cos(theta1)*sin(phi1);//
	D[2][5]=-p1*sin(theta1)*cos(phi1);//
	D[2][6]=-sin(theta2)*sin(phi2);
	D[2][7]=-p2*cos(theta2)*sin(phi2);
	D[2][8]=-p2*sin(theta2)*cos(phi2);
//  H4/dalpha
	D[3][0]=pb/Mb;
	D[3][1]=0.;
	D[3][2]=0.;
	D[3][3]=-p1/M1;
	D[3][4]=0.;
	D[3][5]=0.;
	D[3][6]=-p2/M2;
	D[3][7]=0.;
	D[3][8]=0.;
//	D.Print();
	return;
}


KinFit::KinFit(double p[9],double V[81], double inMb, double inM1, double inM2) {
	Mb=inMb;
	M1=inM1;
	M2=inM2;
	rndm = new TRandom3();
	parameters=9;
	constraints=4;
	const int cparameters=parameters;
	const int cconstraints=constraints;
	const int nparameters=9;
	TMatrixD alpha(nparameters,1);//current parameters
	TMatrixD alpha0(nparameters,1);
	TMatrixD random(nparameters,1);

	const int maxloops=100;
	for(int i=0;i<nparameters;i++) {
		alpha[i][0]=p[i];//fill alpha matrix from the p array input
	}
	TMatrixD v(nparameters,nparameters);//covariance matrix
	TMatrixD valpha(nparameters,nparameters);
	alpha0=alpha;
//
	for(int i=0;i<nparameters*nparameters;i++) {
		int x=floor(i%nparameters);
		int y=i%nparameters;
		valpha[x][y]=V[i];
	}
	valpha.SetMatrixArray(V,"F");
	arrayrow=4;
	TMatrixD D(arrayrow,nparameters);
	CalcConstraint(D,alpha0);//do first calculation of constraints
	bool notconverged=true;
	int loops=0;
	double chisqrd=0,chisqrdlast=0.;
	while(notconverged && loops<maxloops) {
//	Calculate VD
		double Etotal=0.,E1=0.,E2=0.,E3=0.;
		for(int i=0;i<3;i++) {
			Etotal+=(alpha[3*i][0]*alpha[3*i][0])/(2.*Mb)+(alpha[3*i+1][0]*alpha[3*i+1][0])/(2.*M1)+(alpha[3*i+2][0]*alpha[3*i+2][0])/(2.*M2);
			E1+=(alpha[i][0]*alpha[i][0])/(2.*Mb);
			E2+=(alpha[i+3][0]*alpha[i+3][0])/(2.*M1);
			E3+=(alpha[i+6][0]*alpha[i+6][0])/(2.*M2);
		}
		TMatrixD DT(nparameters,arrayrow);
		DT.Transpose(D);//transpose of D (4 x 9)
		TMatrixD Vtemp2 = (D*valpha*DT);//(4 x 9) x (9 x 4) = (4 x 4)
		TMatrixD VD=Vtemp2.Invert();//Eq.13
//	Calculate lambda
		TMatrixD d(arrayrow,1);
		d=Evaluate(alpha0);
		TMatrixD delta_alpha(nparameters,1);
		delta_alpha = alpha-alpha0;//delta alpha_0
		TMatrixD temp1(arrayrow,1);
		temp1 = (d+D*delta_alpha);
		TMatrixD lambda(arrayrow,1);
		lambda=VD*temp1;//Eq.14
		alpha=alpha0-weighting*valpha*DT*lambda;
		delta_alpha = alpha-alpha0;//delta alpha_0
		TMatrixD diffalpha=weighting*valpha*DT*lambda;
                //Update covariance matrix
		valpha=valpha-valpha*DT*VD*D*valpha*weighting;
		chisqrd=0.;
		for(int i=0;i<9;i++) {
			chisqrd+=(valpha[i][i]*valpha[i][i]);
		}
		if(loops==0) firstchisqrd=chisqrd;
		alpha0=alpha;
		loops++;
		if(fabs(chisqrd-chisqrdlast)<1e-9) notconverged=false;//good enough fit
		chisqrdlast=chisqrd;
		CalcConstraint(D,alpha0);//update constraints as these depend on alpha matrix
	}
	for(int i=0;i<nparameters;i++) p[i]=alpha[i][0];//update the p array with the optimized parameters
	delete rndm;
	return;
}


TMatrixD KinFit::Evaluate(TMatrixD alpha0) {
	TMatrixD d(arrayrow,1);
	double pb,p1,p2,thetab,theta1,theta2,phib,phi1,phi2;
	pb=alpha0[0][0];
	thetab=alpha0[1][0];
	phib=alpha0[2][0];
	p1=alpha0[3][0];
	theta1=alpha0[4][0];
	phi1=alpha0[5][0];
	p2=alpha0[6][0];
	theta2=alpha0[7][0];
	phi2=alpha0[8][0];

	d[0][0]=pb*sin(thetab)*cos(phib)-p1*sin(theta1)*cos(phi1)-p2*sin(theta2)*cos(phi2);//px
	d[1][0]=pb*cos(thetab)-p1*cos(theta1)-p2*cos(theta2);//py
	d[2][0]=pb*sin(thetab)*sin(phib)-p1*sin(theta1)*sin(phi1)-p2*sin(theta2)*sin(phi2);//pz
	d[3][0]=(pb*pb)/(2.*Mb)-(p1*p1)/(2.*M1)-(p2*p2)/(2.*M2);//ebeam - erecoillight - erecoilheavy
//	d.Print();
	return d;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
