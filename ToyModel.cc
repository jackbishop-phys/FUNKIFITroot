//////////////////////////////////////////////////////////////////////////////////
//	FUNKIFIT_ROOT v1.0
//	September 2019
//	ROOT version written by J. Bishop
//	email: jackedwardbishop@gmail.com
//	--
//	This ROOT macro demonstrates the interfacing to the kinfit.h kinematic
//	fitting program.
//	This toy model generates some events, smears the energies/angles and plots
//	the spread before and after the KF process.
//////////////////////////////////////////////////////////////////////////////////


#include "kinfit.h"
#include "TRandom3.h"
TRandom3 *rnd3;
TH2F *hEtheta;
TH2F *hthetatheta;
TH1F *hECM1,*hECM2,*hthetaerrl,*hthetaerrh,*hthetaerrKFl,*hthetaerrKFh,*hphierr,*hphierrKF;
double d2r=45./atan(1.);
void Event();
void ToyModel() {
	rnd3 = new TRandom3();
	hEtheta = new TH2F("hEtheta","hEtheta",100,0,90,100,0,20);
	hthetatheta = new TH2F("hthetatheta","hthetatheta",100,0,90,100,0,90);
	hECM1 = new TH1F("hECM1","Centre of mass - Before",1000,0,5);
	hECM2 = new TH1F("hECM2","Centre of mass - Before",1000,0,5);
	hthetaerrl = new TH1F("hthetaerrl","Theta err",200,-10,10);
	hthetaerrh = new TH1F("hthetaerrh","Theta err",200,-10,10);
	hthetaerrKFl = new TH1F("hthetaerrKFl","Theta err",200,-10,10);
	hthetaerrKFh = new TH1F("hthetaerrKFh","Theta err",200,-10,10);
	hphierr = new TH1F("hphierr","Phi err",100,-10,10);
	hphierrKF = new TH1F("hphierrKF","Phi err",100,-10,10);
	const int numb_events=1000;
	for(int i=0;i<numb_events;i++) {
		if(i%100==0) cout<<i<<" events done"<<endl;
		Event();
	}
	TCanvas *c1 = new TCanvas("c1","ECM",0,0,800,600);
	hECM2->Draw();
	hECM2->SetLineColor(kRed);
	hECM1->Draw("SAME");
        TLegend *leg= new TLegend(0.1,0.7,0.4,0.9);
        leg->SetTextSize(0.05);
        leg->AddEntry(hECM1,"Before");
        leg->AddEntry(hECM2,"After");
        leg->Draw("SAME");
	TCanvas *c2 = new TCanvas("c2","Theta",0,0,800,600);
	hthetaerrKFh->Draw();
	hthetaerrKFl->SetLineColor(kRed);
	hthetaerrKFl->Draw("SAME");
	hthetaerrKFh->SetLineColor(kMagenta);
	hthetaerrl->SetLineColor(kBlack);
	hthetaerrl->Draw("SAME");
	hthetaerrh->SetLineColor(kBlue);
	hthetaerrh->Draw("SAME");
        TLegend *leg2= new TLegend(0.1,0.7,0.4,0.9);
        leg2->SetTextSize(0.05);
        leg2->AddEntry(hthetaerrl,"Before - light");
        leg2->AddEntry(hthetaerrh,"Before - heavy");
        leg2->AddEntry(hthetaerrKFl,"After - light");
        leg2->AddEntry(hthetaerrKFh,"After - heavy");
        leg2->Draw("SAME");
	return;
}
void Event() {
//      Binary reaction mechanism for beam + (Zt,At) -> (A1,Z1,Ex1) + (A2,Z2,Ex2)
        double Mt=4.,M1=4,M2=6,Mb=6;
        double benergy=10;//(MeV) - beam energy
////Do some kinematics
        double CM_theta=acos(-1.+2.*rnd3->Rndm());//0->pi relative to the particle direction
        double CM_thetah=4.*atan(1.)-CM_theta;//pi-cm_theta_light
        double CM_psi=rnd3->Rndm()*4.*2.*atan(1.);//0->2pi
        double CM_psih=8.*atan(1.)-CM_psi;//2pi-cm_psi_light
///Rotate so theta/psi are relative to the original beam direction
        double dir[3]={sin(CM_theta)*sin(CM_psi),sin(CM_theta)*cos(CM_psi),cos(CM_theta)};
//Set Q-value
        double Q_value=0;
        double E_cm = (Mt*benergy/(M1+M2))+Q_value;//new CM energy MeV
        double p_1 = sqrt(2.*Mt*E_cm*(1.*M2/(M1+M2)));//E_1 = m2/(m1+m2) * E_t
        double p_new_1[4],p_new_2[4];
	for(int i=0;i<3;i++) {
		p_new_1[i] = p_1*dir[i];// new momentum of scattered Be in COM
		p_new_2[i] = -p_1*dir[i];// new momentum of scattered Be in COM
	}
        p_new_1[2]+=sqrt(2.*Mb*benergy)*(1.*M1/(M1+M2));
        p_new_2[2]+=sqrt(2.*Mb*benergy)*(1.*M2/(M1+M2));
//
	p_new_1[3]=sqrt(p_new_1[0]*p_new_1[0]+p_new_1[1]*p_new_1[1]+p_new_1[2]*p_new_1[2]);
	p_new_2[3]=sqrt(p_new_2[0]*p_new_2[0]+p_new_2[1]*p_new_2[1]+p_new_2[2]*p_new_2[2]);
	double theta[2],phi[2];
	theta[0]=acos(p_new_1[2]/p_new_1[3]);
	theta[1]=acos(p_new_2[2]/p_new_2[3]);
	phi[0]=atan2(p_new_1[1],p_new_1[0]);
	phi[1]=atan2(p_new_2[1],p_new_2[0]);
	double energy[2];
	double thetatrue[2],phitrue[2];
	for(int i=0;i<2;i++) {
		thetatrue[i]=theta[i];
		phitrue[i]=phi[i];
	}
//	Do the smearing here! It is important to see that if the sigma here isn't the same as the
//	sigma fed into the covariance matrix V then the convergence will not be as good!
//	Smear angles
	theta[0]+=rnd3->Gaus(0,1/d2r);
	phi[0]+=rnd3->Gaus(0,1/d2r);
	theta[1]+=rnd3->Gaus(0,.5/d2r);
	phi[1]+=rnd3->Gaus(0,.5/d2r);
	hthetaerrl->Fill((thetatrue[0]-theta[0])*d2r);
	hthetaerrh->Fill((thetatrue[1]-theta[1])*d2r);
//	Smear energies
	energy[0]=((p_new_1[3]*p_new_1[3])/(2.*M1))+rnd3->Gaus(0,0.2);
	energy[1]=((p_new_2[3]*p_new_2[3])/(2.*M2))+rnd3->Gaus(0,0.2);
	hEtheta->Fill(theta[0]*d2r,energy[0]);
	hEtheta->Fill(theta[1]*d2r,energy[1]);
	hthetatheta->Fill(theta[0]*d2r,theta[1]*d2r);
	double V[81],p[9];//Defines the covariance matrix and the variables
//
	p[0]=sqrt(2.*Mb*(benergy+rnd3->Gaus(0,0.0*(M1+M2)/Mb)));//pbeam
	p[0]=sqrt(2.*Mb*benergy);//pbeam
	p[1]=0.;//theta
	p[2]=0.;//phi
	p[3]=sqrt(2.*M1*energy[0]);//plight
	p[4]=theta[0];//theta
	p[5]=phi[0];//phi
	p[6]=sqrt(2.*M2*energy[1]);//pheavy
	p[7]=theta[1];//theta
	p[8]=phi[1];//phi
	double px,py,pz;
	px=p[3]*sin(theta[0])*cos(phi[0])+p[6]*sin(theta[1])*cos(phi[1]);
	py=p[3]*sin(theta[0])*sin(phi[0])+p[6]*sin(theta[1])*sin(phi[1]);
	pz=-p[0]+p[3]*cos(theta[0])+p[6]*cos(theta[1]);
	hECM1->Fill(((p[0]*p[0])/(2.*Mb))*Mt/(Mt+Mb));
        for(int i=0;i<9;i++) {
                for(int j=0;j<9;j++) {
                        V[9*i+j]=0.000001;//off diagonal elements have no co-correlation here but a small number here does no effect the algorithm
                        if(i==j && (i==0)) V[9*i+j]=TMath::Power(0.05,2.);//beam momentum - 50 keV
                        if(i==j && (i==3)) V[9*i+j]=TMath::Power(0.2,2.);//light momentum - 200 keV
                        if(i==j && (i==6)) V[9*i+j]=TMath::Power(0.2,2.);//heavy momentum - 200 keV
                        if(i==j && (i==4 || i==5)) V[9*i+j]=TMath::Power(1./d2r,2.);//1.0 degrees
                        if(i==j && (i==7 || i==8)) V[9*i+j]=TMath::Power(.5/d2r,2.);//0.5 degrees
                }
        }
//	Perform the kinematic fitting! 9 parameters, 4 constraint equations, P is the parameter array, V is the covariance array
        KinFit *kinematic = new KinFit(p,V,Mb,M1,M2);
	px=p[3]*sin(theta[0])*cos(phi[0])+p[6]*sin(theta[1])*cos(phi[1]);
	py=p[3]*sin(theta[0])*sin(phi[0])+p[6]*sin(theta[1])*sin(phi[1]);
	pz=-p[0]+p[3]*cos(theta[0])+p[6]*cos(theta[1]);
	theta[0]=p[4];
	phi[0]=p[5];
	theta[1]=p[7];
	phi[1]=p[8];
	hthetaerrKFl->Fill((thetatrue[0]-theta[0])*d2r);
	hthetaerrKFh->Fill((thetatrue[1]-theta[1])*d2r);
	p_new_1[0]=p[3]*sin(theta[0])*cos(phi[0]);
	p_new_1[1]=p[3]*sin(theta[0])*sin(phi[0]);
	p_new_1[2]=p[3]*cos(theta[0]);
	p_new_2[0]=p[6]*sin(theta[1])*cos(phi[1]);
	p_new_2[1]=p[6]*sin(theta[1])*sin(phi[1]);
	p_new_2[2]=p[6]*cos(theta[1]);
	energy[0]=(p[3]*p[3])/(2.*M1);
	energy[1]=(p[6]*p[6])/(2.*M2);
	hECM2->Fill((energy[0]+energy[1])*Mt/(M1+M2));
	return;
}
