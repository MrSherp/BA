#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdlib.h> 
#include <vector>   
#include <map> 
#include <iostream>
#include "FEM.h"

using namespace std;

Vector gradient(Vector v){
	
	return v;
}

Vector div(Vector p){
	
	return p;
}

Vector projK(Vector p, Vector uL, Vector uR, double eta, int N){
	
	Vector tmp(N*N);
	tmp.setZero();
	tmp.copy(uL);
											//fehlt -w(uR)
	
	for(int i=0; i<(N*N); i++){
		
		if(eta < abs(p[i*2])){
			
			p[i*2]=eta;
		}
		else{}
		
		if(p[i*2+1] < -abs(tmp[i])){
			
			p[i*2+1] = -abs(tmp[i]);
		}
	}
	
	return p;
}

Vector projC(Vector v, int N){
	
	for(int i=0; i < v.size(); i++){
		
		v[i*N] = 1.;
		v[i*N+N-1] = 0.;
	}
	
	return v;
}

void algo(Vector Omega, int N, Vector uL, Vector uR, double eta) { 	
	
	double h;									//Gitterweite
	h = (Omega[1]-Omega[0])/N;
	
	Vector v(N*N);
	Vector phi(N*N*2);
	Vector w(N*N);
	Vector diff(N*N);
		
	//Startvektoren fehen noch
	
	double t;
	double s;
		
	if (sqrt(8.)*t*s/h >= 1)
	{
		cout << "Die Zeitschritte sind zu groß" << endl;
		t= h/sqrt(8.)/s/2.;
	}
		
	diff.setAll(1.);
		
	Vector tmp1(2*N*N); Vector tmp2(N*N); Vector tmp3(N*N);	
	tmp1.setZero(); tmp2.setZero(); tmp3.setZero();	
	
	while(diff.norm() >= 0.00001) {
			
		tmp1.copy(phi);
		tmp1.addMultiple(gradient(w),s);
		phi.copy(projK(tmp1,uL,uR,eta,N));		//phi_n+1
		
		tmp3.copy(v);
		diff.copy(v);
		
		tmp2.copy(v);						
		tmp2.addMultiple(div(phi),-t);	
		v.copy(projC(tmp2, N));					//v_n+1

		tmp3.operator*=(-1.);				
		tmp3.addMultiple(v,2.);
		w.copy(tmp3);							//w_n+1
		
		diff.addMultiple(v,-1.);
	}
		
	v.print();
	phi.print();		
		
}




int main()

{
	// Eingabe h, Omega festlegen, Eingabe u, Eingabe eta
	
	//algo(Omega, N, uL, uR, eta);
	
	//test
	Vector v(8);
	v.setAll(5.);
	v.pushBack(2.5);
	v.print();
	
	Vector w(9);
	w.copy(projC(v,3));
	w.print();
	v.print();
	
	//test end
	
	
	printf("Blah!\n");
	cout << "Nochmal blah" << endl;
	
	
	
	
	
	return 0;
	
}		