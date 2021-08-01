package femSolver;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.log10;

import fem.Model;
import fem.TimeFunction;
import io.Writer;
import math.Eigen;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class StaticLinearMagSolver{
	int stepNumb;
	boolean usePrev=true;

	public StaticLinearMagSolver(){	}

	public Vect solve(Model model,int step,Vect x_init){

		
		this.stepNumb=step;

		if(!usePrev || x_init.length==0){
			x_init=new Vect(model.numberOfUnknowns);

		}

		if(model.numberOfUnknowns==0)
			return new Vect(model.numberOfUnknowns);
	
	SpMat L=new SpMat();

	Vect x=new Vect(model.numberOfUnknowns);

	model.solver.terminate(false);

	

if(step==0){
	model.setMagMat();

}
	
model.magMat.setRHS(model);



/*  int [] ndedge=new int[1+model.numberOfNodes];
	for(int i=1;i<=model.numberOfEdges;i++){
		int n0=model.edge[i].node[0].id;
		ndedge[n0]=i;
	//	util.pr(i+" "+n0);
		
	}*/
	
	  int [] bedg=new int[1+model.numberOfEdges];
	int nb=0;
	for(int i=1;i<=model.numberOfEdges;i++){
		/*
		//if(!model.edge[i].node[0].onBound[2] &&!model.edge[i].node[0].onBound[3] &&model.edge[i].node[0].getCoord(0)<1e-2){
			if(!model.edge[i].node[0].onBound[2] &&!model.edge[i].node[0].onBound[3] &&model.edge[i].node[0].getCoord(0)<1e-2){
				bedg[nb]=i;
			bedg[nb]=i;
			nb++;
		}
		if(!model.edge[i].node[0].onBound[2] &&!model.edge[i].node[0].onBound[3] &&model.edge[i].node[0].getCoord(0)>4e-2){
			bedg[nb]=i;
			nb++;
		}	*/	
		for(int j=0;j<4;j++){
			if(/*!model.edge[i].node[0].onBound[2] &&!model.edge[i].node[0].onBound[3] &&*/model.edge[i].node[0].onBound[j]){
				bedg[nb]=i;
				nb++;
			}		
				
		}
	
/*		if(model.edge[i].node[0].getCoord(0)<-.0001){
			bedg[nb]=i;
			nb++;
		}		
			if(model.edge[i].node[0].getCoord(0)>.0493){
				bedg[nb]=i;
				nb++;
			}	*/	
				
		
	}
util.pr(nb);
	//	util.pr("|RHS|="+model.RHS.norm());


	SpMat  Ks1=model.Hs.deepCopy();

	int neq=Ks1.nRow;
/*	int nb=16;
	int neq2=neq+2*nb;
	int[][] nnb=new int[nb][2];
	for(int i=0;i<nb;i++){
		nnb[i][0]=1+28*i;
		nnb[i][1]=28+28*i;
	}*/
	
	int neq2=neq+nb;
	SpMat Ks=new SpMat(neq2,neq2);
	for(int i=0;i<neq;i++){
		Ks.row[i]=new SpVect(neq2,Ks1.row[i].nzLength);//Ks1.row[i].augh(new SpVect(2*nb,0));
		for(int j=0;j<Ks1.row[i].nzLength;j++){
			Ks.row[i].index[j]=Ks1.row[i].index[j];
			Ks.row[i].el[j]=Ks1.row[i].el[j];
		}

	}

	
	for(int i=0;i<nb;i++){
	//	util.pr(i);
		
	//	Ks.row[i].index[Ks1.row[i].nzLength]=neq+i;
	//	Ks.row[i].el[Ks1.row[i].nzLength]=1;
		
		int index=model.edgeUnknownIndex[bedg[i]]-1;
	//int ed2=ndedge[nnb[i][1]];
		
		Ks.row[i+neq]=new SpVect(neq2,2);
//		Ks.row[i+neq+nb]=new SpVect(neq2,2);
		Ks.row[i+neq].index[0]=index;
		//Ks.row[i+neq+nb].index[0]=model.edgeUnknownIndex[ed2];
		Ks.row[i+neq].index[1]=i+neq;
		//Ks.row[i+neq+nb].index[1]=i+neq+nb;
		
		Ks.row[i+neq].el[0]=1;
		//Ks.row[i+nb+neq].el[0]=1;
		Ks.row[i+neq].el[1]=0;
	//	Ks.row[i+nb+neq].el[1]=0;

	}
	
	
//	Ks1.shownzA();
	///Ks.shownzA();

	Ks.size();
//Ks.shownzA();
//	Ks.matForm().show();
//	util.pr(9999);
	/*//model.RHS.show();
	//Ks.shownz();
	//Ks.plot();
	Mat M1=Ks.matForm();
	for(int i=0;i<M1.nRow;i++)
		for(int j=i+1;j<M1.nCol;j++)
			M1.el[i][j]=M1.el[j][i];

	Mat M=new Mat(model.numberOfEdges,model.numberOfEdges);
	M.eye();
	
	for(int i=0;i<model.numberOfEdges;i++){
		for(int j=0;j<=i;j++){
	int rind=model.edgeUnknownIndex[i+1]-1;
	int cind=model.edgeUnknownIndex[j+1]-1;
	if(rind>=0 && cind>=0)
		M.el[i][j]=M1.el[rind][cind];
		}
	}
	
	SpMat Ms=new SpMat(M);
//	Ms.plot();
	
	for(int i=1;i<=model.numberOfEdges;i++){
	
		if(model.edge[i].edgeKnown) util.pr(i);
	}
	
	Vect vs=new Vect(model.numberOfEdges);
	for(int i=0;i<model.numberOfEdges;i++){

	int rind=model.edgeUnknownIndex[i+1]-1;
	if(rind>=0)
	vs.el[i]=model.RHS.el[rind];
	}*/
//Ms.shownz();
	//M.show();
	
//	vs.show();
	Vect rhs=model.RHS.aug(new Vect(nb));
	
	//rhs=Ks.smul(new Vect().ones(neq2));
	
	Vect Ci=Ks.scale(rhs);
	
	x_init=new Vect(neq2);
////	Mat M=Ks.matForm();
	//M.show();
//	for(int i=0;i<M.nRow;i++)
//		for(int j=i+1;j<M.nCol;j++) M.el[i][j]=M.el[j][i];
	//
//	util.hshow(M.size());
	
	//Eigen eg1=new Eigen(M);
	//eg1.lam.show();

//model.writer.writeArray(M.el, "D:\\JavaWorks\\FEM problems\\EM\\edge element\\mat.txt");
	//util.pr("|RHS|="+model.RHS.norm());
	x_init.timesVoid(Ci.inv());
		
	L=Ks.ichol();

	Vect xag=null;


	if(rhs.abs().max()>1e-8){
			if(!usePrev || model.xp==null){
				xag=model.solver.ICCG(Ks,L,rhs,model.errCGmax,model.iterMax,x_init);
				//x=model.solver.CG(Ks, model.RHS,model.errCGmax,model.iterMax,x_init);
				//xag.show();
			}
			else{
				xag=model.solver.ICCG(Ks,L, rhs,model.errCGmax,model.iterMax,model.xp);
				//x=model.solver.err0ICCG(Ks,L, model.RHS,1e-2*model.errCGmax,model.iterMax,model.xp);	
			
			}
	}else{
		
		x=new Vect(model.numberOfUnknowns);
		
		model.solver.totalIter++;
		model.solver.errs.add(0.);
		model.solver.totalIter++;
		model.solver.errs.add(log10(model.errCGmax));
		model.solver.errs.add(0.);
		if(model.hasBunif) model.scaleKnownEdgeAL(0);
	}



		model.xp=xag.deepCopy();


		xag.timesVoid(Ci);
		
		for(int i=0;i<neq;i++)
			x.el[i]=xag.el[i];
		

		//util.pr("|x|="+x.norm());

boolean unif=false;

if(unif){
x.zero();
Vect u=new Vect(0,0,1);

double Bx=1;
double By=0;
double Bz=0;
if(model.dim==3)
	Bz=model.unifB.el[2];
double Ax,Ay,Az;
double x1,y,z;
Vect A=new Vect(3);
int[] edgeDirs=new int[1+model.numberOfEdges];

for(int i=1;i<=model.numberOfElements;i++){
	boolean[] edgeDir=model.element[i].getEdgeReverse();
	int[] ne=model.element[i].getEdgeNumb();
	for(int j=0;j<model.nElEdge;j++)
		if(edgeDir[j])
		edgeDirs[ne[j]]=-1;
		else	edgeDirs[ne[j]]=1;
}

for(int i=1;i<=model.numberOfEdges;i++){


	Vect edgeVect=model.edge[i].node[1].getCoord().sub(model.edge[i].node[0].getCoord());

	Vect center=model.edge[i].node[1].getCoord().add(model.edge[i].node[0].getCoord()).times(.5);
	
	x1=center.el[0];
	y=center.el[1];

	Az=y*Bx-x1*By;
	if(model.dim==3){
		z=center.el[2];
		Ax=x1*By;
		Ay=y*Bz;
		A=new Vect(Ax,Ay,Az);
	}else{
		A=new Vect(0,0,Az);
	}


	double a=edgeVect.dot(A);
	

	if(model.edgeUnknownIndex[i]>0)
		x.el[model.edgeUnknownIndex[i]-1]=a;
		else{
		model.edge[i].setA(a);
		}




}
	

	}
	

/*for(int i=1;i<=model.numberOfEdges;i++){

//util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")= "+model.edge[i].getA());
util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")");

}*/

boolean edgeNodes=false;
if(edgeNodes)
for(int i=1;i<=model.numberOfElements;i++){

	int[] ne=model.element[i].getEdgeNumb();
	//util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")= "+model.edge[i].getA());
	util.pr("element "+i);
	util.pr(" edges:");
	for(int j=0;j<model.nElEdge;j++)
		util.ph(ne[j]+", ");
	util.pr("");
	
	int[] nn=model.element[i].getVertNumb();
	//util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")= "+model.edge[i].getA());
	util.pr(" nodes:");
	for(int j=0;j<model.nElVert;j++)
		util.ph(nn[j]+",");
	util.pr("");
	
	
	//util.pr("element "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")");


	}
	

	model.setSolution(x);	

	
		System.out.println("Bmax ( linear analysis): "+model.Bmax);
		


		return x;



}




}
