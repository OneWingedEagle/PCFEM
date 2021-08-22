package femSolver;


import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.log10;
import static java.lang.Math.sin;

import fem.Model;
import math.Complex;
import math.Eigen;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatComp;
import math.SpVect;
import math.SpVectComp;
import math.Vect;
import math.VectComp;
import math.util;

public class ACMagSolver {


	public ACMagSolver(){	}

	public VectComp solve(Model model, int step ){
		
		VectComp xc=null;
		VectComp  b=null;
		Vect pcrhs=null;
		model.solver.terminate(false);
		Vect rhsImag=null;
		
		boolean dense=true;
		
		int ntr=2;

		model.setMagMat();

		model.magMat.setRHS(model);

		
		int[] indxT=new int[2*model.RHS.length];
		int[] indxIncid=new int[2*model.RHS.length];
		for(int i=1;i<=model.numberOfEdges;i++){		
			int indx=model.edgeUnknownIndex[i]-1;		
			if(indx>=0 && model.edge[i].exit)
				indxT[indx]=1;
			
			if(indx>=0 && model.edge[i].incident)
				indxIncid[indx]=1;
		}
		
		//int NF=100;
		double pcw=0;
		double f1=model.freq;
		double f2=model.freq2;
		int NF=0;
		double df=0;
		if(model.fdiv==0){
			NF=1;
		}else{
			NF=model.fdiv+1;
			df=(f2-f1)/model.fdiv;	
		}
		
	

		Vect T=new Vect(NF);
		Vect R=new Vect(NF);
		Vect ff=new Vect(NF);
		for(int kf=0;kf<NF;kf++){
			
		model.freq=f1+kf*df;
		 ff.el[kf]=	model.freq;
		util.pr("frequency ="+model.freq);
		double  w=2.*Math.PI*model.freq;
		double  wr=-1./w;
		int photonic=model.photonic;
		int [][] bedg=new int[1+model.numberOfEdges][2];
		
		
		util.pr("photonic:  "+photonic);
		

		SpMatComp Ks=null;
		SpMatComp Ks1;
		int nb=0;
		int neq=0;
		if(photonic==0){
			Ks=new SpMatComp(model.Hs,model.Ss.timesNew(w));
		}
		else{
			int nb1=0;
			int nb2=0;
			for(int i=1;i<=model.numberOfEdges;i++){
				//if(model.edge[i].edgeKnown) continue;
				
				if(model.edge[i].incident) nb1++;
				if(model.edge[i].exit) nb2++;
	

			}
					
					util.pr(nb1);
					util.pr(nb2);

		double width=model.spaceBoundary[1]-model.spaceBoundary[0];
		double a=width;
		double eps1=1;
		double k1=Math.sqrt(eps1)*w/a;
		pcw=k1;
		double pcw2=pcw*pcw;
		model.pcw=pcw;
		
			
			util.pr(pcw);
			
		model.magMat.setRHS_SCAT(model,1);
		//model.RHS.zero();
		rhsImag=model.RHS.deepCopy();
	//	rhsImag.show();
		model.magMat.setRHS_SCAT(model,0);
		//model.RHS.zero();
		
		 Ks=new SpMatComp(model.Hs.addSmallerNew(model.Ss.timesNew(-pcw2)),model.Ss.timesNew(0));
								
	}
		 
 
		if(model.analysisMode==2){
	
			for(int i=0;i<model.numberOfVarNodes;i++){
			
				SpVectComp spv=new SpVectComp(model.Ps.row[i]);
				spv=spv.augh(new SpVectComp(model.Qs.row[i]).times(new Complex(0,wr)));
				Ks.row[i+model.numberOfUnknownEdges]=spv.deepCopy();
			}



		}else if(model.photonic==2){
			
			SpMatComp  Ks2=new SpMatComp(Ks.nRow+ntr,Ks.nRow+ntr);
			
			for(int i=0;i<Ks.nRow;i++){
				SpVectComp spv=new SpVectComp(Ks2.nRow,Ks.row[i].nzLength);
				for(int k=0;k<spv.nzLength;k++){
					spv.el[k]=Ks.row[i].el[k];
					spv.index[k]=Ks.row[i].index[k];
				}
				Ks2.row[i]=spv.deepCopy();
			}
			

	
			
			Complex jpcw=new Complex(0,pcw);
			Complex smallcomp=new Complex(1e-12,1e-12);
			Complex img=new Complex(0,pcw);

			 pcrhs=new Vect(model.RHS.length);
			Vect rt=new Vect(model.RHS.length);
			model.magMat.addRHS_SCAT(model,pcrhs,rt);
		
		//	model.RHS=	model.RHS.add(rhs);
			
		//	model.Qs.shownzA();
			int nuned=model.numberOfUnknownEdges;

			
			int [] map=new int[ntr];
			for(int i=1;i<=model.numberOfEdges;i++){
			 int nd=model.edge[i].node[0].id;
	
			 int ind=-1;
			 if(model.edge[i].incident) ind=0;
			 else if(model.edge[i].exit) ind=1;
			 if(ind<0) continue;
			 map[ind]=i;
			}

				for(int i=0;i<ntr;i++){
					int edgnum=map[i];
					int indx=model.edgeUnknownIndex[map[i]]-1;		
					
					//model.edge[map[i]].node[0].getCoord().hshow();

					SpVectComp spv=new SpVectComp(model.Ps.row[i].times(-1),model.RHS.length);
					
			
					for(int k=0;k<spv.nzLength;k++)
					{

						if(spv.index[k]==indx){
							spv.el[k]=spv.el[k].add(new Complex(rt.el[i+nuned],0).times(jpcw));
							break;
						}
					}

					SpVectComp spvq=new SpVectComp(model.Qs.row[i].times(-1),model.RHS.length).times(jpcw);;
					spv=spv.addGeneral(spvq);
					spv.hshow();
					Ks2.row[i+nuned]=spv.times(1);

		
				//	Ks.row[i+nuned].shownz();

				}


				Ks=Ks2.deepCopy();
				
				Ks.show();
		
				 
		}
		
		
		if(photonic==0 || photonic==2)
			b=new VectComp(model.RHS,rhsImag);
		else
		 b=new VectComp(model.RHS.aug(new Vect(nb)));

	
		 if(model.photonic==2){
				int nuned=model.numberOfUnknownEdges;
				for(int i=0;i<model.numberOfVarNodes;i++){
			
			
				//b.el[i+nuned]=b.el[i+nuned].add(new Complex(0,pcw*pcrhs.el[i+nuned]));
				}

			 }
		 

	//	Ks.shownz();

		//Ks.diagSym().show();
		Ks.setSymHerm(1); // symmetric but not Hermitian

		int m=b.length;
		
		
		model.Ci=Ks.scale(b);



		SpMatComp Ls=Ks.ichol(1.02);
		
		Ls.setSymHerm(1); // symmetric but not Hermitian




//b.show();

		if(b.norm()>1e-8){
			xc=model.solver.COICCG(Ks,Ls,b,model.errCGmax,model.iterMax,new VectComp(m),1,false);
			//xc=model.solver.COCG(Ks,b,model.errCGmax,model.iterMax,new VectComp(m),1,false);
		}
		else{
			xc=new VectComp(m);
			model.solver.totalIter++;
			model.solver.errs.add(0.);
			model.solver.totalIter++;
			model.solver.errs.add(log10(model.errCGmax));
			model.solver.errs.add(0.);
			if(model.hasBunif) model.scaleKnownEdgeAL(0);
		}
			

		xc.timesVoid(model.Ci);	
		
		
		
		
	//	xc.show();	
		double y0=model.spaceBoundary[2];
		double y1=model.spaceBoundary[3];
		double L=y1-y0;
		
		for(int i=0;i<model.numberOfUnknownEdges;i++){	
			
			int iedge=model.unknownEdgeNumber[i+1];
			double y=model.edge[iedge].node[0].getCoord(1)-y0;
		
			double E0r=cos(-pcw*y)*(1-1*y/L);
			double E0m=sin(-pcw*y)*(1-1*y/L);
			xc.el[i].re+=E0r;
			xc.el[i].im+=E0m;
			
		}
		

		//xc.show();

		Complex av2=new Complex(0.,0);
		Complex av1=new Complex(0.,0);
		int cc2=0;
		int cc1=0;
		for(int i=0;i<b.length;i++){
			if(indxT[i]==1) {
				av2=av2.add(xc.el[i]);
				cc2++;
			}
			if(indxIncid[i]==1) {
				av1=av1.add(xc.el[i]);
				cc1++;
			}
		
		}

			
		if(cc1>0) av1=av1.times(1./cc1);
		
		if(cc2>0) av2=av2.times(1./cc2);
		
		av2.show();
		
	//	av.show();
		double pcw2=pcw*pcw;
			
		T.el[kf]=av2.norm2();///(av1.norm2());

	//	R.el[kf]=av2.norm2();
	//	xc.show();
		}

		util.plot(ff,T);
		//util.plot(ff,R);

		for(int i=0;i<ff.length;i++){
			util.pr(i+"  "+ff.el[i]+"  "+T.el[i]);

		}

		return xc;



	}


}
