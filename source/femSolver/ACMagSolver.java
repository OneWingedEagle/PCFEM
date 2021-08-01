package femSolver;


import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.log10;

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
		
		boolean dense=false;

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
		int NF=model.fdiv;
		double df=0;
		if(NF==0){
			NF=1;
		}else{
			df=(f2-f1)/NF;	
		}
		
	

		Vect T=new Vect(NF);
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
		
		if(photonic==1){
		int nb1=0;
		int nb2=0;
		for(int i=1;i<=model.numberOfEdges;i++){
			if(model.edge[i].edgeKnown) continue;
			
			if(model.edge[i].incident) bedg[nb1++][0]=i;
			if(model.edge[i].exit)  bedg[nb2++][1]=i;


		}
		}
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
				if(model.edge[i].edgeKnown) continue;
				
				if(model.edge[i].incident) nb1++;
				if(model.edge[i].exit) nb2++;
	

			}
					
					util.pr(nb1);
					util.pr(nb2);

		double length=model.spaceBoundary[3]-model.spaceBoundary[2];
		double a=length;
		double eps1=1;
		double k1=Math.sqrt(eps1)*w/a;
		pcw=k1;
		double pcw2=pcw*pcw;
		model.pcw=pcw;
		if(photonic==2){
		 Ks=new SpMatComp(model.Hs.addSmallerNew(model.Ss.timesNew(-pcw2)),model.Ss.timesNew(0));
		}else {
			
			 Ks1=new SpMatComp(model.Hs.addSmallerNew(model.Ss.timesNew(-pcw2)),model.Ss.timesNew(0));





				 neq=Ks1.nRow;

				 nb=nb1+nb2;
					int neq2=neq+nb;
					Ks=new SpMatComp(neq2,neq2);
					for(int i=0;i<neq;i++){
						Ks.row[i]=new SpVectComp(neq2,Ks1.row[i].nzLength);//Ks1.row[i].augh(new SpVect(2*nb,0));
						for(int j=0;j<Ks1.row[i].nzLength;j++){
							Ks.row[i].index[j]=Ks1.row[i].index[j];
							Ks.row[i].el[j]=Ks1.row[i].el[j];
						}

					}

					
					for(int i=0;i<nb1;i++){
	
						int index=model.edgeUnknownIndex[bedg[i][0]]-1;
						
						Ks.row[i+neq]=new SpVectComp(neq2,2);
						Ks.row[i+neq].index[0]=index;
						Ks.row[i+neq].index[1]=i+neq;
						
						Ks.row[i+neq].el[0]=new Complex(1,0);
						Ks.row[i+neq].el[1]=new Complex(1e-12,0);

					}
					
					for(int i=0;i<nb2;i++){
						
						int index=model.edgeUnknownIndex[bedg[i][1]]-1;
						
						Ks.row[i+neq+nb1]=new SpVectComp(neq2,2);
						Ks.row[i+neq+nb1].index[0]=index;
						Ks.row[i+neq+nb1].index[1]=i+neq+nb1;
						
						Ks.row[i+neq+nb1].el[0]=new Complex(1,0);
						Ks.row[i+neq+nb1].el[1]=new Complex(1e-12,0);

					}

	
}
						
				}
		 
 
		if(model.analysisMode==2){
	
			for(int i=0;i<model.numberOfVarNodes;i++){
			
				SpVectComp spv=new SpVectComp(model.Ps.row[i]);
				spv=spv.augh(new SpVectComp(model.Qs.row[i]).times(new Complex(0,wr)));
				Ks.row[i+model.numberOfUnknownEdges]=spv.deepCopy();
			}



		}else if(model.photonic==2){
			
			Complex imag=new Complex(0,1);
			Complex jpcw=new Complex(0,pcw);
			Complex smallcomp=new Complex(0,1e-12);

			 pcrhs=new Vect(model.RHS.length);
			Vect rt=new Vect(model.RHS.length);
			model.magMat.addRHS_SCAT(model,pcrhs,rt);
		
		//	model.RHS=	model.RHS.add(rhs);
			
		//	model.Qs.shownzA();
			int nuned=model.numberOfUnknownEdges;

			boolean test=false;
			
			int [] map=new int[model.numberOfVarNodes];
			for(int i=1;i<=model.numberOfEdges;i++){
			 int nd=model.edge[i].node[0].id;
	
			 int ind=model.nodeVarIndex[nd]-1;
			 if(ind<0) continue;
			 map[ind]=i;
			}
			
			if(test){

			for(int i=0;i<model.numberOfVarNodes;i++){
			
				int index=model.edgeUnknownIndex[map[i]]-1;					
				Ks.row[i+nuned]=new SpVectComp(model.numberOfUnknowns,2);

				Ks.row[i+nuned].index[0]=index;
				Ks.row[i+nuned].index[1]=i+nuned;
				Ks.row[i+nuned].el[0]=new Complex(1,0);
				Ks.row[i+nuned].el[1]=new Complex(1e-12,0);
				}
			}else{

				
				//model.Qs.shownzA();;
	//		model.Ps.shownzA();;
				for(int i=0;i<model.numberOfVarNodes;i++){
					int edgnum=map[i];
					int indx=model.edgeUnknownIndex[map[i]]-1;		
					
					//model.edge[map[i]].node[0].getCoord().hshow();

					SpVectComp spv=new SpVectComp(model.Ps.row[i].times(1),model.RHS.length);
					
/*					spv.show();
					
					int  indx12=model.edgeUnknownIndex[201]-1;		
					int  indx11=model.edgeUnknownIndex[199]-1;	
					int  indx22=model.edgeUnknownIndex[1]-1;		
					int  indx21=model.edgeUnknownIndex[4]-1;		
		*/
					
					for(int k=0;k<spv.nzLength;k++)
					{

						if(spv.index[k]==indx){
						//	spv.el[k]=spv.el[k].add(new Complex(0,-rt.el[i+nuned]).times(1));
							spv.el[k]=spv.el[k].add(new Complex(0,-rt.el[i+nuned]).times(pcw));
						}
					}

					SpVectComp spvq=new SpVectComp(model.Qs.row[i].times(-1),model.RHS.length).times(jpcw);;

		
					if(!dense){
					spv=spv.addGeneral(spvq);

					spv.extend(1);
				    spv.index[spv.nzLength-1]=i+nuned;
					spv.el[spv.nzLength-1]=smallcomp;
					}
			
					Ks.row[i+nuned]=spv.times(1);
		
				//	Ks.row[i+nuned].shownz();

				}

			}
		
				 
		}
		
		
		if(photonic==0 || photonic==2)
			b=new VectComp(model.RHS);
		else
		 b=new VectComp(model.RHS.aug(new Vect(nb)));
		
	
		 if(model.photonic==2){
				int nuned=model.numberOfUnknownEdges;
				for(int i=0;i<model.numberOfVarNodes;i++){
			
				//	b.el[i+nuned]=b.el[i+nuned].add(new Complex(1,0));
				//	b.el[i+nuned]=b.el[i+nuned].add(new Complex(0,pcrhs.el[i+nuned]));
					b.el[i+nuned]=b.el[i+nuned].add(new Complex(pcw*pcrhs.el[i+nuned],0));
				}

			 }
	
		if(dense){
			int nund=model.numberOfUnknownEdges;
			Mat [] KK= Ks.matForm();
			
			Mat W=new Mat(Ks.nRow,nund);
			Mat H=new Mat(Ks.nRow,Ks.nRow);
			
			for(int i=0;i<Ks.nRow;i++)
				for(int j=0;j<nund;j++){
					W.el[i][j]=KK[0].el[i][j];
					H.el[i][j]=KK[0].el[i][j];
				}
		//	TT.show();
			
		
			util.hshow(W.size());

			
			Mat Wt=W.transp();
			util.hshow(Wt.size());

			Mat C =Wt.mul(W).times(1./W.nRow);
			util.hshow(C.size());

			 Eigen eg2=new Eigen(C);
			 
			 Mat Q=eg2.V;
				util.hshow(Q.size());
 
			 Mat Phi=W.mul(Q);
				util.hshow(Phi.size());
 
			 Mat PhiT=Phi.transp();
			 util.hshow(PhiT.size());
			 Mat Mr=PhiT.mul(H.mul(Phi));
			 
	
			// Mr.show();

		//	M.show();
			Vect br=new Vect(b.length);
			for(int k=0;k<b.length;k++)
				br.el[k]=b.el[k].re;
			
			Vect bt =PhiT.mul(br);
			
			MatSolver ms=new MatSolver();

			Vect xr=ms.gaussel(Mr, bt);
		
			xc=new VectComp(b.length);
			
			Vect xr2=Phi.mul(xr);
			

			
			for(int k=0;k<bt.length;k++)
				xc.el[k]=new Complex(xr2.el[k],0);
			
		//	xc.show();
			
		}else
		{
		
	

	//	KK[0].show();
		
		//model.writer.writeMat(KK[0], model.resultFolder+"\\matr.txt");
	//	model.writer.writeMat(KK[1], model.resultFolder+"\\matm.txt");
	//	SpMatComp Ka=new SpMatComp(Ks.nRow,Ks.nRow);
	//	Ka.symm=false;

		Ks.size();
		
	//	Ks.shownz();

		//Ks.diagSym().show();
		Ks.setSymHerm(1); // symmetric but not Hermitian

	//	Ks.shownz();


		
		//model.writer.writeArray(b.el, model.resultFolder+"\\rhs.txt");



		int m=b.length;
		
		double [] rr=new double[m];
		for(int k=0;k<m;k++){
			rr[k]=b.el[k].im;
		}
	//	model.writer.writeArray(rr, model.resultFolder+"\\rhs.txt");
		
		model.Ci=Ks.scale(b);



		SpMatComp Ls=Ks.ichol();
		
		Ls.setSymHerm(1); // symmetric but not Hermitian






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
		
		}

		
	//	xc.show();
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
		//av.show();
			
	if(cc1>0) av1=av1.times(1./cc1);
		
		if(cc2>0) av2=av2.times(1./cc2);
		
	//	av.show();
		
		T.el[kf]=av2.norm2()/(av1.norm2());


	/*
		for(int i=0;i<b.length;i++){
			util.pr(i);

			xc.el[i]=xc1.el[i];
			util.pr(i+" "+model.unknownEdgeNumber[i]);
		}
*/
	//	xc.show();
		}

		util.plot(ff,T);
		T.show();

		return xc;



	}


}