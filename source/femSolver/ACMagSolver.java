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
		
		boolean dense=true;

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

		double length=model.spaceBoundary[1]-model.spaceBoundary[0];
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
					int neq1=neq+nb;
					int neq2=neq+nb;
					Ks=new SpMatComp(neq1,neq2);
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
			
			Complex jpcw=new Complex(0,pcw);
			Complex smallcomp=new Complex(1e-12,1e-12);

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
				
			}else{

				
				//model.Qs.shownzA();;
	//		model.Ps.shownzA();;
				for(int i=0;i<model.numberOfVarNodes;i++){
					int edgnum=map[i];
					int indx=model.edgeUnknownIndex[map[i]]-1;		
					
					//model.edge[map[i]].node[0].getCoord().hshow();

					SpVectComp spv=new SpVectComp(model.Ps.row[i].times(1),model.RHS.length);
					
			
					for(int k=0;k<spv.nzLength;k++)
					{

						if(spv.index[k]==indx){
							//spv.el[k]=spv.el[k].add(new Complex(-rt.el[i+nuned],0).times(pcw));
							spv.el[k]=spv.el[k].add(new Complex(0,-rt.el[i+nuned]).times(pcw));
						}
					}

					SpVectComp spvq=new SpVectComp(model.Qs.row[i].times(-1),model.RHS.length).times(jpcw);;
					spv=spv.addGeneral(spvq);
		
					if(!dense){


					spv.extend(1);
				    spv.index[spv.nzLength-1]=i+nuned;
					spv.el[spv.nzLength-1]=smallcomp.times(1);


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
			
				//	b.el[i+nuned]=b.el[i+nuned].add(new Complex(pcrhs.el[i+nuned],0));
				//	b.el[i+nuned]=b.el[i+nuned].add(new Complex(0,pcrhs.el[i+nuned]));
				
				b.el[i+nuned]=b.el[i+nuned].add(new Complex(0,pcw*pcrhs.el[i+nuned]));
				}

			 }
		 
	
		if(dense){
			
			
			boolean eigen=false;
			
			if(eigen){
				
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
			
			}else{


				boolean bb=false;
				if(bb){
					int nund = model.numberOfUnknownEdges;
					Mat[] KK = Ks.matForm();

					Mat W = new Mat(Ks.nRow, nund);

					for (int i = 0; i < Ks.nRow; i++)
						for (int j = 0; j < nund; j++) {
							W.el[i][j] = KK[0].el[i][j];
						}
					//W.show();

					util.hshow(W.size());

					Mat Wt = W.transp();
					util.hshow(Wt.size());

					Mat C = Wt.mul(W);
					// C.show();
					util.hshow(C.size());

					Vect br = new Vect(b.length);
					for (int k = 0; k < b.length; k++)
						br.el[k] = b.el[k].re;

					Vect bt = Wt.mul(br);

					MatSolver ms = new MatSolver();

					Vect xr = ms.gaussel(C, bt);
					xr.show();
					
					xc=new VectComp(b.length);
					
			
					
					for(int k=0;k<bt.length;k++)
						xc.el[k]=new Complex(xr.el[k],0);
					
					util.plot(xr);
					
				}
				
				int dof=model.numberOfUnknownEdges;

				
				SpMatComp Ks2=new SpMatComp(Ks.nRow);
				Ks2.symm=false;
		
				int [] nzz2=new int[Ks.nRow];
				for(int i=0;i<Ks.nRow;i++){
					int nz1=Ks.row[i].nzLength;
					nzz2[i]=nz1;
					Ks2.row[i]=new SpVectComp(dof,nz1);
					for(int j=0;j<nz1;j++){
						int cl=Ks.row[i].index[j];
						if(cl<dof){
						
							Ks2.row[i].el[j]=Ks.row[i].el[j].deepCopy();
							Ks2.row[i].index[j]=Ks.row[i].index[j];
						}
					}
				}
				
				
				for(int i=0;i<dof;i++){
					int nz1=Ks.row[i].nzLength;
					for(int j=0;j<nz1-1;j++){
						int cl=Ks.row[i].index[j];
						if(cl<dof){
							
							Ks2.row[cl].extend(1);

							Ks2.row[cl].el[nzz2[cl]]=Ks.row[i].el[j].deepCopy();
							Ks2.row[cl].index[nzz2[cl]]=i;;
							nzz2[cl]++;
							}
						}
					}
						
				//Ks2.shownz();

			//	Ks2.shownz();

			//	Ks2.shownz();
				
				//Ks2.shownz();
				SpMatComp Ks2t=Ks2.transpose(100);
				
			//	Ks2t.shownz();
			

				SpMatComp Cs = new SpMatComp(dof, dof); 
				
				for (int i = 0; i < Ks2t.nRow; i++) {
					if (Ks2t.row[i].nzLength > 0) {
						SpVectComp spv =null;// new SpVect(dof, 100);

						boolean spv1_filled=false;
						SpVectComp spv1 =null;// contact.constraint_matrix_N_trp[contId].row[i].deepCopy();
					

						int kx = 0;

						for (int j = 0; j <= i; j++) {
							if (Ks2t.row[j].nzLength > 0) {
							
								if(!spv1_filled){
									
									 spv = new SpVectComp(dof, 100);
									 
									 spv1 = Ks2.row[i].deepCopy();

									spv1_filled=true;
								}

								Complex dot = spv1.dot(Ks2.row[j]);

								if (dot.norm() == 0)
									continue;

								spv.index[kx] = j;
								spv.el[kx++] = dot;

							}
						
						}
						
						if(spv!=null){
							spv.trim(kx);
							Cs.row[i] = spv.deepCopy();
						}
					}
				}
				

				b =Ks2t.amul(b);

				Ks=Cs.deepCopy();

				Ks.setSymHerm(1); // symmetric but not Hermitian


				int m=b.length;

				
				model.Ci=Ks.scale(b);



				SpMatComp Ls=Ks.ichol(1.05);
				
				Ls.setSymHerm(1); 




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
			
			}
			
		//	xc.show();
			
		}else
		{
		
		
		Ks.shownz();

		//Ks.diagSym().show();
		Ks.setSymHerm(1); // symmetric but not Hermitian

		int m=b.length;
		
		
		model.Ci=Ks.scale(b);



		SpMatComp Ls=Ks.ichol(1.2);
		
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
	//	av2.show();
			
	if(cc1>0) av1=av1.times(1./cc1);
		
		if(cc2>0) av2=av2.times(1./cc2);
		
	//	av.show();
		double pcw2=pcw*pcw;
		
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

		for(int i=0;i<ff.length;i++){
			util.pr(i+"  "+ff.el[i]+"  "+T.el[i]);

		}

		return xc;



	}


}
