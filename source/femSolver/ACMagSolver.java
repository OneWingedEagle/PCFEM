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
		
		util.pr("============== "+model.numberOfVarNodes);
		//model.R0=0;
	//	model.T0=1;//Math.sqrt(1-model.R0*model.R0);
		
		boolean dense=false;

		model.setMagMat();

		model.magMat.setRHS(model);
		
		model.RHS=model.RHS.sub(model.RHS_boundary);
		//model.RHS_boundary.show();
		
		int[] indxT=new int[2*model.RHS.length];
		int[] indxIncid=new int[2*model.RHS.length];
		for(int i=1;i<=model.numberOfEdges;i++){		

			int indx=model.edgeUnknownIndex[i]-1;		
			if(indx>=0 && model.edge[i].exit){
				indxT[indx]=1;
			}
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
	
		
		util.pr("photonic:  "+photonic);

		SpMatComp Ks=null;
		int nb=0;
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
			
		model.magMat.setRHS_IncidentField(model,1);

		rhsImag=model.RHS.deepCopy();

		model.magMat.setRHS_IncidentField(model,0);

	///	rhsImag.zero();
	///	model.RHS.zero();
		
			Ks=new SpMatComp(model.Hs.addSmallerNew(model.Ss.timesNew(pcw2)),model.Ss.timesNew(0));
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
			model.magMat.addRHS_Continuity(model,pcrhs,rt);
					
		//	model.Qs.shownzA();
			int nuned=model.numberOfUnknownEdges;

			
			int [] map=new int[model.numberOfVarNodes];
			
			if(model.numberOfVarNodes>0) {
			for(int i=1;i<=model.numberOfEdges;i++){
			 int nd=model.edge[i].node[0].id;
	
			 int ind=model.nodeVarIndex[nd]-1;
			 if(ind<0) continue;
			 map[ind]=i;
			}
			}
			
				
				//model.Qs.shownzA();;
	//		model.Ps.shownzA();;
				for(int i=0;i<model.numberOfVarNodes;i++){
	
					int indx=model.edgeUnknownIndex[map[i]]-1;		
					
					SpVectComp spv=new SpVectComp(model.Ps.row[i].times(-1),model.RHS.length);
					
			
					for(int k=0;k<spv.nzLength;k++)
					{

						if(spv.index[k]==indx){
							spv.el[k]=spv.el[k].add(new Complex(-rt.el[i+nuned],0).times(pcw));
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
		
		
		if(photonic==0 || photonic==2)
		//	b=new VectComp(model.RHS);
		b=new VectComp(model.RHS,rhsImag);
		else
		 b=new VectComp(model.RHS.aug(new Vect(nb)));
		
		//b.show();
	
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
				//	xr.show();
					
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
				
			//	b.show();

				b =Ks2t.amul(b);

				Ks=Cs.deepCopy();
			//	b.show();
			//	Ks.show();

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
		
		
		//Ks.shownz();

		//Ks.diagSym().show();
		Ks.setSymHerm(1); // symmetric but not Hermitian

		int m=b.length;
		
		
		model.Ci=Ks.scale(b);



		SpMatComp Ls=Ks.ichol(1.05);
		
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
		
		double y0=model.spaceBoundary[2];
		double y1=model.spaceBoundary[3];
		double L=y1-y0;
		
	//	double R0=model.R0;
	//	double T0=model.T0;
		

		
		for(int i=1;i<=model.numberOfEdges;i++){		
			double y=model.edge[i].node[0].getCoord(1)-y0;
	
			
			if(model.edge[i].map>0) continue;
			
			double E0r=cos(pcw*y);//*(1-1*y/L);
			double E0m=sin(pcw*y);//*(1-1*y/L);
			
/*			double R0r= R0*cos(-pcw*y)*(1-y/L);
			double R0m=R0*sin(-pcw*y)*(1-y/L);
			
			double T0r=T0*cos(pcw*(y-L))*(y/L);
			double T0m=T0*sin(pcw*(y-L))*(y/L);*/
			
			
			int indx=model.edgeUnknownIndex[i]-1;		
			if(indx>=0){
				xc.el[indx].re+=E0r;//+R0r+T0r;
				xc.el[indx].im+=E0m;//+R0m+T0m;
			}
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
		av2.show();
		

			
	if(cc1>0) av1=av1.times(1./cc1);
		
	if(cc2>0) av2=av2.times(1./cc2);
	
/*	double RR=Math.pow(av1.norm(),2);;
	double TT=av2.norm2();
	double sum2=RR;//+TT;

	double factor2=model.E0*model.E0/sum2;
	double factor=Math.sqrt(factor2);
	
	xc.timesVoid(factor);
	av2=av2.times(factor);
	av1=av1.times(factor);
		*/
	//	av.show();
		double pcw2=pcw*pcw;
		
		T.el[kf]=av2.norm2();///(av1.norm2());


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
