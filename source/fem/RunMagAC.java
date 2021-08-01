package fem;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.awt.Color;
import java.io.File;
import java.text.DecimalFormat;

import femSolver.StaticElectricSolver;
import main.Main;
import math.Vect;
import math.VectComp;
import math.util;

public class RunMagAC {

	private DecimalFormat formatter=new DecimalFormat("0.00");

	private Vect xp;
	int[] nn=new int[1];

	public static void main(String[] args){
		new Main();
	}

	public void runMag(Model model, Main main){
		

			String folder=model.resultFolder;


			int nTsteps=model.nTsteps;
		
			//Vect T=new Vect(nTsteps);
			Vect[] R_L=new Vect[nTsteps];
			Vect freqs=new Vect(nTsteps);
			int ix=0;
			
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;

		//	if(model.AC) nEnd=nBegin;

			
			
			String currentFolder=model.resultFolder;
			String fluxFolder=model.resultFolder;
			String dispFolder=model.resultFolder;
			
			if(model.saveForce){
		
			
					dispFolder = folder+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveFlux){
					fluxFolder = folder+"\\fluxes";
				
				File dfolder = new File(fluxFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveJe){
				currentFolder = folder+"\\currents";
			
			File dfolder = new File(currentFolder);
			if(dfolder.exists())
				util.deleteDir(dfolder);
			dfolder.mkdir();

		}

		
			
			VectComp xc=new VectComp();

			main.gui.lbX[0].setText("rot ang. ");
			main.gui.lbX[1].setText("Trq. ");


		
			
			model.setMagBC();

			model.solveCoils();
			
			Vect Tr=new Vect(model.nTsteps);


			for(int step=nBegin;step<=nEnd;step+=inc){
				
				if(step>nBegin)
				model.setMagBC();

				double t0=0;
		
				main.gui.tfX[0].setText((step)+"/"+nEnd);

				model.setJ0();	
				
				if(ix>0){
			//	model.freq*=Math.pow(10, 1);	
				model.freq+=.1/5;
				
				}
				freqs.el[ix]=model.freq;
					
							
							if(step==nBegin){
								model.writer.reportData(model);
							}

							xc=model.femSolver.solveMagAC(model, step-nBegin);	
				

							int m=xc.length;
							Vect vr1=new Vect(m);
							Vect vm1=new Vect(m);

							Vect vr2=new Vect(m);
							Vect vm2=new Vect(m);
							for(int k=0;k<m;k++){
								if(k<model.numberOfUnknownEdges){
									vr1.el[k]=xc.el[k].re;
									vm1.el[k]=xc.el[k].im;
									
									vr2.el[k]=-xc.el[k].im*2*PI;
									vm2.el[k]=xc.el[k].re*2*PI;
								}
								else{
									vr2.el[k]=xc.el[k].re;
									vm2.el[k]=xc.el[k].im;
								}
								
							}
							
							double totalLossRe=0;
							double totalLossIm=0;
							double totalEnergyRe=0;
							double totalEnergyIm=0;

						


								model.setSolution(vr2);	
								model.setJe();
								
								if(model.saveJe){
									String JeFile =  currentFolder+"\\JeRe"+step+".txt";
				
									model.writeJe(JeFile);
									
									if(step==nBegin)
										model.writeMesh(currentFolder+"\\bun"+step+".txt");
						
								}
								
								boolean append=false;
								if( step!=nBegin)
									append=true;
								
								totalLossRe=model.writer.outputLoss(model,model.resultFolder+"\\outputs.txt",step,0,append);		
				
								
								model.setSolution(vr1);	
								
								double vrt=0;
								
								double vmt=0;
								
								totalEnergyRe=model.writer.outputEnergies(model,model.resultFolder+"\\outputs.txt",step,-90,true);		


								if(model.saveFlux){
							
									String fluxFile = fluxFolder+"\\fluxRe"+step+".txt";
									model.writeB(fluxFile);
									
									String aFile = fluxFolder+"\\Ar"+step+".txt";
									for(int i=1;i<=model.numberOfEdges;i++){
									//	model.edge[i].node[0].setDeformable(true);
										model.edge[i].node[0].T=model.edge[i].A;
										
										if(model.edge[i].node[0].id==273)
											vrt=model.edge[i].A;
									}
									model.writeNodalScalar(aFile);
									
									
									String forceFile = fluxFolder+"\\ArNodal"+step+".txt";
									
									for(int i=1;i<=model.numberOfEdges;i++){
										model.edge[i].node[0].setDeformable(true);
				
											model.edge[i].node[0].F=new Vect(model.edge[i].A,0);
										}
									model.writeNodalField(forceFile, 1);
									
									if(step==nBegin)
										model.writeMesh(fluxFolder+"\\bun"+step+".txt");
								}

					
								
								model.setSolution(vm2);	
								model.setJe();	

								if(model.saveJe){
										String JeFile =  currentFolder+"\\JeIm"+step+".txt";
						
										model.writeJe(JeFile);
						
								}
								totalLossIm=model.writer.outputLoss(model,model.resultFolder+"\\outputs.txt",step,-90,true);		

								model.setSolution(vm1);	
								
								for(int i=1;i<=model.numberOfEdges;i++){
			
										if(model.edge[i].node[0].id==273)
											vmt=model.edge[i].A;
									}
								
								Tr.el[ix]=new Vect(vrt,vmt).norm();
								
								totalEnergyIm=model.writer.outputEnergies(model,model.resultFolder+"\\outputs.txt",step,0,true);			

								if(model.saveFlux){
						
									String fluxFile = fluxFolder+"\\fluxIm"+step+".txt";
									model.writeB(fluxFile);
									
									String aFile = fluxFolder+"\\Am"+step+".txt";
									for(int i=1;i<=model.numberOfEdges;i++){
									//	model.edge[i].node[0].setDeformable(true);
										model.edge[i].node[0].T=model.edge[i].A;
										
									//	if(model.edge[i].node[0].id==273)
										//	vrt=model.edge[i].A;
									}
									model.writeNodalScalar(aFile);
									
									String forceFile = fluxFolder+"\\AmNodal"+step+".txt";
									
									for(int i=1;i<=model.numberOfEdges;i++){
										model.edge[i].node[0].setDeformable(true);
				
											model.edge[i].node[0].F=new Vect(model.edge[i].A,0);
										}
									model.writeNodalField(forceFile, 1);

								}
								

							
			
						R_L[ix]=new Vect((totalLossRe+totalLossIm),2*(totalEnergyRe+totalEnergyIm));

		

					main.gui.tfX[1].setText(this.formatter.format(model.TrqZ));

					if(model.solver.terminate) break;

					ix++;
				}
			
			model.writer.writeR_L(model,model.resultFolder+"\\outputs.txt",freqs,R_L);		

			util.plot(freqs,Tr);
			Tr.show();
		
				Vect errs=new Vect(model.solver.totalIter);

				for(int i=0;i<errs.length;i++)
					errs.el[i]=model.solver.errs.get(i);
				if(errs.length>0)
				util.plot("error",errs.el,"ICCG Convergence");
			//	util.plot(errs);

	

	}

}
