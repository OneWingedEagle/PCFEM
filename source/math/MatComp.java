package math;

import static java.lang.Math.*;

import java.awt.Color;

import io.Loader;

import java.util.Arrays;
import java.util.Random;

import javax.swing.JFrame;

import org.math.plot.Plot2DPanel;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.QRDecomposition;

public class MatComp extends Mat {
	public Complex[][] el=null;
	public int nRow;
	public int nCol;
	public static void main(String[] Args){
		
	}



	public MatComp(){};


	public  MatComp(Complex[][] array){
		nRow=array.length;
		nCol=array[0].length;
		el=array;
	}

	public  MatComp(int nRow, int nCol){
		this.nRow=nRow;
		this.nCol=nCol;
		this.el=new Complex[nRow][nCol];
	}

	public  MatComp(int[] dim){
		this.nRow=dim[0];
		this.nCol=dim[1];
		this.el=new Complex[dim[0]][dim[1]];
	}

	public  MatComp(int nCol){
		this.nRow=1;
		this.nCol=nCol;
		this.el=new Complex[nRow][nCol];
	}



	public void set(Complex[][] A){
		this.el=A;
		this.nRow=A.length;;
		this.nCol=A[0].length;
	}

	public void set(Complex[] v){
		int J;
		J=v.length;
		this.el[0]=v;
		this.nRow=1;
		this.nCol=J;
	}

	public int[] size(){
		int[] size=new int[2];
		size[0]=this.nRow;
		size[1]=this.nCol;
		return size;
	}

	public void sz(){
		util.show(this.size());
	}


	public MatComp deepCopy(){
		int I=this.nRow;
		int J=this.nCol;
		MatComp A=new MatComp(I,J);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				A.el[i][j]=this.el[i][j];
		return A;
	}	




	public void ones(int I,int J){
		el=new Complex[I][J];
		this.nRow=I;
		this.nCol=J;
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				el[i][j].re=1;

	}

	public void eye(int I){
		this.nRow=I;
		this.nCol=I;
		el=new Complex[I][I];
		for(int i=0;i<I;i++)
			el[i][i].re=1;	
	}

	public MatComp add(MatComp A){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;

		if(I1!=I2 || J1!=J2) throw new IllegalArgumentException("Array dimensions do not agree");

		MatComp C=new MatComp(I1,J1);

		for(int i=0;i<I1;i++)
			for(int j=0;j<J1;j++)
				C.el[i][j]=this.el[i][j].add(A.el[i][j]);
		return C;
	}


	public MatComp sub(MatComp A){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;

		if(I1!=I2 || J1!=J2) throw new IllegalArgumentException("Array dimensions do not agree");

		MatComp C=new MatComp(I1,J1);

		for(int i=0;i<I1;i++)
			for(int j=0;j<J1;j++)
				C.el[i][j]=this.el[i][j].sub(A.el[i][j]);
		return C;
	}
	


	public MatComp mul(MatComp A){

		return mul(A,A.nCol);
	}

	public MatComp mul(MatComp A,int K){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;
		if(J1!=I2 || K>J2) throw new IllegalArgumentException("Matrix dimensions do not agree");
		Complex s;
		MatComp C=new MatComp(I1,K);
		for(int i=0;i<I1;i++)
			for(int j=0;j<K;j++){
				s=new Complex(0,0);
				for(int n=0;n<J1;n++)
				
					s=s.add(this.el[i][n].times(A.el[n][j]));
				C.el[i][j]=s;
			}
		return C;						
	}

	public VectComp mul(VectComp u){
		int I,J,L;
		I=this.nRow;
		J=this.nCol;
		L=u.length;

		if(J!=L) throw new IllegalArgumentException("Matrix dimensions do not agree");
		Complex s;
		VectComp v=new VectComp(I);
		for(int i=0;i<I;i++){		
			s=new Complex(0,0);
			for(int n=0;n<L;n++)
				//if(this.el[i][n]!=0 && u.el[n]!=0)
				s=s.add(this.el[i][n].times(u.el[n]));
				//s=s+this.el[i][n]*u.el[n];
			v.el[i]=s;
		}
		return v;						
	}
	

	public VectComp getColVectComp(int j){

		VectComp v=new VectComp(this.nRow);
		for(int i=0;i<this.nRow;i++){			
			v.el[i]=this.el[i][j];
		}
		return v;						
	}



	public VectComp rowVectComp(int i){


		return new VectComp(this.el[i]);						
	}

}
