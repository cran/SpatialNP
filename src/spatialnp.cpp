#include <cmath>
//#include <R.h>
//#include <Rdefines.h>
//#include <Rinternals.h>


extern "C" {

using namespace std;

  double **prepmat(double *X, int n, int k)
  //be careful to free this memory!!!
  {
    int i;
    int j;
    double **Y = new double* [n];
    for (i=0; i<n; i++) Y[i]=new double [k];
    for (i=0;i<n; i++) 
      for (j=0;j<k;j++)
    Y[i][j]=X[j*n+i];
    return Y;
  }

  double *matrix_prod1(double *X, double *Y, int n, int p, int m)
  // n*p x p*m 
  {
    int i;
    int j;
    int k;
    int l=0;
    
    double *Z = new double [n*m];  
 
    for(i=0;i<n*m;i++){
       Z[i] = 0;
    } 
   
    for(j=0;j<m;j++){
     for(i=0;i<n;i++){ 
      for(k=0;k<p;k++){
       Z[l] += X[k*n+i]*Y[j*p+k];
      }
     l++;
     }
    } 
    return Z;
  }

  double *row_sums(double *X, int n, int k)
  {
    int i;
    int j;
    double *Z = new double [n];

    for(i=0;i<n;i++){
     Z[i]=0;
     for(j=0;j<k;j++){ 
      Z[i] += X[j*n+i];
     }
    } 
    return Z;
  }

  double *pair_diff_inc(double *X, int n, int k, int num)
  {
    int i;
    int j;
    int m; 
    int l=0; 
    double *Z = new double [((n-num)*num+num*(num-1)/2)*k];
    for(m=0;m<k;m++){ 
     for(i=0;i<(n-num);i++){
      for(j=(i+1);j<(i+num+1);j++){
       Z[l] = X[m*n+i]-X[m*n+j];
       l ++;
      }
     }
     for(i=(n-num);i<(n-1);i++){
      for(j=(i+1);j<n;j++){
       Z[l] = X[m*n+i]-X[m*n+j];
       l ++;
      }
     }
    } 
    return Z;
  }
 
    void pairdiffc(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int e;
    int r=0;
    for (i=0; i<(n-1); i++)
      for (j=(i+1); j<n; j++)
    for (e=0;e<k;e++) {
      result[r]=X[e*n+i]-X[e*n+j];
      r++;}
  }
  
  void pairprod(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int e;
    int r=0;
    for (i=0; i<(n-1); i++)
      for (j=(i+1); j<n; j++)
    for (e=0;e<k;e++) {
      result[r]=X[e*n+i]*X[e*n+j];
      r++;}
  }

  void pairsumc(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int e;
    int r=0;
    for (i=0; i<(n-1); i++)
      for (j=(i+1); j<n; j++)
    for (e=0;e<k;e++) {
      result[r]=X[e*n+i]+X[e*n+j];
      r++;}
  }

  void norming(double *X, int *nk, double *result)
  {
    int i;
    int j;
    int n=nk[0]; 
    int k=nk[1]; 
    for (i=0; i<n; i++)
      result[i]=0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<k; j++) 
    result[i]+=X[j*n+i]*X[j*n+i];
      result[i]=sqrt(result[i]);
    }
  }

  void sum_of_sign_outers(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int m;  
    int p=0;
    double *r; 
    r = new double[n];

    double **signs=prepmat(X,n,k);
    //compute norms:
    for (i=0; i<n; i++)
      r[i]=0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<k; j++) 
	r[i]+=(signs[i][j]*signs[i][j]);
    r[i]=sqrt(r[i]);
    }
    
    //compute signs:
    for(i=0; i<n; i++) 
      for(j=0; j<k; j++)
	signs[i][j]=signs[i][j]/r[i];

   //prepare the result:
    for (i=0;i<(k*k);i++) result[i]=0.0;

    //compute the sum of outer products:
    for(j=0;j<k;j++)
      for(m=0;m<k;m++)
	{ 
	  for(i=0;i<n;i++)
	    result[p]+= (signs[i][j]*signs[i][m]);
	  p++;
	}
    for(i=0;i<n;i++)   
     delete [] signs[i]; 
    delete [] signs;
    delete [] r;
  }

  void sum_of_diff_sign_outers(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int m;  
    int p=0;
    double r; 
    double *u;
    u = new double[k];
    double **Y = new double* [k];
    for (i=0; i<k; i++) Y[i]=new double [k];

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
	Y[i][j]=0.0;

    for(i=0;i<(n-1);i++)
      for(j=(i+1);j<n;j++) {
	//go over the pairs
	r=0;
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = X[m*n+i]-X[m*n+j];
	  r+=(u[m]*u[m]);
	}
	r=sqrt(r);
	for(m=0;m<k;m++)
	  //compose the sign
	  u[m]=u[m]/r;
	for(m=0;m<k;m++)
	  for(p=0;p<k;p++)
	    if(p<(m+1)) 
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=(u[m]*u[p]);
      } // j, i

    for(m=0;m<(k-1);m++)
      for(p=(m+1);p<k;p++)
	//fill up the matrix
	Y[m][p]=Y[p][m];

    i=0;
    for(m=0;m<k;m++)
      for(p=0;p<k;p++) {
	//put to the result  
	result[i]=Y[m][p];
	i++;
      }
  
   
    delete [] u;
    for(i=0;i<k;i++)
      delete [] Y[i];
    delete [] Y;
  }

   void sum_of_rank_outers(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int m;  
    int p=0;
    double r; 
    double *uij;
    double *Ri;
    uij = new double[k];
    Ri = new double[k];
    double **Y = new double* [k];
    for (i=0; i<k; i++) Y[i]=new double [k];

    for(i=0;i<k;i++)
     for(j=0;j<k;j++)
      Y[i][j]=0.0;

   for(i=0;i<n;i++) {
    for(m=0;m<k;m++)
     Ri[m]=0.0;
    //go over the ranks
    for(j=0;j<n;j++) {
     if(j!=i) {
      r=0.0;
      // compute the difference and its squared norm
      for(m=0;m<k;m++) {
       uij[m] = X[m*n+i]-X[m*n+j];
       r+=(uij[m]*uij[m]);
      }
      // compute the sign
      // and add to the sum
      r=sqrt(r);
      for(m=0;m<k;m++)
       Ri[m]+=(uij[m]/r);
     } //if
    } // j
    // divide by n
    for(m=0;m<k;m++)
     Ri[m]=Ri[m]/n;
    for(m=0;m<k;m++)
     for(p=0;p<k;p++)
      if(p<(m+1)) 
       //compute the outer product elements
       //and sum to the sum matrix
       Y[m][p]+=(Ri[m]*Ri[p]);
   } // i
 
   for(m=0;m<(k-1);m++)
    for(p=(m+1);p<k;p++)
     //fill up the matrix
      Y[m][p]=Y[p][m];

   i=0;
   for(m=0;m<k;m++)
    for(p=0;p<k;p++) {
     //put to the result  
    result[i]=Y[m][p];
    i++;
    }

    
   delete [] uij;
   delete [] Ri;
   for(i=0;i<k;i++)
     delete [] Y[i];
   delete [] Y;
  }

  void spatial_ranks(double *X, int *nk, double *result)
{
  int n=nk[0]; 
  int k=nk[1]; 
  int i;
  int j;
  int m;
  double d; 
  double **data=prepmat(X,n,k);

  double **ranks = new double * [n];
    for (i=0; i<n; i++) ranks[i]=new double [k];

  for(i=0;i<n;i++)
    for(m=0;m<k;m++)
      ranks[i][m]=0.0;

  double *temp = new double[k];

  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      if(i!=j) {
    //compute the difference and its norm:
    for(m=0;m<k;m++) 
      temp[m]=data[i][m]-data[j][m];
    d=0.0;
    for(m=0;m<k;m++)
      d+=(temp[m]*temp[m]);
    d=sqrt(d);
    //add to the sum of signs:
    for(m=0;m<k;m++)
      ranks[i][m]+=(temp[m]/d);
      }
    }
  }

  j=0;
  for(i=0;i<n;i++)
    for(m=0;m<k;m++) {
      result[j]=ranks[i][m]/n;
      j++;
    }

  for(i=0;i<n;i++) {
    delete [] data[i];
    delete [] ranks[i];
  } 
  delete [] data;
  delete [] ranks;
  delete [] temp;
}


  void signed_ranks(double *X, int *nk, double *result)
{
  int n=nk[0]; 
  int k=nk[1]; 
  int i;
  int j;
  int m;
  double dm=1.0; // m as in minus
  double dp=1.0; // p as in plus
  double **data=prepmat(X,n,k);

  double **ranks = new double * [n];
    for (i=0; i<n; i++) ranks[i]=new double [k];

  for(i=0;i<n;i++){
    for(m=0;m<k;m++)
      ranks[i][m]=0.0;
  }
  double *tempm = new double[k]; // see double dm above
  double *tempp = new double[k];

  for(m=0;m<k;m++)
   tempm[m]=0.0;
 
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      if(i!=j) {
    //compute the difference and its norm:
    for(m=0;m<k;m++) 
      tempm[m]=data[i][m]-data[j][m];
    dm=0.0;
    for(m=0;m<k;m++)
      dm+=(tempm[m]*tempm[m]);
    dm=sqrt(dm);
      }
    //compute the sum and its norm:
    for(m=0;m<k;m++) 
      tempp[m]=data[i][m]+data[j][m];
    dp=0.0;
    for(m=0;m<k;m++)
      dp+=(tempp[m]*tempp[m]);
    dp=sqrt(dp);
    //add to the sum of signs:
    for(m=0;m<k;m++){
      ranks[i][m]+=(tempm[m]/dm+tempp[m]/dp)/2;
      tempm[m]=0;
      tempp[m]=0;
    }
    }
  }

  j=0;
  for(i=0;i<n;i++)
    for(m=0;m<k;m++) {
      result[j]=ranks[i][m]/n;
      j++;
    }

  for(i=0;i<n;i++) {
    delete [] data[i];
    delete [] ranks[i];
  } 
  delete [] data;
  delete [] ranks;
  delete [] tempm;
  delete [] tempp;
}

   bool issame(double *X, int k)
  {
    int i;
    bool result=true;
    for(i=0;i<k;i++){
     if(X[i]!=X[k+i])
      result=false;
    } 
    return result;
  }
 
  void spat_med_step(double *X, int *nk, double *y, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    double sum;
    sum=0;  
    double *Y; 
    Y = new double[2*k];     
    double *r; 
    r = new double[n];
    double *sumsigns;
    sumsigns = new double[k];
    
     for (j=0; j<k; j++)
      sumsigns[j]=0.0;

    double **vectors=prepmat(X,n,k);

     for(i=0; i<n; i++){
	 for(j=0;j<k;j++){
           Y[j]=vectors[i][j];
           Y[k+j]=y[j];
         }  
	 if(issame(Y,k)){
           y[0]+=0.00001;
         } 
       } 

    double **signs = new double* [n];
    for (i=0; i<n; i++) signs[i]=new double [k];
     for (i=0; i<n; i++) 
      for (j=0; j<k; j++) 
	signs[i][j]=vectors[i][j]-y[j];
                            
                                                
    //compute norms:
    for (i=0; i<n; i++)
      r[i]=0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<k; j++) 
	r[i]+=(signs[i][j]*signs[i][j]);
      r[i]=sqrt(r[i]);
    }
    
    //compute signs:
    for(i=0; i<n; i++) 
      for(j=0; j<k; j++)
	signs[i][j]=signs[i][j]/r[i];

    //compute the sum of signs:
  for(i=0;i<n;i++) 
    sum += 1/r[i];  

  for (i=0;i<k;i++)
    result[i]=0.0;

  for(j=0;j<k;j++){
    for(i=0;i<n;i++)
     sumsigns[j]+=signs[i][j];
   result[j]=y[j]+sumsigns[j]/sum;
  }
   for(i=0;i<n;i++)   
     delete [] signs[i]; 
    delete [] signs;
    for(i=0;i<n;i++)   
     delete [] vectors[i];
    delete [] vectors;
    delete [] Y; 
    delete [] r; 
    delete [] sumsigns;
  }

 void center_step(double *X, int *nk, double *y, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    double sum;
    sum=0;   
    double *Y; 
    Y = new double[2*k];
    double *r; 
    r = new double[n];
    double *sumsigns;
    sumsigns = new double[k];
    
     for (j=0; j<k; j++)
      sumsigns[j]=0.0;

    double **vectors=prepmat(X,n,k);

    for(i=0; i<n; i++){
	 for(j=0;j<k;j++){
           Y[j]=vectors[i][j];
           Y[k+j]=y[j];
         }  
	 if(issame(Y,k)){
           y[0]+=0.00001;
         } 
       }

    double **signs = new double* [n];
    for (i=0; i<n; i++) signs[i]=new double [k];
     for (i=0; i<n; i++) 
      for (j=0; j<k; j++) 
	signs[i][j]=vectors[i][j]-y[j];
                                                                      
    //compute norms:
    for (i=0; i<n; i++)
      r[i]=0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<k; j++) 
	r[i]+=(signs[i][j]*signs[i][j]);
      r[i]=sqrt(r[i]);
    }
    
    //compute signs:
    for(i=0; i<n; i++) 
      for(j=0; j<k; j++)
	signs[i][j]=signs[i][j]/r[i];

    //compute the sum of signs:
  for(i=0;i<n;i++) 
    sum += 1/r[i];  

  for (i=0;i<k;i++)
    result[i]=0.0;

  for(j=0;j<k;j++){
    for(i=0;i<n;i++)
     sumsigns[j]+=signs[i][j];
   result[j]=sumsigns[j]/sum;
  }
   for(i=0;i<n;i++)   
    delete [] signs[i]; 
    delete [] signs;
   for(i=0;i<n;i++)   
    delete [] vectors[i];
    delete [] vectors;
    delete [] r; 
    delete [] sumsigns;
    delete [] Y;
  }

  void hl_center_step(double *X, int *nk, double *y, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int m;
    double sum=0;
    double *sr;
    sr = new double[k];
    double r=0; 
    double *sumsignranks;
    sumsignranks = new double[k];
       
    for (i=0; i<k; i++){
     sumsignranks[i]=0.0;
     result[i]=0.0;
    }   
    
    for (i=0; i<n; i++) 
     for (j=i; j<n; j++) {
       for (m=0; m<k; m++) { 
         sr[m]=X[m*n+i]+X[m*n+j]-2*y[m];
        r+=sr[m]*sr[m];
        }
    r=sqrt(r);
        
    for(m=0; m<k; m++)
     sr[m]=sr[m]/r;
 
    sum += 1/r;  
  
    r=0;  
 
    for(m=0;m<k;m++)
     sumsignranks[m]+=sr[m];
    }
    
    for(m=0;m<k;m++)  
     result[m]=0.5*sumsignranks[m]/sum;

   
    delete [] sr;     
    delete [] sumsignranks;
  }


   void sum_of_diff_sign_select(double *X, int *nk, int *num, double *result)
  {
    int n=nk[0]; 
    int k=nk[1];
    int nu=num[0]; 
    int i;
    int j;
    int m;  
    int p=0;
    double r; 
    double *u;
    u = new double[k];
    double **Y = new double* [k];
    for (i=0; i<k; i++) Y[i]=new double [k];

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
	Y[i][j]=0.0;

    double **data=prepmat(X,n,k);

    // n-nu first
    for(i=0;i<(n-nu);i++)
      for(j=1;j<(nu+1);j++) {
	//go over the pairs
	r=0;
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = data[i][m]-data[i+j][m];
	  r+=(u[m]*u[m]);
	}
	r=sqrt(r);
	for(m=0;m<k;m++)
	  //compose the sign
	  u[m]=u[m]/r;
	for(m=0;m<k;m++)
	  for(p=0;p<k;p++)
	    if(p<(m+1)) 
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=(u[m]*u[p]);
      } // j, i
    
    // the rest
    for(i=(n-nu);i<(n-1);i++)
      for(j=(i+1);j<n;j++) {
	//go over the pairs
	r=0;
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = data[i][m]-data[j][m];
	  r+=(u[m]*u[m]);
	}
	r=sqrt(r);
	for(m=0;m<k;m++)
	  //compose the sign
	  u[m]=u[m]/r;
	for(m=0;m<k;m++)
	  for(p=0;p<k;p++)
	    if(p<(m+1)) 
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=(u[m]*u[p]);
      } // j, i


    for(m=0;m<(k-1);m++)
      for(p=(m+1);p<k;p++)
	//fill up the matrix
	Y[m][p]=Y[p][m];

    i=0;
    for(m=0;m<k;m++)
      for(p=0;p<k;p++) {
	//put to the result  
	result[i]=Y[m][p];
	i++;
      }
  
    for(i=0;i<n;i++)   
      delete [] data[i]; 
    delete [] data;
    delete [] u;
    for(i=0;i<k;i++)
      delete [] Y[i];
    delete [] Y;
  }

   

  void symm_huber(double *X, double *V, int *nk, double *cs, double *sigs, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int l;
    int m;  
    int p=0;
    double csq=cs[0];
    double sigsq=sigs[0];
    double rd=0.0; 
    double *u;
    u = new double[k];
    double *r;
    r = new double[k];
    double w;
    double **Y = new double* [k];
    for (i=0; i<k; i++) Y[i]=new double [k];

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
	Y[i][j]=0.0;
 
    for(i=0;i<k;i++)
       r[i]=0.0;

    for(i=0;i<(n-1);i++) 
      for(j=(i+1);j<n;j++) {
	//go over the pairs
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = X[m*n+i]-X[m*n+j];
          //compute the Mahalanobis distance
          for(l=0;l<k;l++) 
           r[l] += u[m]*V[l*k+m];
	}


       for(m=0;m<k;m++)
        rd += r[m]*u[m];

       if(rd<=csq){ w=1/sigsq;
       }else w=(csq/rd)/sigsq;
 
       for(m=0;m<k;m++)
        r[m]=0.0;
       
       rd=0.0;
        
	for(m=0;m<k;m++){
	  for(p=0;p<k;p++){
	    if(p<(m+1)) {
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=w*u[m]*u[p];
            }
          }
        }   
      } // j, i

    for(m=0;m<(k-1);m++)
      for(p=(m+1);p<k;p++)
	//fill up the matrix
	Y[m][p]=Y[p][m];

    i=0;
    for(m=0;m<k;m++)
      for(p=0;p<k;p++) {
	//put to the result  
	result[i]=Y[m][p];
	i++;
      }
  
    delete [] r;
    delete [] u;
    for(i=0;i<k;i++)
      delete [] Y[i];
    delete [] Y;
  }

  

  void outer(double *u, double *v, int k, double *uvT)
  {
    int i;
    int j;

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
    uvT[i*k+j]=u[i]*v[j];
  }

  void outer2(double *u, int k, double *uut)
  {
    int i;
    int j;
  
    for(i=0;i<k;i++)
      for(j=i;j<k;j++) {
        uut[i*k+j]=u[i]*u[j];
        if(j>i)
         uut[j*k+i]=uut[i*k+j];
      }
  }

  void touij(double *xi, double *xj, int k, double *uij)
  {
    int i;
    double s=0.0;
    for(i=0;i<k;i++) {
      uij[i]= xi[i]-xj[i];
      s+=(uij[i]*uij[i]);
    }
    s=sqrt(s);
    for(i=0;i<k;i++) {
      uij[i]=uij[i]/s;
    }
  } 

  void Q2internals(double *X, int *nk, double *result)
  {
    int n=nk[0];
    int k=nk[1];
    int k2=k*k;
    int k4=k2*k2;
    int i;
    int j;
    int m;
    int ii;
    double **data=prepmat(X,n,k);
    double *uijuik = new double [k4];
    double *uij2 = new double [k2];
    double *uik2 = new double [k2];
    double *uij = new double [k];
    double *uik = new double [k];
    double *S2 = new double [k2];
    double *ave = new double [k4];


    for(i=0;i<(k2);i++) S2[i]=0.0;
    for(i=0;i<k4;i++) ave[i]=0.0;
    for (i=0;i<(n-1);i++) {
      for(j=(i+1);j<n;j++) {
        //compute u_ij
        touij(data[i],data[j],k,uij);
        //compute vec(u_ij u_ij^T)
        outer2(uij,k,uij2);
        //sum to S2
        for(ii=0;ii<k2;ii++)
          S2[ii]+=uij2[ii];

      for (m=0;m<n;m++) {
    if(m!=i) {
      touij(data[i],data[m],k,uik);
      outer2(uik,k,uik2);
      outer(uij2,uik2,k2,uijuik);
      for(ii=0;ii<k4;ii++)
        ave[ii]+=uijuik[ii];
    }
    if(m!=j) {
      touij(data[j],data[m],k,uik);
      outer2(uik,k,uik2);
      outer(uij2,uik2,k2,uijuik);
      for(ii=0;ii<k4;ii++)
        ave[ii]+=uijuik[ii];
    }
      }
    }
  }
  //scale S2 correctly, there were n*(n-1)/2 terms to sum
  for(i=0;i<k2;i++) {
    S2[i]=S2[i]/(n*(n-1)/2);
    result[i]=S2[i];
  }
  //scale ave correctly, there were n*(n-1)^2 terms to sum
  for(i=0;i<k4;i++) {
    ave[i]=ave[i]/(n*(n-1)*(n-1));
    result[k2+i]=ave[i];
  }


  for(i=0;i<n;i++) {
    delete [] data[i];
  }
  delete [] data;
  delete [] uij2;
  delete [] uik2;
  delete [] uik;
  delete [] uij;
  delete [] uijuik;
  delete [] S2;
  delete [] ave;
}

  

}
