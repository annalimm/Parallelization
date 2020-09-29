#include<cstdio>
#include<cstdlib>
#include<cctype>
#include<cmath>
#include<ctime>
#include <sys/time.h>
#include <iostream>
#include <memory>
#include <unistd.h>
#include <thread>
#include <signal.h>

#include "1.hpp"
#define MAX 2147483647
using namespace std;


double get_full_time(){
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)(t.tv_sec+(t.tv_usec)/1000000.0);
}


void* gauss(void* q){
    int d;

    ARGS_gauss* pars=(ARGS_gauss*)q;
    //printf(" first: %d / %d\n",pars->thread_num,pars->p);
    synchronize(pars->p);


    pars->time=get_full_time();
    //printf(" first\n");
    //printf(" 111111111 %d is\n",pars->thread_num);
    d=Solve(pars->n,pars->a,pars->b,pars->x,pars->ind,pars->thread_num,pars->p,pars->er);
    synchronize(pars->p);
    pars->time=get_full_time()-pars->time;
    //printf(" 111111111 %d is\n",pars->thread_num);
    d+=0;

    return 0;
}

void* newyaz(void* q){

    ARGSS_newyaz* pars=(ARGSS_newyaz*)q;
    synchronize(pars->p);
    //printf(" thread_num is");
    //printf(" thread_num is %d \n",thread_num);
    pars->time=get_full_time();
    pars->sum_4tread =newyaz(pars->n,pars->a0,pars->b0,pars->x,pars->thread_num/*,pars->p*/,pars->sum_4tread);
    synchronize(pars->p);
    pars->time=get_full_time()-pars->time;
    return 0;
}
/*
void* norm(void* q){

    ARGSS_norm* pars=(ARGSS_norm*)q;
    synchronize(pars->p);
    //printf(" thread_num is");
    //printf(" thread_num is %d \n",thread_num);
    pars->time=get_full_time();
    pars->sum_4tread =norm_paral(pars->n,pars->a0,pars->thread_num,pars->p,pars->sum_4tread);
    synchronize(pars->p);
    pars->time=get_full_time()-pars->time;
    return 0;
}*/



int main(int argc, char* argv[]){
    int n, p;

    if ((argc!=3)&&(argc!=4)){
        cout<<"error program initialization";
        return -1;
    }else{
    if (sscanf(argv[1],"%d",&p)!=1){
        cout<<"wrong thread format";
        return -1;
    }
    if (p<=0){
        cout<<"thread value error";
        return -1;
    }
    if (sscanf(argv[2],"%d",&n)!=1){
        cout<<"wrong matrix size n";
        return -1;
    }
    if (n<=0){
        cout<<"n must be positive";
        return -1;
    }

    if((double)n > sqrt((double)MAX))
    {
      cout<<"too big size of matrix"<<endl;
      return -1;
    }
    double* a=new double[n*n];
	double* a0=new double[n*n];
	double* b=new double [n];
	double* b0=new double [n];
	double* x=new double [n];
    int* ind=new int [n];
    
    pthread_t* threads_gauss = new pthread_t[p];
    ARGS_gauss* argue_gauss = new ARGS_gauss[p];

    pthread_t* threads_newyaz = new pthread_t[p];
    ARGSS_newyaz* argue_newyaz = new ARGSS_newyaz[p];

   // pthread_t* threads_norm = new pthread_t[p];
   // ARGSS_norm* argue_norm = new ARGSS_norm[p];


	if (!a || !a0)
	{
	    cout<<"no memory"<<endl;
	    return -1;
	}
	if (!b || !b0)
	{
	    cout<<"no memory"<<endl;
	    return -1;
	}
	if (!x)
	{
	    cout<<"no memory"<<endl;
	    return -1;
	}
        if (!threads_gauss || !threads_newyaz /*|| !threads_norm*/)
	{
	    cout<<"no memory"<<endl;
	    return -1;
	}
        if (!argue_gauss || !argue_newyaz /*|| !argue_norm*/)
	{
	    cout<<"no memory"<<endl;
	    return -1;
	}

    if (argc==3)
    {
        matrix_input_formula(n,a,b);
    }
    else if (argc==4)
    {
      if ((matrix_input_file(n,a,b,argv[3]))==-1)
      {
          cout<<"matrix data error"<<endl;
          return -1;
      }
    }
    
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            a0[i*n+j]=a[i*n+j];   
}}
            for(int i=0;i<n;i++)
                b0[i]=b[i];
 

    int *er=new int[1];
    if(!er){
        cout<<"no mem";
        return -1;
    }

    er[0]=0;
 
    
    for (int i=0;i<p;i++){
        argue_gauss[i].n=n;
        argue_gauss[i].a=a;
        argue_gauss[i].b=b;
        argue_gauss[i].x=x;
        argue_gauss[i].ind=ind;
        argue_gauss[i].thread_num=i;
        argue_gauss[i].p=p;
        argue_gauss[i].er=er;

    }

    for (int i=0;i<p;i++){
        argue_newyaz[i].n=n;
        argue_newyaz[i].a0=a0;
        argue_newyaz[i].b0=b0;
        argue_newyaz[i].x=x;
        argue_newyaz[i].ind=ind;
        argue_newyaz[i].thread_num=i;
        argue_newyaz[i].sum_4tread=0.0;
        argue_newyaz[i].p=p;
        argue_newyaz[i].er=er;

    }
 
    for (int i=0;i<p;i++){
        if (pthread_create(threads_gauss+i,0,gauss,argue_gauss+i)){
            cout<<"can't create thread"<< i<<endl;
            delete[] a;
            delete[] b;
            delete[] a0;
            delete[] b0;
            delete[] x;
            delete[] threads_newyaz;
            delete[] threads_gauss;
            delete[] argue_newyaz;
            delete[] argue_gauss;
            delete[] er;
                        delete[] ind;
            return -1;
        }
    }

 
    for (int i=0;i<p;i++){
        if (pthread_join(threads_gauss[i],0)){
            cout<<"can't wait thread"<< i<<endl;

            delete[] a;
            delete[] b;
            delete[] a0;
            delete[] b0;
            delete[] x;
            delete[] threads_newyaz;
            delete[] threads_gauss;
            delete[] argue_newyaz;
            delete[] argue_gauss;
            delete[] er;
                                    delete[] ind;
            return -1;
        }
    }
 
    if (er[0]!=0){
        cout<<"singular matrix"<<endl;
        delete[] a;
        delete[] b;
        delete[] a0;
        delete[] b0;
        delete[] x;
        delete[] threads_newyaz;
        delete[] threads_gauss;
        delete[] argue_newyaz;
        delete[] argue_gauss;
        delete[] er;
                                delete[] ind;
        return -1;
    }


    cout<<"Время гаусс: "<<argue_gauss[0].time<<endl;

    cout<<endl;


 
    for (int i=0;i<p;i++){
        if (pthread_create(threads_newyaz+i,0,newyaz,argue_newyaz+i)){
            cout<<"can't create thread["<< i<<"]"<<endl;
            delete[] a;
            delete[] b;
            delete[] a0;
            delete[] b0;
            delete[] x;
            delete[] threads_newyaz;
            delete[] threads_gauss;
            delete[] argue_newyaz;
            delete[] argue_gauss;
            delete[] er;
                                    delete[] ind;
            return -1;
        }
    }

 
    for (int i=0;i<p;i++){
        if (pthread_join(threads_newyaz[i],0)){
            cout<<"can't wait thread"<< i<<endl;

            delete[] a;
            delete[] b;
            delete[] a0;
            delete[] b0;
            delete[] x;
            delete[] threads_newyaz;
            delete[] threads_gauss;
            delete[] argue_newyaz;
            delete[] argue_gauss;
            delete[] er;
                                    delete[] ind;
            return -1;
        }
    }
 

  double sum=0.0;
  for(int i=0;i<p;i++)
      sum+=argue_newyaz[i].sum_4tread;


    cout<<"\nВремя  невязка: "<<argue_newyaz[0].time<<endl;
   printf(" Невязка: %le\n", sqrt(sum));
 
   //printf("previous Невязка: %le\n",newyaz_nepar(n,a0,b0,x));




/*

    for (int i=0;i<p;i++){
        argue_newyaz[i].n=n;
        argue_newyaz[i].a0=a0;
        argue_newyaz[i].b0=b0;
        argue_newyaz[i].x=x;
        argue_newyaz[i].ind=ind;
        argue_newyaz[i].thread_num=i;
        argue_newyaz[i].sum_4tread=0.0;
        argue_newyaz[i].p=p;
        argue_newyaz[i].er=er;

    }


    for (int i=0;i<p;i++){
        if (pthread_create(threads_norm+i,0,norm,argue_norm+i)){
            cout<<"can't create thread"<< i<<endl;
            return -1;
        }
    }

 
    for (int i=0;i<p;i++){
        if (pthread_join(threads_norm[i],0)){
            cout<<"can't wait thread"<< i<<endl;

            delete[] a;
            delete[] b;
            delete[] a0;
            delete[] b0;
            delete[] x;
            delete[] threads_newyaz;
            delete[] threads_gauss;
                        delete[] threads_norm;
            delete[] argue_newyaz;
                        delete[] argue_norm;
            delete[] argue_gauss;
            delete[] er;

            
            return -1;
        }
    }
 
    if (er[0]!=0){
        cout<<"singular matrix"<<endl;
        delete[] a;
        delete[] b;
        delete[] a0;
        delete[] b0;
        delete[] x;
        delete[] threads_newyaz;
        delete[] threads_gauss;
                    delete[] threads_norm;
        delete[] argue_newyaz;
                    delete[] argue_norm;
        delete[] argue_gauss;
        delete[] er;
            
        return -1;
    }

  double max=0.0;
  for(int i=0;i<p;i++)
      if(argue_norm[i].sum_4tread> max)
        max = argue_norm[i].sum_4tread;


    cout<<"\nВремя  норма: "<<argue_norm[0].time<<endl;
   printf("Невязка: %le\n", max);
 
   printf("previous Ннорма: %le\n",norm(n,a0));

*/
   delete[] a;
   delete[] b;
   delete[] a0;
   delete[] b0;
   delete[] x;
   delete[] threads_newyaz;
   delete[] threads_gauss;
   delete[] argue_newyaz;
   delete[] argue_gauss;
   delete[] er;
   delete[] ind;
}
    return 0;
}

