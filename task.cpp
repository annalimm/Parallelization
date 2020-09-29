#include<cstdio>
#include<cstdlib>
#include<cctype>
#include<cmath>
#include<ctime>
#include<math.h>
#include <sys/time.h>
#include <iostream>
#include <memory>
#include <unistd.h>
#include <pthread.h>

#include "1.hpp"
using namespace std;


void synchronize(int total_threads){
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;

    pthread_mutex_lock(&mutex);

    threads_in++;
    if (threads_in >= total_threads){
            threads_out = 0;
            pthread_cond_broadcast(&condvar_in);
    } else
            while (threads_in < total_threads)
                    pthread_cond_wait(&condvar_in,&mutex);

    threads_out++;
    if (threads_out >= total_threads){
            threads_in = 0;
            pthread_cond_broadcast(&condvar_out);
    } else
            while (threads_out < total_threads)
                    pthread_cond_wait(&condvar_out,&mutex);

    pthread_mutex_unlock(&mutex);
}



int Solve(int n, double* m, double* b, double* x, int* ind, int thread_num, int p, int* er){
    double elmax;
    int indmax;
    int s;

    double buf;
    double norma;

    for (int i = 0; i < n; i++)
        ind[i] = i;

    norma = norm(n,m);

    for (int k = 0; k < n; k++){
        synchronize(p);                     //собираем в одной точке

        if (thread_num == 0){
            indmax = k;
            elmax = fabs(m[k*n+k]);

            for (int j = k + 1; j < n; j++)
               if (fabs(m[k*n+j]) > elmax){
                    elmax = fabs(m[k*n+j]);
                    indmax = j;
               }

            s = ind[k];
            ind[k] = ind[indmax];
            ind[indmax] = s;

            for (int i = 0; i < n; i++){
                buf = m[i*n+k];
                m[i*n+k] = m[i*n+indmax];
                m[i*n+indmax] = buf;
            }


            if (fabs(m[k*n+k]) <= (2e-16)*norma){
                er[0] = 1;
            }
            else
            {
              buf = 1.0/m[k*n+k];

              for (int j = k; j<n; j++){
                m[k*n+j] *= buf;
              }
              b[k]*=buf;

            }
         }//нашли коэффициенты для всех

        synchronize(p);

        if (er[0] == 1)
            return -1;


        //for (int i = thread_num*4; i < thread_num*4 ; i += 1){    //распараллеливание
        for (int i = k + 1 + thread_num; i < n; i += p){    //распараллеливание
        //k + 1 + thread_num + p
            buf = m[i*n+k];
            for (int j = k; j<n; j++){
                m[i*n+j] -= m[k*n+j]*buf;
            }

            b[i] -= b[k]*buf;
        }

    }
        synchronize(p);

        if (thread_num == 0){
            for (int i = n - 1; i >= 0; i--){
                for (int j = i - 1; j >= 0; j--) {
                    b[j] = b[j] - b[i] * m[j * n + i] / m[i * n + i];
                    m[j * n + i] = 0;
                }
            b[i] = b[i] / m[i * n + i];
            m[i * n + i] = 1;
            }

            for (int i = 0; i < n; i++) {
                x[ind[i]] = b[i];
            }
        }

        return 0;

}



double newyaz(int n, double*a, double *b, double* x, int thread_num/*, int p*/, double sum_4tread){
    double c;
    double ans;
    ans=0;

    if (thread_num== 0){
        for (int i = thread_num; i < 4 ; i += 1){
            //берём блоки по p тредов
           // for (int i = thread_num; i < n; i += p){
                c=0.0;
                for (int j=0;j<n;j++){
                    c+=a[i*n+j]*x[j];
                }
                c-=b[i];
                //printf("c %d %le\n",i ,c);
                ans+=c*c;
            }
    }else{
for (int i = (thread_num*4); i < thread_num*4+4 ; i += 1){
   // for (int i = thread_num; i < n; i += p){
        c=0.0;
        for (int j=0;j<n;j++){
            c+=a[i*n+j]*x[j];
        }
        c-=b[i];
        //printf("c %d %le\n",i ,c);
        ans+=c*c;
    }}
    sum_4tread = ans;
   // printf("sum_num %le\n",ans);
    return sum_4tread; //.. от каждого треда получить число и просуммировать в нулевой
}






double newyaz_nepar(int n, double*a, double *b, double* x){
    double c;
    double ans;
    ans=0;
    for (int i=0;i<n;i++)
    {
        c=0.0;
        for (int j=0;j<n;j++)
        {
            c+=a[i*n+j]*x[j];
        }
        c-=b[i];
        //printf("c %d %le\n",i ,c);
        ans+=c*c;
    }
    return sqrt(ans);
}


double err(int n, double *x){
    double ans=0;

    for (int i=0;i<n;i++){
        if (!(i%2)) ans+=(x[i]-1)*(x[i]-1);
        else ans+=x[i]*x[i];
    }
    return sqrt(ans);
}

double norm(int n, double *m)
{
    double sum=0;
    double maximum=0;
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            sum+=fabs(m[i*n+j]);
        }
        if (sum>maximum)
        {
            maximum=sum;
        }
        sum=0;
    }

    return maximum;
}

double norm_paral(int n, double *m, int p,  int thread_num,double sum_4tread)
{
    double sum=0;
    double maximum=0;
    for (int i = thread_num; i < n; i += p)
    {
        for (int j=0;j<n;j++)
        {
            sum+=fabs(m[i*n+j]);
        }
        if (sum>maximum)
        {
            maximum=sum;
        }
        sum=0;
    }
    sum_4tread = maximum;
    return sum_4tread;
}

int  matrix_input_file(int k,double *m, double* b, char* filestr)
{
    double buf;
    FILE* file;

    file=fopen(filestr,"rt");
    if(!file)
    {
        cout<<"file open error";
        return -1;
    }          

    for (int i=0;i<k;i++)
    {

        buf=0;
        b[i]=0;
        for (int j=0;j<k;j++)
        {
            if (fscanf(file,"%lf",&m[i*k+j])!=1)
            {
                fclose(file);
                return -1;
            }
            if (!(j%2)) buf+=m[i*k+j];   
        }

        b[i]=buf;
    }
    fclose(file);
    return 0;
}

void solution_print(int n, double *a, double *b, double* x, int ind)      //int m для вывода первых m чисел решения
{
//    if (m>n) m=n; печатать только первые m элементов решения...

  int i,j;

  if(ind == 1)
  {
    cout<<"Матрица и правая часть:"<<endl;
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
        cout<<a[i*n+j]<<" ";
      cout<<b[i]<<endl;
    }
  }


  if(ind == 2)
  {
    cout<<"Решение: (";
    for (i=0;i<n-1;i++)
        cout<<x[i]<<", ";
    cout<<x[n-1]<<")"<<endl;
  }

}


double f(int i, int j,int n){
    double ii, jj;
double res = 0;

    for (i=0;i<=(n - n/3);i++){
        for ( j=0;j<=(n/3);j++){
    ii = (double)i;
    jj = (double)j;
     res = fabs(ii-jj);}}


    for (int i= n - 1;i<=(n - n/3 + 1);i++)
//for (int i=(n - n/3 + 1);i<n;i++)

    { int k = 0;  res = 0;
    for (int j=(n/3+1);j<n;j++)
    {int l = 0;
        if (k == l){
         res=7;}else{ res = 0;}
        l++;
        //if (!(j%2)) buf+=a[i*n+j];
    }
    k++;

}

    return res;
}


/*
void matrix_input_formula(int n, double* a,double* b)
{
    double buf=0.0;
    for (int i=0;i<=(n - n/3);i++)
    {
        buf=0.0;
        for (int j=0;j<=(n/3);j++)
        {
            a[i*n+j]=f(i,j);
            //if (!(j%2)) buf+=a[i*n+j];
        }
}

        for (int i= n - 1;i<=(n - n/3 + 1);i++)
    //for (int i=(n - n/3 + 1);i<n;i++)

        { int k = 0;
        for (int j=(n/3+1);j<n;j++)
        {int l = 0;
            if (k == l){
            a[i*n+j]=7;}
            l++;
            //if (!(j%2)) buf+=a[i*n+j];
        }
        k++;


    }

        for (int i = 0; i < n;i++){
            buf=0.0;
            for (int j = 0; j < n;i++){
                if (!(j%2)) {buf+=a[i*n+j];}

                //printf("%f ",a[i*n+j]);
            }   b[i]=buf;
            printf("/n" );
        }

}*/

void matrix_input_formula(int n, double* a,double* b)
{
    double buf=0.0;
    for (int i=0;i<n;i++)
    {
        buf=0.0;
        for (int j=0;j<n;j++)
        {
            a[i*n+j]=f(i,j,n);
            if (!(j%2)) buf+=a[i*n+j];
        }
        b[i]=buf;
    }



    for (int i = 0; i < n;i++){

        for (int j = 0; j < n;i++){


            printf("%f ",a[i*n+j]);
        }
        printf("/n" );
    }
}






