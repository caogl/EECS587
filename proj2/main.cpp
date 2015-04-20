#include<iostream>
#include<float.h>
#include<cmath>
#include<omp.h>
#include<stack>
#include<utility>
using namespace std;

stack<pair<double, double> > stack1;
double epsilon=0.000001;
int s=12;
double sum(double x);
double M=-DBL_MAX;
int activeThread=0;

omp_lock_t maxLock;
omp_lock_t stackLock;

int main()
{
        stack1.push(make_pair(1, 100));
        double cpu1=omp_get_wtime();
        omp_init_lock(&maxLock);
        omp_init_lock(&stackLock);

        #pragma omp parallel 
        {
                while(true)
                {
                        omp_set_lock(&stackLock);
                        //cout<<"thread "<<omp_get_thread_num()<<" is running"<<endl;
                        if(stack1.empty())
                        {
                                if(activeThread==0)
                                {
                                        omp_unset_lock(&stackLock);
                                        break;
                                }
                                else
                                        omp_unset_lock(&stackLock);
                        }
                        else
                        {
                                activeThread++;
                                double c=stack1.top().first;
                                double d=stack1.top().second;
                                //cout<<"Now c is: "<<c<<" and d is: "<<d<<endl;
                                //cout<<"Now M is: "<<M<<endl;
                                stack1.pop();
                                omp_unset_lock(&stackLock);

                                double f_c=sum(c);
                                double f_d=sum(d);

                                omp_set_lock(&maxLock);
                                if(f_c>M || f_d>M)
                                {
                                        M=max(f_c, f_d);
                                }
                                if(((f_c+f_d+s*(d-c))/2)>M+epsilon) // split the interval
                                {
                                        omp_set_lock(&stackLock);
                                        stack1.push(make_pair(c, (c+d)/2));
                                        stack1.push(make_pair((c+d)/2, d));
                                        activeThread--;
                                        omp_unset_lock(&stackLock);
                                }
                                else
                                {
                                        omp_set_lock(&stackLock);
                                        activeThread--;
                                        omp_unset_lock(&stackLock);
                                }
                                omp_unset_lock(&maxLock);
                        }
                }
        }
        double cpu2=omp_get_wtime();
        double time=cpu2-cpu1;
        cout<<"The runtime is: "<<time<<endl;
        cout<<"The final result is: "<<M<<endl;
        return 0;
}

double sum(double x)
{
        double result=0;
        for(int i=100; i>=1; i--)
        {
                double tmp=0;
                for(int j=i; j>=1; j--)
                {
                        tmp+=pow(x+j, -3.1);
                }
                result+=(sin(x+tmp)/pow(1.2, i));
        }
        return result;
}

