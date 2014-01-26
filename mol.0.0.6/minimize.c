/*
Copyright (c) 2009-2012, Structural Bioinformatics Laboratory, Boston University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- Neither the name of the author nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific
  prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include _MOL_INCLUDE_
#include "minimize.h"
#include "lbfgs.h"


//Powell/dirPowell is below


void bracket(double* orig, double* dir, double step,
             int ndim, void* prms,
             void (*egfun)(int , double* , void* , double*, double* ),
             double *fb,  double *la, double *lb, double *lc)
{
         int mindim=ndim;
         double* new=_mol_malloc(mindim*sizeof(double) );

         double fa,fc,fnew,lnew;
         double G = 1.618034;
         double TooSmall=1E-10;

         int i, j;
        
         double parabp1, parabp2;

         *la=0;
         *lb=step;
         egfun(ndim, orig, prms, &fa, NULL);
         for (i =0; i<mindim ; i++)
                     new[i] = orig[i]+step*dir[i];
         egfun(ndim, new, prms, fb, NULL);
         if(fa<*fb)
         {
            double temp;
            temp = fa;
            fa = *fb;
            *fb = temp;
            temp =*la;
            *la=*lb;
            *lb=temp;
         }

         *lc=*lb+G*(*lb-*la);
//       printf("la %f lb %f lc %f\n", *la, *lb, *lc);
         for  (i =0; i<mindim ; i++)
             new[i] = orig[i]+*lc*dir[i];
            
         egfun(ndim, new, prms, &fc, NULL);

         while( (fc<=*fb))
         {
//       printf("lb=%f", *lb);

             parabp1 = (*lb -*la)*(*fb-fc);
             parabp2 = (*lb -*lc)*(*fb-fa);
             if( fabs(parabp1-parabp2)>TooSmall)
             {
                 lnew=*lb+((*lb-*lc)*parabp2-(*lb-*la)*parabp1)/(2*(parabp1-parabp2));
                 if(((*lb-lnew)*(lnew-*lc))>0)
                 {
                     for(j =0; j<mindim ; j++)
                        new[j] = orig[j] + lnew*(dir[j]);
                     egfun(ndim, new, prms, &fnew, NULL);
                     if(fnew<fc)
                     {
                        *la=*lb;
                        *fb=fnew;
                        *lb=lnew;
                        free(new);
                        return;
                     }
                     else if(fnew>*fb)
                     {
                        *lc=lnew;
                        free(new);
                        return;
                     }
                     lnew=*lc+G*(*lc-*lb);
                     for(j =0; j<mindim ; j++)
                         new[j] = orig[j] + lnew*(dir[j]);
                     egfun(ndim, new, prms, &fnew, NULL);
                 }
                 else
                 {
                     for(j =0; j<mindim ; j++)
                            new[j] = orig[j] + lnew*(dir[j]);
                     egfun(ndim, new, prms, &fnew, NULL);
                     if(fnew<fc)
                     {
                         *lb=*lc;
                         *lc=lnew;
                         lnew=*lc+G*(*lc-*lb);
                         *fb=fc;
                         fc=fnew;
                         for(j =0; j<mindim ; j++)
                                new[j] = orig[j] + lnew*(dir[j]);
                         egfun(ndim, new, prms, &fnew, NULL);
                     } 
                 }
             }
             else
             {
                 lnew=*lc+G*(*lc-*lb);
                 for(j =0; j<mindim ; j++)
                      new[j] = orig[j] + lnew*(dir[j]);
                 egfun(ndim, new, prms, &fnew, NULL);
             }
             *la=*lb;
             *lb=*lc;
             *lc=lnew;
             fa=*fb;
             *fb=fc;
             fc=fnew;
         }

         free(new);
}
void old_brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double*, double* ),
           double tol, int maxtimes,
           double*  min, double* fmim)
{
           int j;
           double* new=_mol_malloc(ndim*sizeof(double) );
           double s=1E-10;
           double CGOLD=0.3819660;
           int iter;
           double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
           double e=0.0; 
           a=(la < lc ? la : lc); 
           b=(la > lc ? la : lc); 
           x=w=v=lb; 
//           for(j=0; j<ndim; j++)
//                 new[j]=orig[j]+x*dir[j]; 
//           egfun(ndim, new, prms, &fx, NULL); 
           fw=fv=fx=fb;
           for (iter=1;iter<=maxtimes;iter++) 
           { 
                xm=0.5*(a+b);
                tol2=2.0*(tol1=tol*fabs(x)+s);
                if (fabs(x-xm) <= (tol2-0.5*(b-a))) 
                {
                    for(j=0; j<ndim; j++)
                        min[j]=orig[j]+x*dir[j];
                    *fmim=fx; 
                    free(new);
                    return ;
                }
                if (fabs(e) > tol1) 
                {
                    r=(x-w)*(fx-fv);
                    q=(x-v)*(fx-fw);
                    p=(x-v)*q-(x-w)*r;
                    q=2.0*(q-r);
                    if (q > 0.0) p = -p;
                    q=fabs(q);
                    etemp=e;
                    e=d;
                    if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                                                    d=CGOLD*(e=(x >= xm ? a-x : b-x));
                    else 
                    {
                          d=p/q; 
                          u=x+d;
                          if (u-a < tol2 || b-u < tol2)
                                  d=((xm-x) >= 0.0 ? fabs(tol1) : -fabs(tol1));
                    }
                } 
                else 
                {
                    d=CGOLD*(e=(x >= xm ? a-x : b-x));
                }
                u=(fabs(d) >= tol1 ? x+d : x+((d) >= 0.0 ? fabs(tol1) : -fabs(tol1)));
                for(j=0; j<ndim; j++)
                    new[j]=orig[j]+u*dir[j];
                egfun(ndim, new, prms, &fu, NULL);
                if (fu <= fx) 
                { 
                      if(u >= x) a=x; 
                      else b=x; 
                      v=w;w=x;x=u;
                      fv=fw;fw=fx;fx=fu;
                } 
                else 
                {
                    if (u < x) a=u; 
                    else b=u;
                    if (fu <= fw || w == x) 
                    {
                        v=w;
                        w=u;
                        fv=fw;
                        fw=fu;
                    } 
                    else if (fu <= fv || v == x || v == w) 
                    {
                       v=u;
                       fv=fu;
                    }
                } 
           } 
           printf("Too many iterations in brent\n");
           for(j=0; j<ndim; j++)
                        min[j]=orig[j]+x*dir[j];
           *fmim=fx;
           free(new); 
}



void brent(double* orig, double* dir, 
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double*, double* ),
           double tol, unsigned int maxtimes,
           double*  min, double* fmim)
     {
        double lbound1,lbound2;
        unsigned int numIt = 0;
        int doGolden,i;
        int numAct=ndim;
        double lcurmin=lb;
        double lprevmin=lb;
        double lprevmin2=lb;
        double fcurmin=fb;
        double fprevmin=fb;
        double fprevmin2=fb;
        double dMoved=0;
        double dMoved2=0;
        double s= 1E-10;
        double lim;
        double par1,par2,lpartial,lnew,fnew;
        double GOLD=0.381966;
        double*  new=_mol_malloc(numAct*sizeof(double) );
        double denom;
        double lmidpoint;
        double tempdMoved;   
        if( la<lc)
        {
            lbound1=la;
            lbound2=lc;
        }
        else
        {
            lbound1=lc;
            lbound2=la;
        }
        lmidpoint = 0.5*(lbound1+lbound2);

        lim=tol*fabs(lcurmin)+s;
        while(fabs(lcurmin -lmidpoint)>(2*lim - 0.5*(lbound2-lbound1)) && numIt<maxtimes)
        { 
           lmidpoint = 0.5*(lbound1+lbound2);
           lim=tol*fabs(lcurmin)+s;
           doGolden=1;
           if(fabs(dMoved2)>lim) 
           {
              par1=(lcurmin-lprevmin)*(fcurmin-fprevmin2);
              par2=(lcurmin-lprevmin2)*(fcurmin-fprevmin);
              lpartial=((lcurmin-lprevmin)*par1-(lcurmin-lprevmin2)*par2);
              denom=2*(par1-par2);
    
              if(denom>0)
                 lpartial*=-1;
              denom=fabs(denom);
              tempdMoved=dMoved2;
              if(fabs(lpartial)<abs(0.5*denom*tempdMoved)
                 && lpartial>denom* (lbound1-lcurmin)
                 && lpartial<denom*(lbound2-lcurmin))
              {
                 dMoved=lpartial/denom;
                 lnew=lcurmin+dMoved;
                 if(lnew-lbound1<2*lim||lbound2-lnew<2*lim)
                 {
                    doGolden=0;
                    if(lmidpoint-lcurmin !=0)
                        dMoved=lim*(lmidpoint-lcurmin)/fabs(lmidpoint-lcurmin);
                    else
                        dMoved=lim;
                 }
              }
           }

           if( doGolden != 0)
           {
              if(lcurmin>=lmidpoint)
              dMoved2=lbound1-lcurmin;
              else
              dMoved2=lbound2-lcurmin;
              dMoved=GOLD*dMoved2;
           }
      
           if (fabs(dMoved)>lim)
               lnew=lcurmin+dMoved;
           else
           {
              if(dMoved != 0)
              lnew=lcurmin+lim*(dMoved/fabs(dMoved));
              else
              lnew = lcurmin;
           }
           for(i=0; i<numAct; i++)
               new[i] = orig[i] + lnew*dir[i];
           egfun(ndim, new, prms, &fnew, NULL);
           if(fnew<fcurmin)
           {
              if (lnew>=lcurmin)
                 lbound1=lcurmin;
              else
                 lbound2=lcurmin;
                 lprevmin2=lprevmin;
                 fprevmin2=fprevmin;
            
                 lprevmin=lcurmin;
                 fprevmin=fcurmin;
              
                 lcurmin=lnew;
                 fcurmin=fnew;
           }
           else
           {
              if(lnew<lcurmin)
                 lbound1=lnew;
              else
                 lbound2=lnew;
              if(fnew<fprevmin || lprevmin == fcurmin)
              {
                 lprevmin2=lprevmin;
                 fprevmin2=fprevmin;
                 lprevmin=lnew;
                 fprevmin=fnew;
              }
              else if(fnew<fprevmin2 || lprevmin2==lcurmin||lprevmin2==lprevmin)
              {
                 lprevmin2=lnew;
                 fprevmin2=fnew;
              }
           }
           dMoved2=dMoved;
           dMoved=lcurmin-lprevmin;
           numIt++;
        }       
        if(numIt==maxtimes)
             printf("MAXIMUM ITERATIONS REACHED. RETURNING CURRENT BEST\n");
        for(i=0; i<numAct; i++)
           min[i]=orig[i]+lcurmin*dir[i];
  
        *fmim=fcurmin;
        free(new);
}

void dirbrent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, int maxtimes,
           double*  min, double* fmim, double*  grad)
{
       double lbound1,lbound2;
       int numIt=0;
       int doBisect, i, canUse1, canUse2;
       int numAct=ndim;
       double db;
       double lcurmin=lb;
       double lprevmin=lb;
       double lprevmin2=lb;
       double fcurmin=fb;
       double fprevmin=fb;
       double fprevmin2=fb;
       double distMoved=0;
       double distMoved2=0;
       double s= 1E-10;
       double lim;
       double sec1=0,sec2=0,lnew,fnew,dnew;
       double*  new=_mol_malloc(numAct*sizeof(double) );
       double lmidpoint;
       double tempdistMoved;
       double dcurmin;
       double dprevmin;
       double dprevmin2;

       if( la<lc)
       {
           lbound1=la;
           lbound2=lc;
       }
       else
       {
           lbound1=lc;
           lbound2=la;
       }
    
        
       db=0;
       for(i=0; i<numAct; i++)
            db+=grad[i]*dir[i];

       dcurmin= db;
       dprevmin= db;
       dprevmin2= db;

       for(numIt=0; numIt<maxtimes; numIt++ )
       {
           lmidpoint = 0.5*(lbound1+lbound2);
           lim=tol*fabs(lcurmin)+s;

           if    (fabs(lcurmin -lmidpoint)<=(2*lim - 0.5*(lbound2-lbound1))) 
           {
                for(i=0; i<numAct; i++)
                       min[i]=orig[i]+lcurmin*dir[i];
                *fmim=fcurmin;
                free(new);
                return;
           }     

           doBisect=1;
           if(fabs(distMoved2)>lim) 
           {
               canUse1=0;
               canUse2=0;
               if (dprevmin != dcurmin)
               {
                   sec1=(lprevmin-lcurmin)*dcurmin/(dcurmin-dprevmin);
                   canUse1=1;
               }
   
               if (dprevmin2 != dcurmin)
               {
                   sec2=(lprevmin2-lcurmin)*dcurmin/(dcurmin-dprevmin2);
                   canUse2=1;
               }

               tempdistMoved=distMoved2;
               distMoved2=distMoved;

               if (canUse1 !=0)
               {
                   if(((lbound1-lcurmin-sec1) * (lcurmin+sec1-lbound2)<=0)||(dcurmin*sec1>0))
                             canUse1=0;
               }
        
               if (canUse2 !=0)
               {    
                   if(((lbound1-lcurmin-sec2) * (lcurmin+sec2-lbound2)<=0)||(dcurmin*sec2>0))
                                                 canUse2=0;
               }   

               if(canUse2 + canUse1 != 0)
               {
                   if(canUse1+canUse2 == 2)
                   {
                       if(fabs(sec1)<fabs(sec2))
                           distMoved=sec1;
                       else
                           distMoved=sec2;
                   }
                   else if(canUse1 != 0)
                       distMoved=sec1;
                   else if(canUse2 != 0)
                       distMoved=sec2;
                   doBisect=0;
                   if(fabs(distMoved)<= fabs(0.5*tempdistMoved))
                   {
                       lnew=lcurmin+distMoved;
                       if(((lnew-lbound1)<2*lim)||((lbound2-lnew)<2*lim))
                       {
                           if(lcurmin<lmidpoint)
                              distMoved=lim;
                           else 
                              distMoved=-lim;
                              doBisect=0;
                       }
                   }
                   else 
                      doBisect=1;                    
               }    
                     else doBisect=1;
           }
           if( doBisect != 0)
           {
              if(dcurmin>0)
                  distMoved2=lbound1-lcurmin;
              else
                  distMoved2=lbound2-lcurmin;
              distMoved=0.5*distMoved2;
           } 
              
           if (fabs(distMoved)>=lim)
           {
              lnew=lcurmin+distMoved;
              for(i=0; i<numAct; i++)
                  new[i] = orig[i] + lnew*dir[i];
              egfun(ndim, new, prms, &fnew, grad);
           }
           else
           {
              if(distMoved >=0)
                   lnew=lcurmin+lim;
              else
                   lnew = lcurmin-lim;
              for(i=0; i<numAct; i++)
                  new[i] = orig[i] + lnew*dir[i];
              egfun(ndim, new, prms, &fnew, grad);
             
              if (fnew>fcurmin)
              {
                 for(i=0; i<numAct; i++)
                     min[i]=orig[i]+lcurmin*dir[i];
                 *fmim=fcurmin;
                 free(new);
                 return;
              }
           }  
           dnew=0;
           for(i=0; i<numAct; i++)
               dnew+=dir[i]*grad[i];      
           if((fnew<=fcurmin))
           {
               if (lnew>=lcurmin)
                   lbound1=lcurmin;
               else
                   lbound2=lcurmin;

               lprevmin2=lprevmin;
               fprevmin2=fprevmin;
               dprevmin2=dprevmin;

               lprevmin=lcurmin;
               fprevmin=fcurmin;
               dprevmin=dcurmin;

               lcurmin=lnew;
               fcurmin=fnew;
               dcurmin=dnew;
           }
           else
           {
              if(lnew<lcurmin)
                    lbound1=lnew;
              else
                    lbound2=lnew;
              if(fnew<=fprevmin || lprevmin ==lcurmin)
              {
                    lprevmin2=lprevmin;
                    fprevmin2=fprevmin;
                    dprevmin2=dprevmin;
                    lprevmin=lnew;
                    fprevmin=fnew;
                    dprevmin=dnew;
              }
              else if(fnew<=fprevmin2 || lprevmin2==lcurmin||lprevmin2==lprevmin)
              {
                    lprevmin2=lnew;
                    fprevmin2=fnew;
                    dprevmin2=dnew;
              }
           }
       }
       if(numIt==maxtimes)
             printf("MAXIMUM ITERATIONS REACHED. RETURNING CURRENT BEST\n");
       for(i=0; i<numAct; i++)
             min[i]=orig[i]+lcurmin*dir[i];
       *fmim=fcurmin;
       free(new);
}


void powell(double* orig, double* directions, unsigned int maxIt, double tol,
            int ndim, void* prms,
            void (*egfun)(int , double* , void* , double* , double*),
            double* min, double* fmim) 
{
    int 
        j, k, //variables that will be used in for loops
        jmostdec=0,//direction of  the biggest  decrease withing a single directions set
        numAct=ndim;
    unsigned int i;
    double mostdec=0;
    double fbrac,fprev, fprev2=0, val;
    double    la=0,lb=0,lc=0;
    double    test;
    double* pcur=_mol_malloc(numAct*sizeof(double));//point of current minimum
    double* prev=_mol_malloc(numAct*sizeof(double));//previous minimum
    double* newdir=_mol_malloc(numAct*sizeof(double));//average direction moved in one iteration

    
    for(j=0; j<numAct;j++)
    {
        prev[j]=orig[j];
        pcur[j]=orig[j];
    }
    egfun(ndim, pcur, prms, &val, NULL);

    for(i=0; i<maxIt; i++)//loop through once for each direction set(break out if close enough)
    {
        fprev=val;
        mostdec=0;
        jmostdec=0;
        for(j=0; j<numAct ;j++)//loop through once for each direction
        {
            for(k=0;k<numAct;k++)
               newdir[k]=directions[numAct*j+k];//get direction
            fprev2=val;
            bracket(pcur, newdir, 1, ndim, prms, egfun, &fbrac,  &la, &lb, &lc );
            brent(pcur,  newdir,  
                  fbrac, la, lb, lc,
                  ndim, prms, egfun,
                  1E-5, 100,
                  pcur, &val);
            if(fabs(fprev2-val)>mostdec)//get direction of greatest decrease and greatest decrease
            {
                mostdec=fabs(fprev2-val);
                jmostdec=j;
            }
        }

        if( 2*fabs(fprev-val)<=tol*(fabs(fprev)+fabs(val)))//if you reach the desired tolerance end powell
        { 
            for(j=0; j<numAct; j++)
                min[j]=pcur[j];
            
            *fmim=val;
            free(pcur);
            free(prev);
            free(newdir);
            return;
        }


        for(j=0; j<numAct; j++)
            newdir[j]=2*pcur[j]-prev[j];//get new direction

        egfun(ndim, newdir, prms, &fprev2, NULL);

        for(j=0; j<numAct; j++)
        {
            newdir[j]=pcur[j]-prev[j];
            prev[j]=pcur[j];
        }
        
        if(fprev2<fprev)
        {
            test=2.*(fprev-2.*val+fprev2)* _mol_sq(fprev-val-mostdec)-mostdec*_mol_sq(fprev-fprev2);
            if(test<0)
            {//if it makes sence to the direction set
                bracket(pcur, newdir, 1, ndim, prms, egfun,
                        &fbrac,  &la, &lb, &lc );//minimize in new direction
                brent(pcur, newdir,
                      fbrac, la, lb, lc,
                      ndim, prms, egfun,
                      1E-5, 100,
                      pcur, &val);
                for(k=0; k<numAct; k++)
                {
                      directions[jmostdec*numAct+k]=directions[(numAct-1)*(numAct)+k];
                      directions[(numAct-1)*(numAct)+k]= newdir[k];
                }
            }
        }
    }         
    printf("warning maximum number of iterations reached. using current best\n");
        
    for(j=0; j<numAct; j++)
                min[j]=pcur[j];

    *fmim=val;
    free(pcur);
    free(prev);
    free(newdir);
}

void dirMin(double* orig,  unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim)
{
     int numAct = ndim;
     int i,j;
     unsigned int it;

     double* dir=malloc(numAct*sizeof(double)); 
     double* interm=malloc(numAct*sizeof(double)); 
     double* temp=malloc(numAct*sizeof(double));
     double* curmin=malloc(numAct*sizeof(double));
     double* grad = malloc(numAct*sizeof(double));
     double* mvec = malloc(numAct*sizeof(double));

     double val,t1,t2,t3,la,lb,lc,fbrac ;
     double s = 1E-5;
     double fprev;
     double dnorm=0, maxnorm=1.5;

     for(i=0; i<numAct; i++)
            curmin[i]=orig[i];
        
     egfun(ndim, orig, prms, &fprev, grad);

     for(i=0; i<numAct; i++)
     {
        interm[i] =-grad[i] ;
        temp[i]=interm[i];
        dir[i]=interm[i];
     }

     for(it=0; it<maxIt; it++)
     {
        dnorm=0;
        for(i=0; i<numAct; i++)
            dnorm+= dir[i]*dir[i];
        dnorm=sqrt(dnorm);
        if(dnorm>maxnorm)
            for(i=0; i<numAct; i++)
                dir[i]*=maxnorm/dnorm;
                 
        bracket(curmin,  dir, 1,
                ndim, prms, egfun,
                &fbrac,  &la, &lb, &lc);
        
        dirbrent(curmin,dir,  fbrac, la, lb, lc,
                 ndim, prms, egfun,
                 1E-10, 100,
                 mvec, &val, grad);

        if(2*fabs(val-fprev)<= tol*(fabs(fprev)+fabs(val)+s))
        {
           for(j=0; j<numAct; j++)
               min[j]=mvec[j];
           *fmim=val;
           free(dir);
           free(mvec);
           free(interm);
           free(temp);
           return;
        }
    
        for(j=0; j<numAct; j++)
        {
           curmin[j]=mvec[j];
           dir[j]=grad[j];
        }    
        fprev =val;
        t1=0;
        t2=0;
        for(j=0; j<numAct; j++)
        {
            t1+= _mol_sq(interm[j]);
            t2+= (dir[j]+interm[j])*dir[j] ;
        }
      
        if(t1 ==0)
        {
            for(j=0; j<numAct; j++)
                min[j]=curmin[j];
            *fmim=val;
            free(dir);
            free(interm);
            free(temp);
            free(mvec);
            return;
        }
        
        t3=t2/t1;

        for(j=0; j<numAct; j++)
        {
            interm[j] =-dir[j];
            temp[j]=interm[j]+t3*temp[j];
            dir[j]=temp[j];
        }
     }
     printf("maximum num of iterations reached displaying current min\n");
     for(j=0; j<numAct; j++)
                min[j]=curmin[j];
     *fmim=val;
     free(dir);
     free(interm);
     free(temp);
     free(mvec);
     return;
}

void limin(double* orig, double* dir, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim)
{
        double fbrac, la, lb, lc;
        bracket(orig,  dir, 1,
                ndim, prms, egfun,
                &fbrac,  &la, &lb, &lc);
        brent(orig,  dir,
                  fbrac, la, lb, lc,
                  ndim, prms, egfun,
                  tol, maxIt,
                  min, fmim);
}


//Powell/dirPowell is above
 




int lbfgs_new_debug_mode = 0;  
/*
	the value of the lbfgs_new_debug_mode:
		0 : print nothing
		1 : At each step print the iteration number and the value of the function and gnorm and xnorm and the stepsize and the number of function evaluation
		2 : At each step print the gnorm and xnorm and print the g and x vector the stepsize and the number of function evaluation for the corresponding iteration
*/

int progress(void *instance,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fval,const lbfgsfloatval_t fnorm,const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,int n,int k,int ls)

{

    if ( (lbfgs_new_debug_mode > 0 ) && ( k == 1 ) )
	
	{
		printf ("*************************************************\n");
      		printf ("number of variables = %d\n ", n);
	}
	if ( ( lbfgs_new_debug_mode > 0 ) && ( k != 1 ) )
	
	{
		printf("Iteration %d:\n", k);
	        printf (" f =  %.3f   xnorm = %.3f gnorm =  %.3f stepsize = %.3f nfun = %d \n", fval, fnorm, gnorm,step,ls);
		printf("\n");
	} 
	
    if ( lbfgs_new_debug_mode == 2 ) 
    {
	int i;
	printf("vector x : ");
	for ( i = 0 ; i < n ; i++ )
		printf("%f ",x[i]);
	printf("\n");
	printf("vector g :");
	for ( i = 0 ; i < n ; i++ )
		printf("%f ",g[i]);
	printf("\n");
    }

    return 0;
}



/* Minimizes atomgroup with specified parameter set
 minimization types; 0 = LBFGS, 1- Conjugate gradients, 2 -Powell 3 -New implementation of LBFGS
*/
void minimize_ag(mol_min_method min_type, unsigned int maxIt, double tol,struct atomgrp* ag, void* minprms, void (*egfun)(int , double* , void* , double* , double*)){
     int ndim = ag->nactives*3;
     double* xyz =_mol_malloc(ndim*sizeof(double) ); 
     double* minv =_mol_malloc(ndim*sizeof(double) );
	double fmim;
        int i;
        double* directions;
       	ag2array(xyz, ag);
        /*my_en_grad(ndim, xyz, minprms, &fmim, NULL); 
	 printf("Start E-value=%f\n" , fmim );
	*/
	if (min_type == MOL_CONJUGATE_GRADIENTS) dirMin(xyz,maxIt,tol,ndim,minprms, egfun,minv,&fmim);
	if (min_type == MOL_POWELL) {
        directions =_mol_malloc(ndim*ndim*sizeof(double) );
        for(i=0; i<ndim*ndim;i++) directions[i]=0;
        for(i=0; i<ndim; i++) directions[i*ndim+i]=1;
        powell(xyz,directions,maxIt,tol,ndim,minprms, egfun,minv,&fmim);
	free(directions);
	}

	if (min_type == MOL_LBFGS ){
                lbfgs_parameter_t param={
                6, tol, 0, 1e-5,
                maxIt, LBFGS_LINESEARCH_MORETHUENTE, 40,
                1e-20, 1e20, 1e-4, 0.9, 0.9, 1.0e-16,
                0.0, 0, -1,
                };

		ag2array(minv,ag);


                lbfgs_new(ndim,minv,&fmim,egfun,progress,minprms,&param);
	}
	/*
        my_en_grad(ndim, minv, minprms, &fmim, NULL); 
        printf("fmin=%f\n" , fmim );
        */
	array2ag(minv, ag);
	free(minv);
	free(xyz); 
}
