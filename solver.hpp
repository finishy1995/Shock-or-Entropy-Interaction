/**
 *  solver.hpp
 *  the main header file to solve the problem:
 *  Shock/Entropy Wave Interaction
 *
 *  Created by David Wang on 16/11/5.
 *
 */

#ifndef solver_hpp
#define solver_hpp

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#define times(a) (a*a)
#define mem(a) memset(a, 0, sizeof(a))
using namespace std;

const int maxn = 5000;
const int nValue = 201;
const double xLeft = -5;
const double xRight = 5;
const double CFL = 0.6;
const double maxT = 1.8;
const double GAMA = 1.4;
const double rL = 3.857143;
const double rR = 1;
const double uL = 2.629369;
const double uR = 0;
const double pL = 10.33333;
const double pR = 1;
const double urR = 0.2;
const double RK[3][3] = {1, 0.75, 1.0/3.0, 0, 0.25, 2.0/3.0, 1, 0.25, 2.0/3.0};
const string FILEPATH="/Users/cncuser/Desktop/c_code/CFD3/CFD3/";

class solver
{
public:
    
    solver();
    
    void reset();
    
    void solve();
    
    void output(string);
    
    void setKappa(double);
    
    void setLimiter(int);
    
private:
    
    int limiter, n;
    
    double kappa, gama, x1, x2, cfl, t, deltaX, deltaT;
    
    double u[3][maxn+8], u1[3][maxn+8], uLeft[3][maxn+8], uRight[3][maxn+8], sValue[3][maxn+8];
    
    void getDeltaT();
    
    void roeSolve();
    
    void nonMuscl();
    
    void vanLeer();
    
    void vanAlbada();
    
    void minmod();
    
    void superbee();
};

solver::solver()
{
    reset();
    setLimiter(0);
    setKappa(1.0/3.0);
}

void solver::reset()
{
    int i;
    double xpos;
    
    gama = GAMA;
    n = nValue;
    x1 = xLeft;
    x2 = xRight;
    cfl = CFL;
    t = maxT;
    deltaX = (x2-x1)/(n-1);
    mem(u);
    mem(u1);
    mem(uLeft);
    mem(uRight);
    mem(sValue);
    
    for (i=0;i<(n+8);i++)
    {
        xpos = x1+(i-4)*deltaX;
        if (xpos<-4)
        {
            u[0][i] = rL;
            u[1][i] = rL*uL;
            u[2][i] = 0.5*rL*times(uL)+pL/(gama-1.0);
        } else {
            u[0][i] = rR+urR*sin(5.0*xpos);
            u[1][i] = rR*uR;
            u[2][i] = 0.5*rR*times(uR)+pR/(gama-1.0);
        }
    }
}

void solver::solve()
{
    int i, j, k;
    double tNow = 0;
    
    reset();
    for (i=0; i<(n+8); i++)
        for (j=0; j<3; j++)
            u1[j][i] = u[j][i];
    
    while (tNow<t)
    {
        getDeltaT();
        if (tNow+deltaT>t) deltaT = t-tNow;
        tNow += deltaT;
        
        for (k=0; k<3; k++)
        {
            switch (limiter)
            {
                case 0:
                    nonMuscl();
                    break;
                
                case 1:
                    vanLeer();
                    break;
                
                case 2:
                    vanAlbada();
                    break;
                
                case 3:
                    minmod();
                    break;
                
                case 4:
                    superbee();
                    break;
                
                default:
                    cout<<"Please change limiter between 0 and 4.\n0 : No Muscl\n1 : Van Leer Limiter\n2 : Van Albada Limiter\n3 : Minmod Limiter\n4 : Superbee"<<endl;
                    return;
                    break;
            }
            roeSolve();
            sValue[0][0] = 0;
            sValue[1][0] = 0;
            sValue[2][0] = 0;
            sValue[0][n-1] = 0;
            sValue[1][n-1] = 0;
            sValue[2][n-1] = 0;
            for (i=4; i<(n+4); i++)
                for (j=0; j<3; j++)
                    u1[j][i] = RK[0][k]*u[j][i]+RK[1][k]*u1[j][i]+RK[2][k]*deltaT/deltaX*sValue[j][i];
        }
        
        for (i=0; i<(n+8); i++)
            for (j=0; j<3; j++)
                u[j][i] = u1[j][i];
    }
}

void solver::output(string filename)
{
    ofstream out(FILEPATH+filename);
    if (out.is_open())
    {
        out<<x1<<" "<<x2<<" "<<limiter<<" "<<kappa<<setprecision(12)<<fixed<<"\n";
        for(int i=4; i<(n+4); i++)
            out<<u[0][i]<<" ";
        out<<"\n";
        for(int i=4; i<(n+4); i++)
            out<<(u[1][i]/u[0][i])<<" ";
        out<<"\n";
        for(int i=4; i<(n+4); i++)
            out<<u[2][i]<<" ";
        out<<"\n";
        for(int i=4; i<(n+4); i++)
            out<<((u[2][i]-0.5/u[0][i]*times(u[1][i]))*(gama-1))<<" ";
        out<<endl;
        out.close();
    }
    cout<<"Successful"<<endl;
}

void solver::setKappa(double value)
{
    kappa = value;
}

void solver::setLimiter(int value)
{
    limiter = value;
}

void solver::roeSolve()
{
    double l[4], r[4], s[3], lamda[3], sm[3], sFinal[3], f[3][maxn];
    double tempL, tempR;
    double tol = 1e-6;
    
    for (int i=4; i<(n+3); i++)
    {
        l[0] = uLeft[0][i];
        l[1] = uLeft[1][i]/uLeft[0][i];
        l[2] = uLeft[2][i];
        l[3] = (l[2]-0.5*l[0]*times(l[1]))*(gama-1);
        r[0] = uRight[0][i];
        r[1] = uRight[1][i]/uRight[0][i];
        r[2] = uRight[2][i];
        r[3] = (r[2]-0.5*r[0]*times(r[1]))*(gama-1);
        
        s[0] = (sqrt(l[0])*l[1]+sqrt(r[0])*r[1])/(sqrt(l[0])+sqrt(r[0]));;
        tempL = (l[2]+l[3])/l[0];
        tempR = (r[2]+r[3])/r[0];
        s[2] = (sqrt(l[0])*tempL+sqrt(r[0])*tempR)/(sqrt(l[0])+sqrt(r[0]));
        s[1] = sqrt((gama-1)*(s[2]-0.5*times(s[0])));
        
        if (fabs(s[0])>=tol) lamda[0] = fabs(s[0]);
        else lamda[0]=(s[0]*s[0]+tol*tol)/2.0/tol;
        if (fabs(s[0]+s[1])>=tol) lamda[1] = fabs(s[0]+s[1]);
        else lamda[1]=((s[0]+s[1])*(s[0]+s[1])+tol*tol)/2.0/tol;
        if (fabs(s[0]-s[1])>=tol) lamda[2] = fabs(s[0]-s[1]);
        else lamda[2]=((s[0]-s[1])*(s[0]-s[1])+tol*tol)/2.0/tol;
        
        double S[3][3] = {1,1,1,s[0],s[0]+s[1],s[0]-s[1],0.5*s[0]*s[0],s[2]+s[0]*s[1],s[2]-s[0]*s[1]};
        sm[0] = (gama-1)/(s[1]*s[1])*((r[0]-l[0])*(s[2]-s[0]*s[0])+s[0]*(r[1]*r[0]-l[1]*l[0])-(r[2]-l[2]));
        sm[1]=0.5/s[1]*((r[0]-l[0])*(-s[0]+s[1])+(r[1]*r[0]-l[1]*l[0])-s[1]*sm[0]);
        sm[2]=(r[0]-l[0])-sm[0]-sm[1];
        lamda[0] = lamda[0]*sm[0];
        lamda[1] = lamda[1]*sm[1];
        lamda[2] = lamda[2]*sm[2];
        sFinal[0] = S[0][0]*lamda[0]+S[0][1]*lamda[1]+S[0][2]*lamda[2];
        sFinal[1] = S[1][0]*lamda[0]+S[1][1]*lamda[1]+S[1][2]*lamda[2];
        sFinal[2] = S[2][0]*lamda[0]+S[2][1]*lamda[1]+S[2][2]*lamda[2];
        
        f[0][i] = 0.5*(l[0]*l[1]+r[0]*r[1]-sFinal[0]);
        f[1][i] = 0.5*(l[0]*l[1]*l[1]+l[3]+r[0]*r[1]*r[1]+r[3]-sFinal[1]);
        f[2][i] = 0.5*((l[2]+l[3])*l[1]+(r[2]+r[3])*r[1]-sFinal[2]);
    }
    for (int i=5; i<(n+3); i++)
        for (int j=0; j<3; j++)
            sValue[j][i] = f[j][i-1]-f[j][i];
}

void solver::getDeltaT()
{
    double maxValue = 1e-99;
    double temp;
    for(int i=4; i<(n+4); i++)
    {
        temp = sqrt(gama*((gama-1.0)*(u[2][i]-0.5*u[1][i]*u[1][i]/u[0][i]))/u[0][i])+fabs(u[1][i]/u[0][i]);
        if (temp>maxValue) maxValue = temp;
    }
    deltaT = deltaX/maxValue*cfl;
}

void solver::nonMuscl()
{
    for (int i=4; i<(n+3); i++)
    {
        uLeft[0][i] = u1[0][i];
        uLeft[1][i] = u1[1][i];
        uLeft[2][i] = u1[2][i];
        uRight[0][i] = u1[0][i+1];
        uRight[1][i] = u1[1][i+1];
        uRight[2][i] = u1[2][i+1];
    }
}

void solver::vanLeer()
{
    double tol = 1e-16;
    double deltaP[3], deltaM[3];
    double RL[3], RR[3], phiRL[3], phiRLRV[3], phiRR[3], phiRRRV[3];
    int i, j;

    for (i=5; i<(n+2); i++)
    {
        deltaP[0] = u1[0][i+1]-u1[0][i];
        deltaP[1] = u1[1][i+1]-u1[1][i];
        deltaP[2] = u1[2][i+1]-u1[2][i];
        deltaM[0] = u1[0][i]-u1[0][i-1];
        deltaM[1] = u1[1][i]-u1[1][i-1];
        deltaM[2] = u1[2][i]-u1[2][i-1];
        for (j=0; j<3; j++)
        {
            RL[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            phiRL[j] = (RL[j]+fabs(RL[j]))/(1.0+RL[j]);
            phiRLRV[j] = (1.0/RL[j]+fabs(1.0/RL[j]))/(1.0+1.0/RL[j]);
        }
        for (j=0; j<3; j++)
            uLeft[j][i] = u1[j][i]+0.5*deltaM[j]*0.5*((1-kappa)*phiRL[j]+(1+kappa)*RL[j]*phiRLRV[j]);
        
        deltaP[0] = u1[0][i+2]-u1[0][i+1];
        deltaP[1] = u1[1][i+2]-u1[1][i+1];
        deltaP[2] = u1[2][i+2]-u1[2][i+1];
        deltaM[0] = u1[0][i+1]-u1[0][i];
        deltaM[1] = u1[1][i+1]-u1[1][i];
        deltaM[2] = u1[2][i+1]-u1[2][i];
        for (j=0; j<3; j++)
        {
            RR[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            phiRR[j] = (RR[j]+fabs(RR[j]))/(1.0+RR[j]);
            phiRRRV[j] = (1.0/RR[j]+fabs(1.0/RR[j]))/(1.0+1.0/RR[j]);
        }
        for (j=0; j<3; j++)
            uRight[j][i] = u1[j][i+1]-0.5*deltaP[j]*0.5*((1-kappa)*phiRR[j]/RR[j]+(1+kappa)*phiRRRV[j]);
    }
    
    uLeft[0][4] = u1[0][4];
    uLeft[1][4] = u1[1][4];
    uLeft[2][4] = u1[2][4];
    uRight[0][4] = u1[0][5];
    uRight[1][4] = u1[1][5];
    uRight[2][4] = u1[2][5];
    uLeft[0][n+2] = u1[0][n+2];
    uLeft[1][n+2] = u1[1][n+2];
    uLeft[2][n+2] = u1[2][n+2];
    uRight[0][n+2] = u1[0][n+3];
    uRight[1][n+2] = u1[1][n+3];
    uRight[2][n+2] = u1[2][n+3];
}

void solver::vanAlbada()
{
    double tol = 1e-16;
    double deltaP[3], deltaM[3];
    double RL[3], RR[3], phiRL[3], phiRLRV[3], phiRR[3], phiRRRV[3];
    int i, j;
    
    for (i=5; i<(n+2); i++)
    {
        deltaP[0] = u1[0][i+1]-u1[0][i];
        deltaP[1] = u1[1][i+1]-u1[1][i];
        deltaP[2] = u1[2][i+1]-u1[2][i];
        deltaM[0] = u1[0][i]-u1[0][i-1];
        deltaM[1] = u1[1][i]-u1[1][i-1];
        deltaM[2] = u1[2][i]-u1[2][i-1];
        for (j=0; j<3; j++)
        {
            RL[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            phiRL[j] = (RL[j]+RL[j]*RL[j])/(1.0+RL[j]*RL[j]);
            phiRLRV[j] = (1.0/RL[j]+(1.0/RL[j])*(1.0/RL[j]))/(1.0+(1.0/RL[j])*(1.0/RL[j]));
        }
        for (j=0; j<3; j++)
            uLeft[j][i] = u1[j][i]+0.5*deltaM[j]*0.5*((1-kappa)*phiRL[j]+(1+kappa)*RL[j]*phiRLRV[j]);
        
        deltaP[0] = u1[0][i+2]-u1[0][i+1];
        deltaP[1] = u1[1][i+2]-u1[1][i+1];
        deltaP[2] = u1[2][i+2]-u1[2][i+1];
        deltaM[0] = u1[0][i+1]-u1[0][i];
        deltaM[1] = u1[1][i+1]-u1[1][i];
        deltaM[2] = u1[2][i+1]-u1[2][i];
        for (j=0; j<3; j++)
        {
            RR[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            phiRR[j] = (RR[j]+RR[j]*RR[j])/(1.0+RR[j]*RR[j]);
            phiRRRV[j] = (1.0/RR[j]+(1.0/RR[j])*(1.0/RR[j]))/(1.0+(1.0/RR[j])*(1.0/RR[j]));
        }
        for (j=0; j<3; j++)
            uRight[j][i] = u1[j][i+1]-0.5*deltaP[j]*0.5*((1-kappa)*phiRR[j]/RR[j]+(1+kappa)*phiRRRV[j]);
    }
    
    uLeft[0][4] = u1[0][4];
    uLeft[1][4] = u1[1][4];
    uLeft[2][4] = u1[2][4];
    uRight[0][4] = u1[0][5];
    uRight[1][4] = u1[1][5];
    uRight[2][4] = u1[2][5];
    uLeft[0][n+2] = u1[0][n+2];
    uLeft[1][n+2] = u1[1][n+2];
    uLeft[2][n+2] = u1[2][n+2];
    uRight[0][n+2] = u1[0][n+3];
    uRight[1][n+2] = u1[1][n+3];
    uRight[2][n+2] = u1[2][n+3];
}

void solver::minmod()
{
    double tol = 1e-16;
    double deltaP[3], deltaM[3];
    double RL[3], RR[3], phiRL[3], phiRLRV[3], phiRR[3], phiRRRV[3];
    int i, j;
        
    for (i=5; i<(n+2); i++)
    {
        deltaP[0] = u1[0][i+1]-u1[0][i];
        deltaP[1] = u1[1][i+1]-u1[1][i];
        deltaP[2] = u1[2][i+1]-u1[2][i];
        deltaM[0] = u1[0][i]-u1[0][i-1];
        deltaM[1] = u1[1][i]-u1[1][i-1];
        deltaM[2] = u1[2][i]-u1[2][i-1];
        for (j=0; j<3; j++)
        {
            RL[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            if (RL[j]>0)
            {
                phiRL[j] = fmin(RL[j], 1.0);
                phiRLRV[j] = fmin(1.0/RL[j], 1.0);
            } else {
                phiRL[j] = 0;
                phiRLRV[j] = 0;
            }
        }
        for (j=0; j<3; j++)
            uLeft[j][i] = u1[j][i]+0.5*deltaM[j]*0.5*((1-kappa)*phiRL[j]+(1+kappa)*RL[j]*phiRLRV[j]);
        
        deltaP[0] = u1[0][i+2]-u1[0][i+1];
        deltaP[1] = u1[1][i+2]-u1[1][i+1];
        deltaP[2] = u1[2][i+2]-u1[2][i+1];
        deltaM[0] = u1[0][i+1]-u1[0][i];
        deltaM[1] = u1[1][i+1]-u1[1][i];
        deltaM[2] = u1[2][i+1]-u1[2][i];
        for (j=0; j<3; j++)
        {
            RR[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            if (RR[j]>0)
            {
                phiRR[j] = fmin(RR[j], 1.0);
                phiRRRV[j] = fmin(1.0/RR[j], 1.0);
            } else {
                phiRR[j] = 0;
                phiRRRV[j] = 0;
            }
        }
        for (j=0; j<3; j++)
            uRight[j][i] = u1[j][i+1]-0.5*deltaP[j]*0.5*((1-kappa)*phiRR[j]/RR[j]+(1+kappa)*phiRRRV[j]);
    }
    
    uLeft[0][4] = u1[0][4];
    uLeft[1][4] = u1[1][4];
    uLeft[2][4] = u1[2][4];
    uRight[0][4] = u1[0][5];
    uRight[1][4] = u1[1][5];
    uRight[2][4] = u1[2][5];
    uLeft[0][n+2] = u1[0][n+2];
    uLeft[1][n+2] = u1[1][n+2];
    uLeft[2][n+2] = u1[2][n+2];
    uRight[0][n+2] = u1[0][n+3];
    uRight[1][n+2] = u1[1][n+3];
    uRight[2][n+2] = u1[2][n+3];
}

void solver::superbee()
{
    double tol = 1e-16;
    double deltaP[3], deltaM[3];
    double RL[3], RR[3], phiRL[3], phiRLRV[3], phiRR[3], phiRRRV[3];
    int i, j;
    
    for (i=5; i<(n+2); i++)
    {
        deltaP[0] = u1[0][i+1]-u1[0][i];
        deltaP[1] = u1[1][i+1]-u1[1][i];
        deltaP[2] = u1[2][i+1]-u1[2][i];
        deltaM[0] = u1[0][i]-u1[0][i-1];
        deltaM[1] = u1[1][i]-u1[1][i-1];
        deltaM[2] = u1[2][i]-u1[2][i-1];
        for (j=0; j<3; j++)
        {
            RL[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            phiRL[j] = fmax(fmin(2.0*RL[j], 1.0), fmin(RL[j], 2.0));
            phiRLRV[j] = fmax(fmin(2.0/RL[j], 1.0), fmin(1.0/RL[j], 2.0));
            if (phiRL[j]<0) phiRL[j] = 0;
            if (phiRLRV[j]<0) phiRLRV[j] = 0;
        }
        for (j=0; j<3; j++)
            uLeft[j][i] = u1[j][i]+0.5*deltaM[j]*0.5*((1-kappa)*phiRL[j]+(1+kappa)*RL[j]*phiRLRV[j]);
        
        deltaP[0] = u1[0][i+2]-u1[0][i+1];
        deltaP[1] = u1[1][i+2]-u1[1][i+1];
        deltaP[2] = u1[2][i+2]-u1[2][i+1];
        deltaM[0] = u1[0][i+1]-u1[0][i];
        deltaM[1] = u1[1][i+1]-u1[1][i];
        deltaM[2] = u1[2][i+1]-u1[2][i];
        for (j=0; j<3; j++)
        {
            RR[j] = (deltaP[j]+tol)/(deltaM[j]+tol);
            phiRR[j] = fmax(fmin(2.0*RR[j], 1.0), fmin(RR[j], 2.0));
            phiRRRV[j] = fmax(fmin(2.0/RR[j], 1.0), fmin(1.0/RR[j], 2.0));
            if (phiRR[j]<0) phiRR[j] = 0;
            if (phiRRRV[j]<0) phiRRRV[j] = 0;
        }
        for (j=0; j<3; j++)
            uRight[j][i] = u1[j][i+1]-0.5*deltaP[j]*0.5*((1-kappa)*phiRR[j]/RR[j]+(1+kappa)*phiRRRV[j]);
    }
    
    uLeft[0][4] = u1[0][4];
    uLeft[1][4] = u1[1][4];
    uLeft[2][4] = u1[2][4];
    uRight[0][4] = u1[0][5];
    uRight[1][4] = u1[1][5];
    uRight[2][4] = u1[2][5];
    uLeft[0][n+2] = u1[0][n+2];
    uLeft[1][n+2] = u1[1][n+2];
    uLeft[2][n+2] = u1[2][n+2];
    uRight[0][n+2] = u1[0][n+3];
    uRight[1][n+2] = u1[1][n+3];
    uRight[2][n+2] = u1[2][n+3];
}

#endif /* solver_hpp */
