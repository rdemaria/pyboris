#include <math.h>
#include <stdio.h>

#define mu04pi 1e-7
#define clight 299792458

double v_dot_v(const double *v){
   return v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
};

double a_dot_b(const double *a, const double *b){
  return a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
};

void a_cross_b(const double *a, const double *b, double *c){
    c[1] = (a[2]*b[3] - a[3]*b[2]);
    c[2] = (a[3]*b[1] - a[1]*b[3]);
    c[3] = (a[1]*b[2] - a[2]*b[1]);
};

/*
 * x={t,x,y,z};
 * u*m={E,px,py,pz}
 *
 * x[k+1/2] = x[k] + v[k] dt/2
 * v[k+1]   = u'[k+1/2] +q' E[k+1/2]
 * x[k+1]   = x[k+1/2]+ v[k+1] dt/2
 *
 * efield[3] [V/m] 
 * bfield[3] [T]
 * Delta P = dt q *( E + v x B)
 * Delta gamma beta = dt q/m *( E/c + beta x B)
 *
 * */
void boris_drift(double dt, double x[4], double u[4]){
  double gamma;
  x[0] += dt;
  for (int d=1; d<=3; ++d) x[d] += dt * u[d] / u[0] * clight;
}

void boris_kick(double dt, const double x[4], double u[4],
    double q_over_m, const double efield[3], const double bfield[3]){
  double t[4], s[4];
  double u_nph[4], u_nmh[4];
  double u_plus[4], u_minus[4], u_prime[4];
  double u_prime_cross_s[4], u_minus_cross_t[4];
  double gamma_n, gamma_nph;

  double h = q_over_m * dt;

  //printf("%g %g\n",q_over_m,dt);

  for (int d=1; d<=3; ++d) u_nmh[d] = u[d];
  for (int d=1; d<=3; ++d) u_minus[d] = u_nmh[d] + efield[d-1]/clight * 0.5*h;
  gamma_n = sqrt( 1.0 + v_dot_v(u_minus) );
  for (int d=1; d<=3; ++d) t[d] = 0.5*h * bfield[d-1] / gamma_n;
  //printf("%g %g %g %g\n",t[0],t[1],t[2],t[3]);
  for (int d=1; d<=3; ++d) s[d] = 2.0*t[d] / ( 1.0 + v_dot_v(t) );
  a_cross_b(u_minus, t, u_minus_cross_t);
  for (int d=1; d<=3; ++d) u_prime[d] = u_minus[d] + u_minus_cross_t[d];
  a_cross_b(u_prime, s, u_prime_cross_s);
  for (int d=1; d<=3; ++d) u_plus[d] = u_minus[d] + u_prime_cross_s[d];
  for (int d=1; d<=3; ++d) u_nph[d] = u_plus[d] + efield[d-1]/clight * 0.5*h ;
  gamma_nph = sqrt( 1.0 + v_dot_v(u_nph) );
  u[0] = gamma_nph;
  for (int d=1; d<=3; ++d) u[d] = u_nph[d];
  //printf("%g %g %g %g\n",u[0],u[1],u[2],u[3]);
}


void bwire(double cur, const double *pwire, const int nwire,
           const double *pfield, double* bfield){
  int ii;
  double xa,ya,za, xb,yb,zb, xf,yf,zf, xm,ym,zm, dx,dy,dz,
         xaf,yaf,zaf, xbf,ybf,zbf, xmf,ymf,zmf,dzn,
         ta,tb,t1,t2,t3,tatb,dms,daf,dbf,abf;
  for (ii=0; ii<(nwire-1)*3; ii+=3){
    xa=pwire[ii+0]; ya=pwire[ii+1]; za=pwire[ii+2];
    xb=pwire[ii+3]; yb=pwire[ii+4]; zb=pwire[ii+5];
    xf=pfield[0];   yf=pfield[1];   zf=pfield[2];
    dx=(xb-xa); dy=(yb-ya); dz=(zb-za);
    dzn=sqrt(dx*dx+dy*dy+dz*dz);
    dx/=(2*dzn); dy/=(2*dzn); dz/=(2*dzn);
    xaf=xa-xf; yaf=ya-yf; zaf=za-zf;
    xbf=xb-xf; ybf=yb-yf; zbf=zb-zf;
    xmf=xm-xf; ymf=ym-yf; zmf=zm-zf;
    t1=ymf*dz-zmf*dy;
    t2=zmf*dx-xmf*dz;
    t3=xmf*dy-ymf*dx;
    ta =xaf*dx+yaf*dy+zaf*dz;
    tb =xbf*dx+ybf*dy+zbf*dz;
    dms=t1*t1+t2*t2+t3*t3;
    daf=sqrt(xaf*xaf+yaf*yaf+zaf*zaf);
    dbf=sqrt(xbf*xbf+ybf*ybf+zbf*zaf);
    abf=mu04pi*cur*(tb/dbf-ta/daf)/dms;
    bfield[0]+=t1*abf;
    bfield[1]+=t2*abf;
    bfield[2]+=t3*abf;
  };
};
