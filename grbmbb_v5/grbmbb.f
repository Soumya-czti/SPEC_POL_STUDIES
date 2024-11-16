       SUBROUTINE grbmbb(EAR,NE,param,ifl,PHOTAR,PHOTER)
       IMPLICIT NONE

       INTEGER NE,IFL
      
       DOUBLE PRECISION EAR(0:NE),PHOTAR(NE),PHOTER(NE)
       DOUBLE PRECISION PARAM(7)
       
       DOUBLE PRECISION y(20),wy(20),ly(6),lwy(6)
       integer i,j,n,ln
       DOUBLE PRECISION aa,bb,nmin,nmax,r0min,ti,tf,r1,r2,T0,c,beta
       DOUBLE PRECISION aval,nu,Grb_flx,Dl,rd,xi1,alp,gam,Angle, Func
       EXTERNAL Grb_flx, Angle, Func

!#####################GAUSS QUADRATURE ABSCISSAS - 6 POINT########################
      DATA ly/-0.93246951420315205,-0.66120938646626459,
     1        -0.23861918608319693,0.23861918608319693,
     2         0.66120938646626459,0.93246951420315205/,
     3     lwy/0.17132449237916600,0.36076157304813861,
     4         0.46791393457269126,0.46791393457269126,
     5         0.36076157304813861,0.17132449237916600/,ln/6/        

!#####################GAUSS QUADRATURE ABSCISSAS - 20 POINT########################

       DATA y/-0.99312859918509488,     -0.96397192727791381,    
     3 -0.91223442825132595,     -0.83911697182221878,    
     5 -0.74633190646015080,     -0.63605368072651502,    
     7 -0.51086700195082713,     -0.37370608871541955,    
     9 -0.22778585114164507,      -7.6526521133497324E-002,
     1   7.6526521133497324E-002, 0.22778585114164507,    
     3  0.37370608871541955,      0.51086700195082713,    
     5  0.63605368072651502,      0.74633190646015080,    
     7  0.83911697182221878,      0.91223442825132595,    
     9  0.96397192727791381,      0.99312859918509488/,

     1   wy /1.7614007139152264E-002,  4.0601429799713346E-002,
     3   6.2672048334049296E-002,  8.3276741576698288E-002,
     5  0.10193011981724021,      0.11819453196151831,    
     7  0.13168863844917650,      0.14209610931838215,    
     9  0.14917298647259827,      0.15275338713072331,    
     1  0.15275338713072331,      0.14917298647259827,    
     3  0.14209610931838215,      0.13168863844917650,    
     5  0.11819453196151831,      0.10193011981724021,    
     7   8.3276741576698288E-002,  6.2672048334049296E-002,
     9   4.0601429799713346E-002,  1.7614007139152264E-002/,
     1      n/20/
!##################################################################################
      COMMON/srcpar/r1,r2,gam,T0,r0min,Planck,alp,rd,pi
      
      DOUBLE PRECISION Echarge,               !Electron charge
     &       Emass,                 !Electron mass
     &       Pmass,                 !Proton mass
     &       Planck,                !Planck constant
     &       ThXn,                  !Thompson Cross-section
     &       LightVel,              !Velocity of Light
     &       Hubble,                !Hubble Constant
     &       omega_m,               !cosmological parameter
     &       omega_v,               !cosmological parameter
     &       Pi,                    !Pi
     &       Boltzman,              !Boltzman Constant
     &       CMBRTemp,              !Cosmic Microwave Background Temperature (in ergs)
     &       Stefan,                !Stefan Boltzman constant
     &       parsec,                !parsec to cm
     &       gravity,               !Gravitiational Constant
     &       m_sun                  !Solar Mass
           
c  All values are in CGS units
    
      Echarge   =  4.8032068E-10 
      Emass     =  9.1093897E-28
      Pmass     =  1836.152*Emass
      Planck    =  6.6260755E-27
      ThXn      =  0.66524616E-24
      LightVel  =  2.99792458E+10
      Hubble    =  71.0
      omega_m   =  0.27
      omega_v   =  1.0-0.27       !Flat universe omega_k = 0
      Pi        =  ASIN(1.d0)*2.d0
      Boltzman  =  1.380658e-16
      CMBRTemp  =  2.728*Boltzman
      Stefan    =  5.67e-5
      parsec    =  3.08e18
      gravity   =  6.67390e-8
      m_sun     =  1.98843e33

      c = LightVel

      gam = param(1)
      alp = param(2) 
      r0min = param(3) !in log
      T0 = param(4) !in keV
      ti = param(5) !in sec
      tf = param(6) !in sec
      Dl = param(7) !in Mpc

      Dl = Dl* 1.0e6 *parsec
      T0 = T0 * 1.0e3 * 1.602e-12 !temp in erg
      r0min = 10**r0min

      beta = DSQRT(1.d0-1.d0/gam**2)
      aval = 2*beta*c*gam**2
      r1 = DLOG(r0min+aval*ti)
      r2 = DLOG(r0min+aval*tf)
      
      xi1 = 2.0*pi*gam**2   !rd/rph(0)
      rd = r0min*xi1
      
      beta = DSQRT(1.d0-1.d0/gam**2)

      DO j=NE,1,-1

         nmax=EAR(J)*1.6021765e-09/planck !keV to Hz
         nmin=EAR(J-1)*1.6021765e-09/planck !keV to Hz
         PHOTAR(J) = 0.0

         aa = DLOG(nmin)
         bb = DLOG(nmax)

         PHOTAR(j) = 0.d0

         DO i = 1,ln
            nu = DEXP(0.5d0*(ly(i)*(bb-aa)+(bb+aa)))

            PHOTAR(j) = PHOTAR(j) + lwy(i)*nu*Grb_flx(nu,y,wy,n) 
   
         ENDDO
         PHOTAR(j) = PHOTAR(j)*(bb-aa)/2.0/(1.0-beta)**2
         PHOTAR(j) = PHOTAR(j)*4.0*Pi*Planck/c**2/DL**2

      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION Grb_flx(nu,y,wy,n)
      IMPLICIT NONE
      DOUBLE PRECISION nu,r1,r2,gam,T0,alp,r0min,Planck,rd,pi
      DOUBLE PRECISION beta,rr,r,val,temp,mu,aa,bb,fint,phi
      INTEGER n, i,k,m
      DOUBLE PRECISION y(20),wy(20),Angle       
      DOUBLE PRECISION facb,facr,facr0,z,za,zb,zint,ztemp
      DOUBLE PRECISION tau,rth,rphth,sw
      COMMON/srcpar/r1,r2,gam,T0,r0min,planck,alp,rd,pi
      EXTERNAL Angle

      Grb_flx = 0.d0
      beta = DSQRT(1.d0-1.d0/gam**2)
      bb = 1.0
      facb = 1.0-beta
      facr0 = r0min/facb
      
      DO k = 1,n

        r =  DEXP(0.5d0*(y(k)*(r2-r1)+(r2+r1)))            
        rr = r/r0min
        facr = r*facb 
        val = Planck*nu*(rr*facb)**alp/gam/T0 

        temp = Angle(r,rd,beta,pi)
        phi = MIN(1.d0/gam,temp)  
 
        aa = COS(phi) 
        fint = 0.0

        DO i = 1,n
           mu =  0.5d0*(y(i)*(bb-aa)+(bb+aa)) 
           rth = facr/(1.0-beta*mu)
            
           IF (ABS(mu) .EQ. 1.0) THEN
               rphth = 1.0
           ELSE    
               rphth = ACOS(mu)/SQRT(1-mu**2)
           ENDIF    
           rphth = facr0*(rphth-beta)

           sw = rth-rphth
           za = MIN(sw,r0min*mu)
           zb = rth
           tau = mu*rphth !tau = rph0*coz(theta)/z
           zint = 0.0

           DO m = 1,n
              
              z =  0.5d0*(y(m)*(zb-za)+(zb+za))
              ztemp = DEXP(-tau/z)
              zint = zint + wy(m)*ztemp

           ENDDO
           temp = DEXP(val/(1.d0-beta*mu)**alp)-1.0
           temp = temp *zint
           temp = mu/(1.-beta*mu)**2/temp

           fint = fint + wy(i)*temp 
        ENDDO
  
        fint = nu**2*fint*(bb-aa)/2.0 
        Grb_flx = Grb_flx + wy(k)*fint*r**3

      ENDDO

      Grb_flx = Grb_flx/2.0 *(r2-r1)/(DEXP(r2)-DEXP(r1)) !required unit #/cm^2/s/Hz

      RETURN
      END

      DOUBLE PRECISION FUNCTION Angle(r,rd,beta,pi)
      IMPLICIT NONE
      DOUBLE PRECISION th1,th2,th,pi,beta,fac,rd,r,Func
      INTEGER i
      EXTERNAL Func

      !f(x) = 1.0-fac*(x/sin(x)-beta)*(1.0-beta*cos(x))
     
      Angle = 0.0
      fac = rd/r/pi/(1.0-beta)
      
      th1 = DLOG(1.0d-10)
      th2 = DLOG(1.57d0)

      IF((Func(fac,beta,DEXP(th1))*Func(fac,beta,DEXP(th2))).GT.0.0)THEN
        RETURN
      ENDIF
      DO i = 1,30
       th = (th1+th2)/2.d0
       IF((Func(fac,beta,DEXP(th))*Func(fac,beta,DEXP(th1))).LT.0.0)THEN
             th2 = th
       ELSE
             th1 = th
       ENDIF
      ENDDO

      IF (Func(fac,beta,DEXP(th)) .LT. 1.0e-04) THEN
          Angle = DEXP(th)
      ELSE 
          print *,"THETA DIDN'T CONVERGE"
          stop
      ENDIF  
      RETURN
      END


      DOUBLE PRECISION FUNCTION Func(fac,beta,x)
      IMPLICIT NONE
      DOUBLE PRECISION fac,x,beta

      Func = 1.0-fac*(x/sin(x)-beta)*(1.0-beta*cos(x))
     
      RETURN
      END
      
   
