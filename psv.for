C -PSM/FDM2D-
C                    P and SV WAVES
C*********************************************************************** 
C Calculating P and SV wavefields in homogeneous half-space for a      *
C point source by the Fourier Pseudo-spectral Method.                  *
C **********************************************************************
C Last modified: August 6, 1999                                        * 
C Author: Yanbin WANG                                                  *
C         Department of Earth and Planetary Sciences                   *
C         Faculty of Sciences, Kyushu University                       *
C         Hakozaki 6-10-1, Fukuoka, 812-8581, Japan                    *
C Now at: Department of Geophysics, Peking University                  *
C         100871, Beijing, China                                       * 
C Modified to staggered-grid scheme on 16 June 2005.                   *
C Modified to PSM/FDM hybrid method in February 2006                   *
C                        by Xing Wei and Yanbin Wang.                  *
C Modified for Lanzhou basin on 11 January 2011.                       *
C***********************************************************************
C  PARAMETER                                                           *
C***********************************************************************
      PARAMETER( ZERO = 0.0 , TRUE = 1   , FALSE = 0   )                
      PARAMETER( NX =2048 , NY = 1024, NX2 = NX*2, NY2 = NY*2)            
      PARAMETER( DX = 0.0342  , DY = 0.0342   , DT = 1.0E-3     )           
      PARAMETER( NTMAX = 30000, NWRITE =   500            )             
      PARAMETER( AX=DX, AY=DY, AT=0.1/4.0, T0 =AT*2 )
      PARAMETER( NA=0 )
      PARAMETER( I0=292, J0=292)
      PARAMETER( RMXX=1.0, RMXY=0.0, RMYY=-1.0, RMYX=0.0)
      PARAMETER( FXX=0.0, FYY=0.0, FZZ=0.0)
      PARAMETER( DPXX=0.0, DPYY=0.0, DPZZ=0.0)
      PARAMETER( NST =  512, NSSKIP=NX/NST            )
      PARAMETER( NXA =  20  , NYA =  20                )                
      PARAMETER( NSKIP =10 , NTSKP = NTMAX/NSKIP + 1   )  
      parameter(m=8,n=4)               
C***********************************************************************
C  DIMENSION                                                           *
C***********************************************************************
      REAL     SXX  (NX,NY) ,  SXY  (NX,NY)  ,   SYY  (NX,NY)           
      REAL     DEN  (NX2,NY2) ,  RIG  (NX2,NY2)  ,   LAM  (NX2,NY2)           
      REAL     UX   (NX,NY) ,  UY   (NX,NY)                             
      REAL     VX   (NX,NY) ,  VY   (NX,NY) 
      REAL     UP   (NX,NY) ,  US   (NX,NY)                            
      REAL     DXUX (NX,NY) ,  DXUY (NX,NY)                             
      REAL     DYUX (NX,NY) ,  DYUY (NX,NY)                             
      REAL     DXVX (NX,NY) ,  DXVY (NX,NY)                             
      REAL     DYVX (NX,NY) ,  DYVY (NX,NY)                             
      REAL     DXSXX(NX,NY) ,  DXSXY(NX,NY)                             
      REAL     DYSXY(NX,NY) ,  DYSYY(NX,NY)                             
      REAL     GGG  (NX,NY)                  
      REAL     DVP(2,NX), DVS(2,NX), DDEN(2,NX)                           
      REAL     CYWORK(NY)   ,  CXWORK( NX)                            
      REAL     UXALL  ( NST,NTSKP ) , UYALL  ( NST,NTSKP )              
      REAL     GX( NX ), GY( NY )
      REAL     CLOCK  ( 2 )  
      REAL     ND(NX2), Q1D(NX2), Q2D(NX2)                                     
      INTEGER  ISTX( NST )  , ISTY( NST )  , IMAP ( NX,NY )
      real a(m,n),b(m,n),da(m,n),db(m,n)
C********************************************************************** 
C   OUTPUT FILE NAME                                                  * 
C********************************************************************** 
      CHARACTER*40   ONAME, WNAME, RMODEL                                 
C                                                                       
C# output file name                     
      ONAME = 'PSV_lanzhouQ1_10.OUT'                                           
      WNAME = 'PSV_lanzhouQ1_10.WAV'                                            
C**********************************************************************
C   OBSERVATION POSITIONS (STATIONS)                                  *
C**********************************************************************
C# observation points
      DO I = 1, NST                                    
        ISTX ( I ) =   (I-1)*4+1
        ISTY ( I ) =   NA+1                                             
      END DO
C**********************************************************************                                          
C   Velocity Structure                                                *
C**********************************************************************

      do i=1,m
      do j=1,n
      a(i,j)=i+j
      b(i,j)=i*j
      end do
      end do

      CALL   DIFFXSP ( a, da, CXWORK, m, n, 0, 2, DX )                 
      CALL   DIFFXSM ( b, db, CXWORK, m, n, 0, 2, DX )

      do i=1,m
      do j=1,n
      print *,da(i,j),db(i,j)
      end do
      print *,""
      end do

      return

      OPEN(56,FILE='N4096.dat',status='old')
      OPEN(66,FILE='Q14096.dat',status='old')

      DO I=1,NX2
        READ(56,*) ND(I)
        READ(66,*) Q1D(I)
      END DO

      close(56)
      close(66)

      VPB = 1.70
      VSB = 0.85
      ROB = 1.8

      RRIGB = ROB * VSB**2                                              
      RLANB = ROB * VPB**2 - 2.0 * RRIGB
      RDENB = ROB

      DO I=1,NX2
         DO J=1,NY2
              DEPTH=(J-1)*DY/2.0
	      IF(J.LE.NA*2) THEN
	         RIG(I,J)=0.0
	         DEN(I,J)=RDENB
	         LAM(I,J)=0.0
            ELSE IF(J.EQ.NA*2+1) THEN
              	 RIG(I,J)=RRIGB/2.0
	         DEN(I,J)=RDENB/2.0
	         LAM(I,J)=0.0
	    ELSE IF(DEPTH.LE.-Q1D(i)/1000.0) THEN
                 VPB = 1.70
                 VSB = 0.85
                 ROB = 1.8

                 RRIGB = ROB * VSB**2                                              
                 RLANB = ROB * VPB**2 - 2.0 * RRIGB
                 RDENB = ROB
 
		 RIG(I,J)=ROB *VSB**2
               DEN(I,J)=RDENB
               LAM(I,J)=RLANB

	    ELSE IF(DEPTH.LE.-ND(i)/1000.0) THEN
                 VPB = 4.0
                 VSB = 2.1
                 ROB = 2.4

                 RRIGB = ROB * VSB**2                                              
                 RLANB = ROB * VPB**2 - 2.0 * RRIGB
                 RDENB = ROB

		 RIG(I,J)=ROB *VSB**2
               DEN(I,J)=RDENB
               LAM(I,J)=RLANB

	    ELSE IF(DEPTH.LE.15.0) THEN
                 VPB = 5.8
                 VSB = 3.3
                 ROB = 2.7

                 RRIGB = ROB * VSB**2                                              
                 RLANB = ROB * VPB**2 - 2.0 * RRIGB
                 RDENB = ROB

		 RIG(I,J)=ROB *VSB**2
               DEN(I,J)=RDENB
               LAM(I,J)=RLANB

	    ELSE IF(DEPTH.LE.32.0) THEN
                 VPB = 6.4
                 VSB = 3.6
                 ROB = 2.85

                 RRIGB = ROB * VSB**2                                              
                 RLANB = ROB * VPB**2 - 2.0 * RRIGB
                 RDENB = ROB

		 RIG(I,J)=ROB *VSB**2
               DEN(I,J)=RDENB
               LAM(I,J)=RLANB

	    ELSE 
                 VPB = 6.9
                 VSB = 3.9
                 ROB = 3.1

                 RRIGB = ROB * VSB**2                                              
                 RLANB = ROB * VPB**2 - 2.0 * RRIGB
                 RDENB = ROB

		 RIG(I,J)=ROB *VSB**2
               DEN(I,J)=RDENB
               LAM(I,J)=RLANB

	      END IF
         END DO
       END DO
        print *,rig(1986,1807),den(1345,1423),lam(523,645)
       return

      DO 42 I = 1, NX                                                   
      DO 42 J = 1, NY                                                   
         IMAP ( I,J ) = 0                                               
 42   CONTINUE                                                          
                                                                      
C  --  Check Write

      DO 171 I = 1, NST
         IMAP ( ISTX ( I ), ISTY ( I ) ) = 7
 171  CONTINUE

      DO 68 J = 1, NY, 1                                                
         WRITE(6,'(1H ,256I1)' )  (IMAP(I,J), I = 1, NX )            
  68  CONTINUE                   
                                                                       
                                                                        
C********************************************************************** 
C     INITIALIZE                                                      * 
C********************************************************************** 
      KX = NBEGI2( NX )                                                 
      KY = NBEGI2( NY )                                                 
                                                                        
C -- Grid position

      DO 701 I = 1, NX                                                  
  701 GX(I) = DX * I                                                    

      DO 702 I = 1, NY                                                  
  702 GY(I) = DY * I                                                    
                                                                        
      FTMAX = T0 + AT*2                                                 

      CALL   CLEAR ( VX    , NX, NY, ZERO )                             
      CALL   CLEAR ( VY    , NX, NY, ZERO )                             
      CALL   CLEAR ( UX    , NX, NY, ZERO )                             
      CALL   CLEAR ( UY    , NX, NY, ZERO )
      CALL   CLEAR ( SXX    , NX, NY, ZERO )                             
      CALL   CLEAR ( SXY    , NX, NY, ZERO )                             
      CALL   CLEAR ( SYY    , NX, NY, ZERO )
	      
	OPEN( 26, FILE=WNAME, FORM = 'UNFORMATTED', STATUS='UNKNOWN' )    
      OPEN( 16, FILE=ONAME,  FORM = 'UNFORMATTED', STATUS='UNKNOWN' )       
C# snap shots file
      WRITE ( 16 )  NX, NY, NST                                         
      WRITE ( 16 )  GX, GY                                              
      WRITE ( 16 )  DT, NWRITE, NTSKP                                   
      WRITE ( 16 )  AX, AY, AT, T0                                      
      WRITE ( 16 )  I0, J0, ISTX, ISTY                                  
C# wave form file
      WRITE ( 26 )  NSKIP, DT*NSKIP, NST                                
      WRITE ( 26 )  AX,  AY, AT, T0                                     
      WRITE ( 26 )  I0, J0, ISTX, ISTY                                  
C                                                                       
C********************************************************************** 
C       ABSORBING  BOUNDARY CONDITION                                 * 
C********************************************************************** 
      APARA = 0.015                                                     
      DO 209  I = 1,NX                                                  
      DO 209  J = 1,NY                                                  
         IF(  I .LT. NXA ) THEN                                         
           GG = EXP( -( ( APARA * (NXA-I     ) )**2 ) )                 
         ELSE IF ( I .GT. ( NX-NXA+1 ) ) THEN                           
           GG = EXP( -( ( APARA * (I-NX+NXA-1) )**2 ) )
c         ELSE IF ( J. LT. NYA ) THEN
c           GG = EXP( -( ( APARA * (NYA-J     ) )**2 ) )               
         ELSE IF ( J .GT. ( NY-NYA+1 ) ) THEN                           
           GG = EXP( -( ( APARA * (J-NY+NYA-1) )**2 ) )                 
         ELSE                                                           
           GG = 1.0                                                     
         END IF                                                         
         GGG ( I,J ) = GG                                               
 209  CONTINUE                                                          
C                                                                       
C***********************************************************************
C  TIME STEP START                                                     *
C***********************************************************************
C                                                                       
      NTW = 0   
      NTT = 0 
C
                                                    
      DO  800  IT = 1, NTMAX                                            
C                                                                       
          write(6,* ) IT , '/', NTMAX

       WRITE ( 6,'( 3(F6.2,''s ''))') DTIME(CLOCK), CLOCK(1), CLOCK(2)

      NTT = NTT + 1                                                     
      T     =   DT * (IT-1)                                             
      NTW   =   NTW + 1                                                 

      CALL   DIFFXSP ( VX, DXVX, CXWORK, NX, NY, 0, KX, DX )                 
      CALL   DIFFXSM ( VY, DXVY, CXWORK, NX, NY, 0, KX, DX )

      CALL   FINIDYX (VX,DYVX,NX,NY,DX,DY,DT)
      CALL   FINIDYY (VY,DYVY,NX,NY,DX,DY,DT)
c      CALL   DIFFYSP ( VX, DYVX, CYWORK, NX, NY, KY, DY )            
c      CALL   DIFFYSM ( VY, DYVY, CYWORK, NX, NY, KY, DY )            
C                                                                       
                                                                        
      DO 5 I = 1,NX                                                     
      DO 5 J = 1,NY                                                     
         RAM1  = LAM ( I*2,(J-1)*2+1 )                                            
         RIG1  = RIG ( I*2,(J-1)*2+1 )
   	   RIG2  = RIG ((I-1)*2+1,J*2)
         GG    = GGG(I,J)

         SXXT1IJ = (RAM1+2.0*RIG1) * DXVX(I,J) + RAM1 * DYVY(I,J)            
         SYYT1IJ = (RAM1+2.0*RIG1) * DYVY(I,J) + RAM1 * DXVX(I,J)            
         SXYT1IJ = RIG2 * (DXVY(I,J)+DYVX(I,J))

         SXX(I,J) = SXX(I,J)*GG + DT*SXXT1IJ
	   SYY(I,J) = SYY(I,J)*GG + DT*SYYT1IJ
	   SXY(I,J) = SXY(I,J)*GG + DT*SXYT1IJ	   	                           
   5  CONTINUE                                                          

	DO I=1,NX
	   SYY(I,NA+1) = 0.0
c	   SXY(I,NA) = -SXY(I,NA+1)
	END DO
                                                                        
      CALL   DIFFXSM ( SXX, DXSXX, CXWORK, NX, NY, NA, KX, DX )               
      CALL   DIFFXSP ( SXY, DXSXY, CXWORK, NX, NY, NA, KX, DX )               

      CALL   FINIDYY (SXY,DYSXY,NX,NY,DX,DY,DT)
	CALL   FINIDYX (SYY,DYSYY,NX,NY,DX,DY,DT)
c      CALL   DIFFYSM ( SXY, DYSXY, CYWORK, NX, NY, KY, DY )          
c      CALL   DIFFYSP ( SYY, DYSYY, CYWORK, NX, NY, KY, DY )          

      DO 7 I =  1, NX                                                   
      DO 7 J =  1, NY                                                   
         GG   = GGG ( I,J )                                             
         DENVX = DEN ( (I-1)*2+1,(J-1)*2+1 )
	   DENVY = DEN ( I*2, J*2) 
       
         IF ( T .LT. FTMAX ) THEN
            FX1 = RMXX*FXMXX (I,J, I0,J0, DX,DY, AX,AY, T,T0,AT,
     :                        0.0,0.0)
     :          + RMXY*FXMXZ (I,J, I0,J0, DX,DY, AX,AY, T,T0,AT,
     :                        0.0,0.0)
     :          + FXX *FX    (I,J, I0,J0, DX,DY, AX,AY, T,T0,AT )
     :          + DPXX*EXFORCE(I,J, I0,J0, DX,DY, AX,AY, T,T0,AT)
            FY1 = RMYX*FZMXZ (I,J, I0,J0, DX,DY, AX,AY, T,T0,AT,
     :                        -DX/2.0,-DY/2.0)   
     :          + RMYY*FZMZZ (I,J, I0,J0, DX,DY, AX,AY, T,T0,AT,
     :                        -DX/2.0,-DY/2.0)
     :          + FZZ *FZ    (I,J, I0,J0, DX,DY, AX,AY, T,T0,AT )
     :          + DPZZ*EZFORCE(I,J, I0,J0, DX,DY, AX,AY, T,T0,AT)
         ELSE
            FX1 = 0.0
            FY1 = 0.0
         END IF                                                  

         UXT2IJ  = ( DXSXX (I,J) + DYSXY ( I,J ) + FX1 ) / DENVX        
         UYT2IJ  = ( DXSXY (I,J) + DYSYY ( I,J ) + FY1 ) / DENVY        
                                                                        
         VX ( I,J ) =  VX ( I,J ) * GG + DT * UXT2IJ                
         VY ( I,J ) =  VY ( I,J ) * GG + DT * UYT2IJ                
                                                                        
         UX ( I,J ) =  UX   ( I,J ) * GG + DT * VX( I,J )           
         UY ( I,J ) =  UY   ( I,J ) * GG + DT * VY( I,J )           
                                                                        
  7   CONTINUE                                                          

      IF ( NTT .EQ. NSKIP ) THEN                                        
       DO 25 NS = 1,NST                                                 
          NTT = 0                                                       
          ISX  =  ISTX ( NS )                                           
          ISY  =  ISTY ( NS )                                           
          IT1 = IT/NSKIP + 1                                            
          UXALL ( NS,IT1 ) =  UX   ( ISX , ISY )                        
          UYALL ( NS,IT1 ) =  UY   ( ISX , ISY )                        
 25    CONTINUE                                                         
      END IF                                                 
C                                                                       
C***********************************************************************
C  WRITE  SNAP SHOT                                                    *
C***********************************************************************

      IF( NTW .EQ. NWRITE ) THEN                                        
        NTW = ZERO
	  
        CALL   DIFFXSP ( UX, DXUX, CXWORK, NX, NY, 0, KX, DX )                 
        CALL   DIFFXSM ( UY, DXUY, CXWORK, NX, NY, 0, KX, DX )
        CALL   FINIDYX (UX,DYUX,NX,NY,DX,DY,DT)
        CALL   FINIDYY (UY,DYUY,NX,NY,DX,DY,DT)
c        CALL   DIFFYSP ( UX, DYUX, CYWORK, NX, NY, KY, DY )            
c        CALL   DIFFYSM ( UY, DYUY, CYWORK, NX, NY, KY, DY ) 

        DO  9  I = 1, NX                                                 
        DO  9  J = 1, NY                                                 
           UP ( I,J ) =  DXUX ( I,J ) + DYUY ( I,J )                     
           US ( I,J ) =  DXUY ( I,J ) - DYUX ( I,J )                     
  9     CONTINUE                                                         
                                                      
        WRITE( 16 )  T                                                  
        WRITE( 16 )  UP, US                                            
      END IF                                                            
C                                                                       
C***********************************************************************
C  NEXT TIME STEP                                                      *
C***********************************************************************
C                                                                       
 800  CONTINUE                                                          
C                                                                       
C***********************************************************************
C  WRITE  WAVEFORM                                                      
C***********************************************************************
C                                                                       
      WRITE( 26 ) UXALL, UYALL                                          
C                                                                       
      CLOSE( 16 )                                                       
      CLOSE( 26 )                                                       
C                                                                       
      STOP                                                              
      END                                                               
C                                                                       
C***********************************************************************
C  ***********   END OF MAIN PROGRAM  *****************                *
C***********************************************************************
C                                                                       
      INTEGER  FUNCTION NBEGI2( NX )                                    
C     INPUT :  KX = 2**N                                                
C     OUTPUT:  N                                                        
      NXX = NX                                                          
      KX = 0                                                            
      DO 10 I = 1,NX                                                    
        KX = KX + 1                                                     
        NXX = NXX/2                                                     
        IF( NXX.EQ.1 ) GOTO 99                                          
   10 CONTINUE                                                          
   99 NBEGI2   =    KX                                                  
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE   DIFFX ( A, B, WORK, NX, NY, KX, DX )                 
      REAL         A( NX,NY ) , B( NX,NY ) , WORK ( NX )                
C                                                                       
      DO 100 J = 1, NY                                                  
                                                                        
        DO 1 I = 1,NX                                                   
           WORK( I )  =  A ( I,J )                                      
 1      CONTINUE                                                        
                                                                        
        CALL  FFTR ( WORK, KX, ILL )                                    
                                                                        
        DKX = 6.283185 / NX / DX                                        
        NXD2 = NX / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NXD2+1  ) = 0.0                                          
                                                                        
        DO 2 I = 2, NXD2                                                
          WC = -DKX*( I-1 ) * WORK (      I )                           
          WS =  DKX*( I-1 ) * WORK ( NXD2+I )                           
          WORK (      I ) = WS                                          
          WORK ( NXD2+I ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KX, ILL )                                   
                                                                        
        DO 3 I = 1,NX                                                   
          B ( I,J ) =  WORK ( I )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE   DIFFY ( A, B, WORK, NX, NY, KY, DY )                 
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY )                                           
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( J ) =  A ( I,J )                                       
 1      CONTINUE                                                        
                                                                        
        CALL  FFTR ( WORK, KY, ILL )                                    
                                                                        
        DKY = 6.283185 / NY / DY                                        
        NYD2 = NY / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2                                                
          WC = -DKY*( J-1 ) * WORK (      J )                           
          WS =  DKY*( J-1 ) * WORK ( NYD2+J )                           
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY, ILL )                                   
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( J )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE   DIFFYU ( A, B, WORK, NX, NY, KY, DY, NY2 )           
C     symmetric differentiation
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY2 )                                          
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( NY-J+1 ) =  A ( I,J )                                  
          WORK ( NY+J   ) =  A ( I,J )                                  
 1      CONTINUE                                                        
                                                                        
        CALL  FFTR ( WORK, KY+1, ILL )                                  
                                                                        
        DKY = 6.283185 / NY2 / DY                                       
        NYD2 = NY2 / 2                                                  
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2                                                
          WC = -DKY*( J-1 ) * WORK (      J )                           
          WS =  DKY*( J-1 ) * WORK ( NYD2+J )                           
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY+1, ILL )                                 
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( NY+J )                                    
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE   DIFFYS ( A, B, WORK, NX, NY, KY, DY, NY2 )           
C     free surface differentiation
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY2 )                                          
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK (    J ) =  0.0                                          
          WORK ( NY+J ) =  A ( I,J )                                    
 1      CONTINUE                                                        
                                                                        
        CALL  FFTR ( WORK, KY+1, ILL )                                  
                                                                        
        DKY  = 6.283185 / NY2 / DY                                      
        NYD2 = NY2 / 2                                                  
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2                                                
          WC = -DKY*( J-1 ) * WORK (      J )                           
          WS =  DKY*( J-1 ) * WORK ( NYD2+J )                           
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY+1, ILL )                                 
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) = WORK ( NY+J )                                     
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C
      SUBROUTINE   DIFFXSP ( A, B, WORK, NX, NY, NS, KX, DX)                 
      REAL         A( NX,NY ) , B( NX,NY ) , WORK ( NX )              
C                                                                       
      DO 100 J = 1, NY                                                  

        IF (J.LE.NS) THEN
           DO I=1,NX
              A(I,J) = 0.0
           END DO
        ELSE
                                                                        
        DO 1 I = 1,NX                                                   
           WORK( I )  =  A ( I,J )                                      
 1      CONTINUE                                                        
                                                                        
        CALL  FFTR ( WORK, KX, ILL )                                    
                                                                        
        DKX = 6.283185 / NX / DX                                        
        SFT = 6.283185/NX/2.0
        NXD2 = NX / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NXD2+1  ) = 0.0                                          
                                                                        
        DO 2 I = 2, NXD2
          WC = - DKX*(I-1)* (WORK ( I )*COS(SFT*(I-1))+
     &                       WORK(NXD2+I)*SIN(SFT*(I-1)))                                
          WS =   DKX*(I-1)* (WORK ( NXD2+I )*COS(SFT*(I-1))-
     &                       WORK(I)*SIN(SFT*(I-1)))  

          WORK (      I ) = WS                                          
          WORK ( NXD2+I ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KX, ILL )                                   
                                                                        
        DO 3 I = 1,NX                                                   
          B ( I,J ) =  WORK ( I )                                       
 3      CONTINUE                                                        
                                                                        
        END IF
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                              

      SUBROUTINE   DIFFXSM ( A, B, WORK, NX, NY, NS, KX, DX)                 
      REAL         A( NX,NY ) , B( NX,NY ) , WORK ( NX )           
C                                                                       
      DO 100 J = 1, NY                                                  

        IF ( J.LE.NS) THEN
           DO I=1,NX
              A(I,J) = 0.0
           END DO
        ELSE
        
        DO 1 I = 1,NX                                                   
           WORK( I )  =  A ( I,J )                                      
 1      CONTINUE                                                        
                                                                        
        CALL  FFTR ( WORK, KX, ILL )                                    
                                                                        
        DKX = 6.283185 / NX / DX                                        
        SFT = 6.283185/NX/2.0
        NXD2 = NX / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NXD2+1  ) = 0.0                                          
                                                                        
        DO 2 I = 2, NXD2

          WC = - DKX*(I-1)* (WORK ( I )*COS(SFT*(I-1))-
     &                       WORK(NXD2+I)*SIN(SFT*(I-1)))                                
          WS =   DKX*(I-1)* (WORK ( NXD2+I )*COS(SFT*(I-1))+
     &                       WORK(I)*SIN(SFT*(I-1)))  
          
          WORK (      I ) = WS                                          
          WORK ( NXD2+I ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KX, ILL )                                   
                                                                        
        DO 3 I = 1,NX                                                   
          B ( I,J ) =  WORK ( I )                                       
 3      CONTINUE                                                        
                                                                        
        END IF
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                              
      SUBROUTINE   DIFFYSP ( A, B, WORK, NX, NY, KY, DY )              
        
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY )                                           
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( J ) =  A ( I,J )                                    
 1      CONTINUE

        CALL  FFTR ( WORK, KY, ILL )                                    
                                                                        
        DKY = 6.283185 / NY / DY                                        
        SFT = 6.283185/NY/2.0
        NYD2 = NY / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2
          WC = - DKY*(J-1)* (WORK ( J )*COS(SFT*(J-1))+
     &                       WORK(NYD2+J)*SIN(SFT*(J-1)))                    
                
          WS =   DKY*(J-1)* (WORK ( NYD2+J )*COS(SFT*(J-1))-
     &                       WORK(J)*SIN(SFT*(J-1)))  
                
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY, ILL )                                   
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( J )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE   DIFFYSM ( A, B, WORK, NX, NY, KY, DY )              
        
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY )                                           
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( J ) =  A ( I,J )                                          
 1      CONTINUE  

        CALL  FFTR ( WORK, KY, ILL )                                    
                                                                        
        DKY = 6.283185 / NY / DY                                        
        SFT = 6.283185/NY/2.0
        NYD2 = NY / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2
          WC = - DKY*(J-1)* (WORK ( J )*COS(SFT*(J-1))-
     &                       WORK(NYD2+J)*SIN(SFT*(J-1)))                    
                
          WS =   DKY*(J-1)* (WORK ( NYD2+J )*COS(SFT*(J-1))+
     &                       WORK(J)*SIN(SFT*(J-1)))  
        
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY, ILL )                                   
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( J )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                 
C
      SUBROUTINE   DIFFYUSP ( A, B, WORK, NX, NY, KY, DY, NY2 )                 
C     symmetric differentiation
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY2 )     
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( NY-J+1 ) =  A ( I,J )                                  
          WORK ( NY+J   ) =  A ( I,J )                                  
 1      CONTINUE     
                                                                        
        CALL  FFTR ( WORK, KY+1, ILL )                                    
                                                                        
        DKY = 6.283185 / NY2 / DY                                        
        SFT = 6.283185/NY2/2.0
        NYD2 = NY2 / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2
          WC = - DKY*(J-1)* (WORK ( J )*COS(SFT*(J-1))+
     &                       WORK(NYD2+J)*SIN(SFT*(J-1)))                                
          WS =   DKY*(J-1)* (WORK ( NYD2+J )*COS(SFT*(J-1))-
     &                       WORK(J)*SIN(SFT*(J-1)))  
                
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY+1, ILL )                                   
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( NY+J )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE   DIFFYUSM ( A, B, WORK, NX, NY, KY, DY, NY2 )                 
C     symmetric differentiation
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY2 )     
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( NY-J+1 ) =  A ( I,J )                                  
          WORK ( NY+J   ) =  A ( I,J )                                  
 1      CONTINUE     
                                                                        
        CALL  FFTR ( WORK, KY+1, ILL )                                    
                                                                        
        DKY = 6.283185 / NY2 / DY                                        
        SFT = 6.283185/NY2/2.0
        NYD2 = NY2 / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2
          WC = - DKY*(J-1)* (WORK ( J )*COS(SFT*(J-1))-
     &                       WORK(NYD2+J)*SIN(SFT*(J-1)))                                
          WS =   DKY*(J-1)* (WORK ( NYD2+J )*COS(SFT*(J-1))+
     &                       WORK(J)*SIN(SFT*(J-1)))  
                
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY+1, ILL )                                   
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( NY+J )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                      
C
      SUBROUTINE   DIFFYSSP ( A, B, WORK, NX, NY, KY, DY, NY2 )                 
C     free surface differentiation
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY2 )                                          
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( NY-J+1 ) =  - A ( I,J )                                  
          WORK ( NY+J   ) =  A ( I,J )  
c          WORK (    J ) =  -A ( I, J )                                          
c          WORK ( NY+J ) =  A ( I,J )                                    
 1      CONTINUE             
                                                                        
        CALL  FFTR ( WORK, KY+1, ILL )                                    
                                                                        
        DKY = 6.283185 / NY2 / DY                                        
        SFT = 6.283185/NY2/2.0
        NYD2 = NY2 / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2
          WC = - DKY*(J-1)* (WORK ( J )*COS(SFT*(J-1))+
     &                       WORK(NYD2+J)*SIN(SFT*(J-1)))                                
          WS =   DKY*(J-1)* (WORK ( NYD2+J )*COS(SFT*(J-1))-
     &                       WORK(J)*SIN(SFT*(J-1)))  
                
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC    
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY+1, ILL )                                   
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( NY+J )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                          
C
      SUBROUTINE   DIFFYSSM ( A, B, WORK, NX, NY, KY, DY, NY2 )                 
C     free surface differentiation
      REAL         A( NX,NY ) , B( NX,NY )                              
      REAL         WORK( NY2 )                                          
C                                                                       
      DO 100 I = 1, NX                                                  
                                                                        
        DO 1 J = 1, NY                                                  
          WORK ( NY-J+1 ) =  -A ( I,J )                                  
          WORK ( NY+J   ) =  A ( I,J )  
c          WORK (    J ) =  0.0                                          
c          WORK ( NY+J ) =  A ( I,J )                                    
 1      CONTINUE             
                                                                        
        CALL  FFTR ( WORK, KY+1, ILL )                                    
                                                                        
        DKY = 6.283185 / NY2 / DY                                        
        SFT = 6.283185/NY2/2.0
        NYD2 = NY2 / 2                                                   
                                                                        
        WORK (      1  ) = 0.0                                          
        WORK ( NYD2+1  ) = 0.0                                          
                                                                        
        DO 2 J = 2, NYD2
          WC = - DKY*(J-1)* (WORK ( J )*COS(SFT*(J-1))-
     &                       WORK(NYD2+J)*SIN(SFT*(J-1)))                                
          WS =   DKY*(J-1)* (WORK ( NYD2+J )*COS(SFT*(J-1))+
     &                       WORK(J)*SIN(SFT*(J-1)))  
        
          WORK (      J ) = WS                                          
          WORK ( NYD2+J ) = WC                                          
 2      CONTINUE                                                        
                                                                        
        CALL  FFTRI ( WORK, KY+1, ILL )                                   
                                                                        
        DO 3 J = 1, NY                                                  
          B ( I,J ) =  WORK ( NY+J )                                       
 3      CONTINUE                                                        
                                                                        
 100  CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                              
      SUBROUTINE  CLEAR( A, NX, NY, VALUE )                             
      REAL   A( NX,NY )                                                 
      DO 1 I = 1,NX                                                     
      DO 1 J = 1,NY                                                     
           A( I,J ) = VALUE                                             
 1    CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE PADD ( A, NX, NY, NX0 ,NX1, NY0, NY1, VALUE )          
      REAL A( NX,NY )                                                   
      DO 1 I = NX0 , NX1                                                
      DO 1 J = NY0 , NY1                                                
           A( I,J ) = VALUE                                             
  1   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C
C=====================================================================C 
C     SUBROUTINE BITREV PERFORMS BIT REVERSAL REARRANGEMENT OF        C 
C     REAL ARRAY OF SIZE OF THE FORM 2**L                             C 
C     USAGE: CALL BITREV(A,L,ILL)                                     C 
C          A......1-DIMENSIONAL REAL ARRAY OF SIZE 2**L               C 
C                 INPUT....DATA TO BE REARRANGED                      C 
C                 OUTPUT...REARRANGED DATA                            C 
C          L......INTEGER INPUT...BASE 2 LOGARITHM OF DATA SIZE       C 
C          ILL....INTEGER OUTPUT..ERROR CODE 0: NORMAL TERMINATION    C 
C                                            30000: PARAMETER ERROR   C 
C     AUTHOR: I.NINOMIYA (CHUBU UNIVERSITY,KASUGAI,AICHI JAPAN 487)   C 
C     DATA OF CREATION:  1981-04-01                                   C 
C=====================================================================C 
C---------------------------------------------------------------------C 
      SUBROUTINE BITREV(A,L,ICON)                                       
C---------------------------------------------------------------------C 
      DIMENSION A(*),ITEST(20),INC(20)                                  
C-------- PREPARATION ------------------------------------------------C 
      IF(L.LE.0.OR.L.GT.23) GO TO 8000                                  
      NN=0                                                              
      NR=0                                                              
      I=L-1                                                             
      M=2**I                                                            
      K=2                                                               
      IF(I-2) 60,30,10                                                  
   10 I=I-1                                                             
      ITEST(I-1)=M-K                                                    
      K=K+K                                                             
      INC(I-1)=K-ITEST(I-1)                                             
      IF(I-3) 30,10,10                                                  
C-------- MAIN LOOP --------------------------------------------------C 
   20 NR=INC(I)+NR                                                      
   30 MR=M+NR                                                           
      IF(NR-NN) 50,50,40                                                
   40 W=A(NN+1)                                                         
      A(NN+1)=A(NR+1)                                                   
      A(NR+1)=W                                                         
      MN=M+NN                                                           
      W=A(MN+2)                                                         
      A(MN+2)=A(MR+2)                                                   
      A(MR+2)=W                                                         
   50 NN=NN+2                                                           
      W=A(NN)                                                           
      A(NN)=A(MR+1)                                                     
      A(MR+1)=W                                                         
      NR=K+NR                                                           
   60 MR=M+NR                                                           
      IF(NR-NN) 80,80,70                                                
   70 W=A(NN+1)                                                         
      A(NN+1)=A(NR+1)                                                   
      A(NR+1)=W                                                         
      MN=M+NN                                                           
      W=A(MN+2)                                                         
      A(MN+2)=A(MR+2)                                                   
      A(MR+2)=W                                                         
   80 NN=NN+2                                                           
      W=A(NN)                                                           
      A(NN)=A(MR+1)                                                     
      A(MR+1)=W                                                         
      I=1                                                               
      IF(NN-M) 100,110,110                                              
   90 I=I+1                                                             
  100 IF(NR-ITEST(I)) 20,90,90                                          
C-------- EXIT -------------------------------------------------------C 
  110 ICON=0                                                            
      RETURN                                                            
C-------- ERROR HANDLING ---------------------------------------------C 
 8000 ICON=30000                                                        
      RETURN                                                            
      END                                                               
C=====================================================================C 
C     SUBROUTINE FFTR PERFORMS DISCRETE FOURIER ANALYSIS BY           C 
C     RADIX-2 REAL FFT ALGORITHM                                      C 
C     USAGE: CALL FFTR(A,M,ILL)                                       C 
C          A......1-DIMENSIONAL REAL ARRAY OF SIZE 2**M               C 
C                 INPUT....DATA TO BE ANALYZED                        C 
C                 OUTPUT...ANALYZED DATA                              C 
C                          K-TH COSINE COMPONENT IN A(K+1)            C 
C                          K-TH SINE COMPONENT IN A(N/2+K+1),N=2**M   C 
C          M......INTEGER INPUT...BASE 2 LOGARITHM OF DATA SIZE       C 
C          ILL....INTEGER OUTPUT..ERROR CODE 0: NORMAL TERMINATION    C 
C                                            30000: PARAMETER ERROR   C 
C     SLAVE ROUTINE: BITREV(NUMPAC) FOR BIT REVERSAL REARRANGEMENT    C 
C     AUTHOR: I.NINOMIYA (CHUBU UNIVERSITY,KASUGAI,AICHI JAPAN 487)   C 
C     DATA OF CREATION:  1981-04-01                                   C 
C=====================================================================C 
C---------------------------------------------------------------------C 
      SUBROUTINE FFTR(A,M,ILL)                                          
C---------------------------------------------------------------------C 
      DIMENSION A(*)                                                    
      DOUBLE PRECISION DC(32),DS(32),SS                                 
C---------------------------------------------------------------------C 
C     DC(I)=COS(PI/2**(I+1)),DS(I)=SIN(PI/2**(I+1))                     
C---------------------------------------------------------------------C 
      DATA DC( 1)/ 0.707106781186547531D+00/                            
      DATA DC( 2)/ 0.923879532511286752D+00/                            
      DATA DC( 3)/ 0.980785280403230444D+00/                            
      DATA DC( 4)/ 0.995184726672196887D+00/                            
      DATA DC( 5)/ 0.998795456205172391D+00/                            
      DATA DC( 6)/ 0.999698818696204222D+00/                            
      DATA DC( 7)/ 0.999924701839144545D+00/                            
      DATA DC( 8)/ 0.999981175282601137D+00/                            
      DATA DC( 9)/ 0.999995293809576177D+00/                            
      DATA DC(10)/ 0.999998823451701907D+00/                            
      DATA DC(11)/ 0.999999705862882213D+00/                            
      DATA DC(12)/ 0.999999926465717850D+00/                            
      DATA DC(13)/ 0.999999981616429293D+00/                            
      DATA DC(14)/ 0.999999995404107320D+00/                            
      DATA DC(15)/ 0.999999998851026833D+00/                            
      DATA DC(16)/ 0.999999999712756701D+00/                            
      DATA DC(17)/ 0.999999999928189179D+00/                            
      DATA DC(18)/ 0.999999999982047291D+00/                            
      DATA DC(19)/ 0.999999999995511826D+00/                            
      DATA DC(20)/ 0.999999999998877953D+00/                            
      DATA DC(21)/ 0.999999999999719488D+00/                            
      DATA DC(22)/ 0.999999999999929876D+00/                            
      DATA DC(23)/ 0.999999999999982472D+00/                            
      DATA DC(24)/ 0.999999999999995615D+00/                            
      DATA DC(25)/ 0.999999999999998904D+00/                            
      DATA DC(26)/ 0.999999999999999722D+00/                            
      DATA DC(27)/ 0.999999999999999931D+00/                            
      DATA DC(28)/ 0.999999999999999986D+00/                            
      DATA DC(29)/ 0.100000000000000000D+01/                            
      DATA DC(30)/ 0.100000000000000000D+01/                            
      DATA DC(31)/ 0.100000000000000000D+01/                            
      DATA DC(32)/ 0.100000000000000000D+01/                            
      DATA DS( 1)/ 0.707106781186547531D+00/                            
      DATA DS( 2)/ 0.382683432365089768D+00/                            
      DATA DS( 3)/ 0.195090322016128262D+00/                            
      DATA DS( 4)/ 0.980171403295606036D-01/                            
      DATA DS( 5)/ 0.490676743274180141D-01/                            
      DATA DS( 6)/ 0.245412285229122881D-01/                            
      DATA DS( 7)/ 0.122715382857199263D-01/                            
      DATA DS( 8)/ 0.613588464915447527D-02/                            
      DATA DS( 9)/ 0.306795676296597625D-02/                            
      DATA DS(10)/ 0.153398018628476561D-02/                            
      DATA DS(11)/ 0.766990318742704540D-03/                            
      DATA DS(12)/ 0.383495187571395563D-03/                            
      DATA DS(13)/ 0.191747597310703308D-03/                            
      DATA DS(14)/ 0.958737990959773447D-04/                            
      DATA DS(15)/ 0.479368996030668847D-04/                            
      DATA DS(16)/ 0.239684498084182193D-04/                            
      DATA DS(17)/ 0.119842249050697064D-04/                            
      DATA DS(18)/ 0.599211245264242774D-05/                            
      DATA DS(19)/ 0.299605622633466084D-05/                            
      DATA DS(20)/ 0.149802811316901114D-05/                            
      DATA DS(21)/ 0.749014056584715715D-06/                            
      DATA DS(22)/ 0.374507028292384129D-06/                            
      DATA DS(23)/ 0.187253514146195347D-06/                            
      DATA DS(24)/ 0.936267570730980836D-07/                            
      DATA DS(25)/ 0.468133785365490931D-07/                            
      DATA DS(26)/ 0.234066892682745532D-07/                            
      DATA DS(27)/ 0.117033446341372770D-07/                            
      DATA DS(28)/ 0.585167231706863850D-08/                            
      DATA DS(29)/ 0.292583615853431935D-08/                            
      DATA DS(30)/ 0.146291807926715968D-08/                            
      DATA DS(31)/ 0.731459039633579864D-09/                            
      DATA DS(32)/ 0.365729519816789906D-09/                            
C---------------------------------------------------------------------C 
      N2=1                                                              
      IF(M.LE.0) GOTO 110                                               
      IF(M.EQ.1) GOTO 90                                                
C---------------------------------------------------------------------C 
      CALL BITREV(A,M,ILL)                                              
C---------------------------------------------------------------------C 
      N=2**M                                                            
      F=2.0/FLOAT(N)                                                    
      DO 20 L=1,N,2                                                     
      P=A(L+1)                                                          
      A(L+1)=(A(L)-P)*F                                                 
   20 A(L)=(A(L)+P)*F                                                   
      N1=1                                                              
      N2=2                                                              
      N3=4                                                              
      N4=8                                                              
      DO 70 I=1,M-2                                                     
      DO 30 L=1,N,N3                                                    
      L1=L+N2                                                           
      P=A(L)                                                            
      A(L)=P+A(L1)                                                      
   30 A(L1)=P-A(L1)                                                     
      C=DC(I)                                                           
      S=DS(I)                                                           
      SS=DS(I)+DS(I)                                                    
      CO=1.0                                                            
      SO=0.0                                                            
      DO 50 K=1,N1-1                                                    
      N3K=N3-K                                                          
      N2K=N2-K                                                          
      DO 40 J=1,N,N4                                                    
      L0=J+K                                                            
      L2=L0+N3                                                          
      L1=J+N3K                                                          
      L3=L1+N3                                                          
      P=C*A(L2)-S*A(L3)                                                 
      Q=S*A(L2)+C*A(L3)                                                 
      A(L2)=Q-A(L1)                                                     
      A(L3)=Q+A(L1)                                                     
      A(L1)=A(L0)-P                                                     
      A(L0)=A(L0)+P                                                     
      L1=L0+N2                                                          
      L3=L1+N3                                                          
      L0=J+N2K                                                          
      L2=L0+N3                                                          
      P=A(L2)*S-A(L3)*C                                                 
      Q=A(L2)*C+A(L3)*S                                                 
      A(L2)=Q-A(L1)                                                     
      A(L3)=Q+A(L1)                                                     
      A(L1)=A(L0)-P                                                     
      A(L0)=A(L0)+P                                                     
   40 CONTINUE                                                          
      CN=CO-SS*S                                                        
      SN=SS*C+SO                                                        
      CO=C                                                              
      C=CN                                                              
      SO=S                                                              
   50 S=SN                                                              
      DO 60 J=1,N,N4                                                    
      L0=J+N1                                                           
      L2=L0+N3                                                          
      L1=L2-N2                                                          
      L3=L1+N3                                                          
      P=A(L2)*S-A(L3)*C                                                 
      Q=A(L2)*C+A(L3)*S                                                 
      A(L2)=Q-A(L1)                                                     
      A(L3)=Q+A(L1)                                                     
      A(L1)=A(L0)-P                                                     
   60 A(L0)=A(L0)+P                                                     
      N1=N2                                                             
      N2=N3                                                             
      N3=N4                                                             
   70 N4=N4+N4                                                          
      DO 80 J=2,N1                                                      
      P=A(N2+J)                                                         
      A(N2+J)=A(N+2-J)                                                  
   80 A(N+2-J)=P                                                        
   90 A(1)=(A(1)+A(N2+1))*0.5                                           
      A(N2+1)=A(1)-A(N2+1)                                              
C---------------------------------------------------------------------C 
  100 ILL=0                                                             
      RETURN                                                            
C---------------------------------------------------------------------C 
  110 IF(M.EQ.0) GO TO 100                                              
      ILL=30000                                                         
      RETURN                                                            
      END                                                               
C=====================================================================C 
C     SUBROUTINE FFTRI PERFORMS DISCRETE FOURIER SYMTHESIS BY         C 
C     RADIX-2 REAL FFT ALGORITHM                                      C 
C     USAGE: CALL FFTRI(A,M,ILL)                                      C 
C          A......1-DIMENSIONAL REAL ARRAY OF SIZE 2**M               C 
C                 INPUT....DATA TO BE SYNTHESIZED                     C 
C                          K-TH COSINE COMPONENT IN A(K+1)            C 
C                          K-TH SINE COMPONENT IN A(N/2+K+1),N=2**M   C 
C                 OUTPUT...SYNTHESIZED DATA                           C 
C          M......INTEGER INPUT...BASE 2 LOGARITHM OF DATA SIZE       C 
C          ILL....INTEGER OUTPUT..ERROR CODE 0: NORMAL TERMINATION    C 
C                                            30000: PARAMETER ERROR   C 
C     SLAVE ROUTINE: BITREV(NUMPAC) FOR BIT REVERSAL REARRANGEMENT    C 
C     AUTHOR: I.NINOMIYA (CHUBU UNIVERSITY,KASUGAI,AICHI JAPAN 487)   C 
C     DATA OF CREATION:  1981-04-01                                   C 
C=====================================================================C 
C---------------------------------------------------------------------C 
      SUBROUTINE FFTRI(A,M,ILL)                                         
C---------------------------------------------------------------------C 
      DIMENSION A(*)                                                    
      DOUBLE PRECISION DC(32),DS(32),SS                                 
C---------------------------------------------------------------------C 
C     DC(I)=COS(PI/2**(I+1)),DS(I)=SIN(PI/2**(I+1))                     
C---------------------------------------------------------------------C 
      DATA DC( 1)/ 0.707106781186547531D+00/                            
      DATA DC( 2)/ 0.923879532511286752D+00/                            
      DATA DC( 3)/ 0.980785280403230444D+00/                            
      DATA DC( 4)/ 0.995184726672196887D+00/                            
      DATA DC( 5)/ 0.998795456205172391D+00/                            
      DATA DC( 6)/ 0.999698818696204222D+00/                            
      DATA DC( 7)/ 0.999924701839144545D+00/                            
      DATA DC( 8)/ 0.999981175282601137D+00/                            
      DATA DC( 9)/ 0.999995293809576177D+00/                            
      DATA DC(10)/ 0.999998823451701907D+00/                            
      DATA DC(11)/ 0.999999705862882213D+00/                            
      DATA DC(12)/ 0.999999926465717850D+00/                            
      DATA DC(13)/ 0.999999981616429293D+00/                            
      DATA DC(14)/ 0.999999995404107320D+00/                            
      DATA DC(15)/ 0.999999998851026833D+00/                            
      DATA DC(16)/ 0.999999999712756701D+00/                            
      DATA DC(17)/ 0.999999999928189179D+00/                            
      DATA DC(18)/ 0.999999999982047291D+00/                            
      DATA DC(19)/ 0.999999999995511826D+00/                            
      DATA DC(20)/ 0.999999999998877953D+00/                            
      DATA DC(21)/ 0.999999999999719488D+00/                            
      DATA DC(22)/ 0.999999999999929876D+00/                            
      DATA DC(23)/ 0.999999999999982472D+00/                            
      DATA DC(24)/ 0.999999999999995615D+00/                            
      DATA DC(25)/ 0.999999999999998904D+00/                            
      DATA DC(26)/ 0.999999999999999722D+00/                            
      DATA DC(27)/ 0.999999999999999931D+00/                            
      DATA DC(28)/ 0.999999999999999986D+00/                            
      DATA DC(29)/ 0.100000000000000000D+01/                            
      DATA DC(30)/ 0.100000000000000000D+01/                            
      DATA DC(31)/ 0.100000000000000000D+01/                            
      DATA DC(32)/ 0.100000000000000000D+01/                            
      DATA DS( 1)/ 0.707106781186547531D+00/                            
      DATA DS( 2)/ 0.382683432365089768D+00/                            
      DATA DS( 3)/ 0.195090322016128262D+00/                            
      DATA DS( 4)/ 0.980171403295606036D-01/                            
      DATA DS( 5)/ 0.490676743274180141D-01/                            
      DATA DS( 6)/ 0.245412285229122881D-01/                            
      DATA DS( 7)/ 0.122715382857199263D-01/                            
      DATA DS( 8)/ 0.613588464915447527D-02/                            
      DATA DS( 9)/ 0.306795676296597625D-02/                            
      DATA DS(10)/ 0.153398018628476561D-02/                            
      DATA DS(11)/ 0.766990318742704540D-03/                            
      DATA DS(12)/ 0.383495187571395563D-03/                            
      DATA DS(13)/ 0.191747597310703308D-03/                            
      DATA DS(14)/ 0.958737990959773447D-04/                            
      DATA DS(15)/ 0.479368996030668847D-04/                            
      DATA DS(16)/ 0.239684498084182193D-04/                            
      DATA DS(17)/ 0.119842249050697064D-04/                            
      DATA DS(18)/ 0.599211245264242774D-05/                            
      DATA DS(19)/ 0.299605622633466084D-05/                            
      DATA DS(20)/ 0.149802811316901114D-05/                            
      DATA DS(21)/ 0.749014056584715715D-06/                            
      DATA DS(22)/ 0.374507028292384129D-06/                            
      DATA DS(23)/ 0.187253514146195347D-06/                            
      DATA DS(24)/ 0.936267570730980836D-07/                            
      DATA DS(25)/ 0.468133785365490931D-07/                            
      DATA DS(26)/ 0.234066892682745532D-07/                            
      DATA DS(27)/ 0.117033446341372770D-07/                            
      DATA DS(28)/ 0.585167231706863850D-08/                            
      DATA DS(29)/ 0.292583615853431935D-08/                            
      DATA DS(30)/ 0.146291807926715968D-08/                            
      DATA DS(31)/ 0.731459039633579864D-09/                            
      DATA DS(32)/ 0.365729519816789906D-09/                            
C---------------------------------------------------------------------C 
      IF(M.LE.0) GO TO 100                                              
      N=2**M                                                            
      N4=N                                                              
      N3=N/2                                                            
      N2=N3/2                                                           
      N1=N2/2                                                           
      P=A(1)                                                            
      A(1)=P+A(N3+1)                                                    
      A(N3+1)=P-A(N3+1)                                                 
      IF(M.EQ.1) GOTO 90                                                
      DO 20 J=2,N2                                                      
      P=A(N3+J)                                                         
      A(N3+J)=A(N+2-J)                                                  
   20 A(N+2-J)=P                                                        
      DO 80 I=M-2,0,-1                                                  
      DO 40 L=1,N,N3                                                    
      L1=L+N2                                                           
      P=A(L)                                                            
      A(L)=P+A(L1)                                                      
   40 A(L1)=P-A(L1)                                                     
      IF(N1.EQ.0) GO TO 90                                              
      C=DC(I)                                                           
      S=DS(I)                                                           
      SS=DS(I)+DS(I)                                                    
      CO=1.0                                                            
      SO=0.0                                                            
      DO 60 K=1,N1-1                                                    
      N3K=N3-K                                                          
      N2K=N2-K                                                          
      DO 50 J=1,N,N4                                                    
      L0=J+K                                                            
      L2=L0+N3                                                          
      L1=J+N3K                                                          
      L3=L1+N3                                                          
      P=A(L0)-A(L1)                                                     
      Q=A(L2)+A(L3)                                                     
      A(L0)=A(L0)+A(L1)                                                 
      A(L1)=A(L3)-A(L2)                                                 
      A(L2)=C*P+S*Q                                                     
      A(L3)=C*Q-S*P                                                     
      L1=L0+N2                                                          
      L3=L1+N3                                                          
      L0=J+N2K                                                          
      L2=L0+N3                                                          
      P=A(L0)-A(L1)                                                     
      Q=A(L2)+A(L3)                                                     
      A(L0)=A(L0)+A(L1)                                                 
      A(L1)=A(L3)-A(L2)                                                 
      A(L2)=S*P+C*Q                                                     
      A(L3)=S*Q-C*P                                                     
   50 CONTINUE                                                          
      CN=CO-SS*S                                                        
      SN=SS*C+SO                                                        
      CO=C                                                              
      C=CN                                                              
      SO=S                                                              
   60 S=SN                                                              
      DO 70 J=1,N,N4                                                    
      L0=J+N1                                                           
      L2=L0+N3                                                          
      L1=L2-N2                                                          
      L3=L1+N3                                                          
      P=A(L0)-A(L1)                                                     
      Q=A(L2)+A(L3)                                                     
      A(L0)=A(L0)+A(L1)                                                 
      A(L1)=A(L3)-A(L2)                                                 
      A(L2)=P*S+Q*C                                                     
   70 A(L3)=Q*S-P*C                                                     
      N4=N3                                                             
      N3=N2                                                             
      N2=N1                                                             
   80 N1=N1/2                                                           
C---------------------------------------------------------------------C 
   90 CALL BITREV(A,M,ILL)                                              
C---------------------------------------------------------------------C 
      RETURN                                                            
C---------------------------------------------------------------------C 
  100 IF(M.LT.0) GO TO 110                                              
      ILL=0                                                             
      RETURN                                                            
  110 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C**********************************************************************
C Sunroutine for point source terms
C**********************************************************************
      REAL FUNCTION  FXMXZ(I,J,I0,J0,DX,DZ,AX,AZ,T,T0,AT,XS,ZS)  
      X0 = I0 * DX+XS                                                      
      Z0 = J0 * DZ+ZS                                                      
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX =  HERRMAN ( AX, X, X0 )                                      
      FHZ =-DHERRMAN ( AZ, Z, Z0 )                                      
      FHT =  HERRMAN ( AT, T, T0 )                                      
      FXMXZ  = FHX * FHZ * FHT                                          
      RETURN                                                            
      END                                                               
                                                                        
      REAL FUNCTION  FZMXZ (I,J,I0,J0,DX,DZ,AX,AZ,T,T0,AT,XS,ZS)  
      X0 = I0 * DX+XS                                                    
      Z0 = J0 * DZ+ZS                                                     
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX =-DHERRMAN ( AX, X, X0 )                                      
      FHZ =  HERRMAN ( AZ, Z, Z0 )                                      
      FHT =  HERRMAN ( AT, T, T0 )                                      
      FZMXZ  = FHX * FHZ * FHT                                          
      RETURN                                                            
      END                                                               
                                                                        
      REAL FUNCTION  FZMZZ (I,J,I0,J0,DX,DZ,AX,AZ,T,T0,AT,XS,ZS)  
      X0 = I0 * DX+XS                                                      
      Z0 = J0 * DZ+ZS                                                     
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX =  HERRMAN ( AX, X, X0 )                                      
      FHZ =-DHERRMAN ( AZ, Z, Z0 )                                      
      FHT =  HERRMAN ( AT, T, T0 )                                      
c      FHT =  ricker ( AT, T, T0 )
      FZMZZ  = FHX * FHZ * FHT                                          
      RETURN                                                            
      END                                                               
                                                                        
      REAL FUNCTION  FXMXX(I,J,I0,J0,DX,DZ,AX,AZ,T,T0,AT,XS,ZS)  
      X0 = I0 * DX+XS                                                      
      Z0 = J0 * DZ+ZS                                                      
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX =-DHERRMAN ( AX, X, X0 )                                      
      FHZ =  HERRMAN ( AZ, Z, Z0 )                                      
      FHT =  HERRMAN ( AT, T, T0 )                                      
c      FHT =  ricker ( AT, T, T0 ) 
      FXMXX   = FHX * FHZ * FHT                                         
      RETURN                                                            
      END                                                               
                                                                        
      REAL FUNCTION  FX ( I, J, I0, J0, DX, DZ, AX, AZ, T, T0, AT )     
      X0 = I0 * DX                                                      
      Z0 = J0 * DZ                                                      
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX = HERRMAN ( AX, X, X0 )                                       
      FHZ = HERRMAN ( AZ, Z, Z0 )                                       
      FHT = HERRMAN ( AT, T, T0 )                                       
      FX  = FHX * FHZ * FHT                                             
      RETURN                                                            
      END                                                               

      REAL FUNCTION  FZ ( I, J, I0, J0, DX, DZ, AX, AZ, T, T0, AT )     
      X0 = I0 * DX                                                      
      Z0 = J0 * DZ                                                      
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX = HERRMAN ( AX, X, X0 )                                       
      FHZ = HERRMAN ( AZ, Z, Z0 )                                       
      FHT = HERRMAN ( AT, T, T0 )                                       
      FZ  = FHX * FHZ * FHT                                             
      RETURN                                                            
      END                                                               

      REAL FUNCTION  EXFORCE( I, J, I0, J0, DX, DZ, AX, AZ, T, T0, AT ) 
      X0 = I0 * DX                                                      
      Z0 = J0 * DZ                                                      
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX =-DHERRMAN ( AX, X, X0 )                                      
      FHZ =  HERRMAN ( AZ, Z, Z0 )                                      
      FHT =  HERRMAN ( AT, T, T0 )                                      
      EXFORCE  = FHX * FHZ * FHT                                        
      RETURN                                                            
      END                                                               
                                                                        
      REAL FUNCTION  EZFORCE( I, J, I0, J0, DX, DZ, AX, AZ, T, T0, AT ) 
      X0 = I0 * DX                                                      
      Z0 = J0 * DZ                                                      
      X  = I  * DX                                                      
      Z  = J  * DZ                                                      
      FHX =  HERRMAN ( AX, X, X0 )                                      
      FHZ =-DHERRMAN ( AZ, Z, Z0 )                                      
      FHT =  HERRMAN ( AT, T, T0 )                                      
      EZFORCE  = FHX * FHZ * FHT                                        
      RETURN                                                            
      END                                                               
C                                                                       
      REAL FUNCTION  DHERRMAN ( A, X, X0 )                              
      A2  = 2.0 * A                                                     
      T   = X - X0                                                      
      TD  = ( T + A2 ) / A                                              
      IF( T .LE. -A2 ) THEN                                             
          DHERRMAN  = 0.0                                               
          RETURN                                                        
      ELSE IF( (T .GT. -A2) .AND. (T .LE. -A) ) THEN                    
          DHERRMAN  = (       TD    ) / (A2*A)                       
          RETURN                                                        
      ELSE IF( (T .GT. -A ) .AND. (T .LE.  A) ) THEN                    
          DHERRMAN  = (       -TD    + 2.0     )/(A2*A)     
          RETURN                                                        
      ELSE IF( (T .GT.  A ) .AND. (T .LE. A2) ) THEN                    
          DHERRMAN  = (        TD    - 4.0     )/(A2*A)     
          RETURN                                                        
      ELSE                                                              
          DHERRMAN  = 0.0                                               
          RETURN                                                        
      END IF                                                            
      END
C
      REAL FUNCTION  HERRMAN ( A, X, X0 )                               
      A2  = 2.0 * A                                                     
      T   = X - X0                                                      
      TD  = ( T + A2 ) / A                                              
      IF( T .LE. -A2 ) THEN                                             
          HERRMAN  = 0.0                                                
          RETURN                                                        
      ELSE IF( (T .GT. -A2) .AND. (T .LE. -A) ) THEN                    
          HERRMAN  = ( 0.5 * TD**2 ) / A2                               
          RETURN                                                        
      ELSE IF( (T .GT. -A ) .AND. (T .LE.  A) ) THEN                    
          HERRMAN  = ( -0.5 * TD**2 + 2.0 * TD - 1.0 ) / A2             
          RETURN                                                        
      ELSE IF( (T .GT.  A ) .AND. (T .LE. A2) ) THEN                    
          HERRMAN  = (  0.5 * TD**2 - 4.0 * TD + 8.0 ) / A2             
          RETURN                                                        
      ELSE                                                              
          HERRMAN  = 0.0                                                
          RETURN                                                        
      END IF                                                            
      END                                                               

      real function ricker(xs,xp,x,x0)
      pi=acos(-1.0)
      b=pi*(x-xs-x0)/xp
      if((x.gt.x0).and.(x.lt.(x0+2*xs))) then
          ricker=exp(-b**2)*(b**2-0.5)*sqrt(pi)/2.0
      else
          ricker=0.0
      end if
      return
      end

      real function dricker(xs,xp,x,x0)
      pi=acos(-1.0)
      b=pi*(x-xs-x0)/xp
      if((x.gt.x0).and.(x.lt.(x0+2*xs))) then
          dricker=exp(-b**2)*(1.5-b**2)*b*pi*sqrt(pi)/xp
      else 
          dricker=0.0
      end if
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	SUBROUTINE   FINIDYX (A,DYA,NX,NY,DX,DY,DT) 
c     finite-difference method
      REAL         A(NX,NY), DYA(NX,NY)
	REAL         c0 , c1
	c0 = 9.0/8.0
	c1 = 1.0/24.0
      
      DO i=1,NX
      DYA(i,1) = 1.0/DY * (c0*(A(i, 2) - A(i, 1))- c1*A(i, 3))
      DYA(i,NY-1) = 1.0/DY * (c0*(A(i, NY)-A(i, NY-1))+c1*A(i, NY-2))
      DYA(i,NY) = 1.0/DY * (c0*(-A(i, NY)) + c1*A(i, NY-1))
	END DO

      DO j=2,NY-2
	  DO i=1,NX
           DYA(i,j) = 1.0/DY  
     &                     * (c0*(A(i, j+1) - A(i,   j))
     &                     -  c1*(A(i, j+2) - A(i, j-1)))
    
         END DO
      END DO
	END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	SUBROUTINE   FINIDYY (A,DYA,NX,NY,DX,DY,DT) 
c     finite-difference method
      REAL         A(NX,NY), DYA(NX,NY)
	REAL         c0 , c1
	c0 = 9.0/8.0
	c1 = 1.0/24.0
      
      DO i=1,NX
       DYA(i, 1) = 1.0/DY * (c0* A(i, 1) - c1*A(i, 2))
	DYA(i, 2) = 1.0/DY * (c0*(A(i, 2)-A(i, 1)) - c1*A(i, 3))
	DYA(i,NY) = 1.0/DY * (c0*(A(i, NY)-A(i, NY-1))+c1*A(i, NY-2))
      END DO

      DO j=3,NY-1
	  DO i=1,NX
           DYA(i,j) = 1.0/DY  
     &                     * (c0*(A(i,   j) - A(i, j-1))
     &                     -  c1*(A(i, j+1) - A(i, j-2)))
    
         END DO
      END DO
	END
