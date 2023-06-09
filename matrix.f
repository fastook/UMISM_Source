      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(NMAX=29999,NZ=9,NZ1=NZ+1)
      EXTERNAL ASUB,ATSUB
      COMMON /MAT/  AA(NMAX,NZ),KA(NMAX,NZ1),N                   
      DIMENSION B(NMAX),X(NMAX)
      READ(1,*) N
      DO I=1,N
        READ(1,*) KA(I,10)
        DO J=1,KA(I,10)
          READ(1,*) I,KA(I,J),AA(I,J)
        ENDDO
        X(I)=0.D0
      ENDDO
      READ(1,*) (B(I),I=1,N)
      eps=1d-6
      CALL SPARSE(B,N,ASUB,ATSUB,X,RSQ,EPS)                       
      print *,'rsq=',rsq
      PRINT *,(X(I),I=1,N)
      END
      
      SUBROUTINE SPARSE(B,N,ASUB,ATSUB,X,RSQ,EPS)                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      EXTERNAL ASUB,ATSUB                                               
      PARAMETER (NMAX=29999)                                  
      DIMENSION B(N),X(N),G(NMAX),H(NMAX),XI(NMAX),XJ(NMAX)             
      CALL PRECOND(B)
      EPS2=N*EPS**2                                                     
      IRST=0                                                            
1     IRST=IRST+1                                                       
      CALL ASUB(X,XI)                                                   
      print 23,'first asub, xi',(xi(i),i=1,n)
23    format(a,/,(5g13.6))
      RP=0.                                                             
      BSQ=0.                                                            
      DO 11 J=1,N                                                       
        BSQ=BSQ+B(J)**2                                                 
        XI(J)=XI(J)-B(J)                                                
        RP=RP+XI(J)**2                                                  
11    CONTINUE                                                          
      print 24,'rp=',rp
      print 24,'bsq=',bsq
      print 23,'after do 11, xi',(xi(i),i=1,n)
      CALL ATSUB(XI,G)                                                  
      print 23,'first atsub, g',(g(i),i=1,n)
      DO 12 J=1,N                                                       
        G(J)=-G(J)                                                      
        H(J)=G(J)                                                       
12    CONTINUE                                                          
      DO 19 ITER=1,10*N                                                 
      print *,'-------------------------'
        print *,'starting ',iter
        CALL ASUB(H,XI)                                                 
        print 23,'second asub, xi',(xi(i),i=1,n)
        ANUM=0.D0
        ADEN=0.D0
        DO 13 J=1,N                                                     
          ANUM=ANUM+G(J)*H(J)                                           
          ADEN=ADEN+XI(J)**2                                            
13      CONTINUE                                                        
        print 24,'anum=',anum
        print 24,'aden=',aden
        IF(ADEN.EQ.0.)PAUSE 'very singular matrix'                      
        ANUM=ANUM/ADEN                                                  
        print 24,'ratio anum=',anum 
        DO 14 J=1,N                                                     
          XI(J)=X(J)                                                    
          X(J)=X(J)+ANUM*H(J)                                           
14      CONTINUE                                                        
      print 23,'after do 14, x',(x(i),i=1,n)
        CALL ASUB(X,XJ)                                                 
        print 23,'third asub, xj',(xj(i),i=1,n)
        RSQ=0.                                                          
        DO 15 J=1,N                                                     
          XJ(J)=XJ(J)-B(J)                                              
          RSQ=RSQ+XJ(J)**2                                              
15      CONTINUE                                                        
      print 23,'after do 15, xj',(xj(i),i=1,n)
      print 24,'rsq=',rsq
24    format(a,g13.6)
      print 24,'rp=',rp
      print 24,'bsq=',bsq
      print 24,'eps2*bsq=',eps2*bsq
        IF(ITER.GT.1) THEN
          IF(RSQ.EQ.RP.OR.RSQ.LE.BSQ*EPS2)THEN                            
            print *,'normal return',iter,rsq
            RETURN                                                        
          endif                                                           
          IF(RSQ.GT.RP)THEN                                               
            DO 16 J=1,N                                                   
              X(J)=XI(J)                                                  
16          CONTINUE                                                      
            IF(IRST.GE.3) THEN                                            
              print *,'normal return with roundoff',iter,rsq
              RETURN                                                      
            endif                                                         
            GO TO 1                                                       
          ENDIF                                                           
        ENDIF
        RP=RSQ                                                          
        CALL ATSUB(XJ,XI)                                               
        print 23,'second atsub, xi',(xi(i),i=1,n)
        GG=0.                                                           
        DGG=0.                                                          
        DO 17 J=1,N                                                     
          GG=GG+G(J)**2                                                 
          DGG=DGG+(XI(J)+G(J))*XI(J)                                    
17      CONTINUE                                                        
        IF(GG.EQ.0.) THEN                                               
          print *,'normal but rare return',iter,rsq
          RETURN                                                        
        endif                                                           
        GAM=DGG/GG                                                      
        DO 18 J=1,N                                                     
          G(J)=-XI(J)                                                   
          H(J)=G(J)+GAM*H(J)                                            
18      CONTINUE                                                        
19    CONTINUE                                                          
      PAUSE 'too many iterations'                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE ASUB(X,V)                                              
      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(NMAX=29999, NZ=9, NZ1=NZ+1)                                   
      DIMENSION X(NMAX),V(NMAX)                                         
      COMMON /MAT/  AA(NMAX,NZ),index(nmax,NZ1),n                   
      DO I=1,N                                                          
        sum=0.D0
        DO J=1,INDEX(I,NZ1)                                              
          jj=index(i,j)                                                 
          jjj=jj+i-1                                              
          sum=sum+AA(I,J)*X(jj)                                        
        enddo                                                           
        v(i)=sum                                                        
      enddo                                                             
      END                                                               
      SUBROUTINE ATSUB(X,V)                                             
      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(NMAX=29999, NZ=9, NZ1=NZ+1)                                   
      DIMENSION X(NMAX),V(NMAX)                                         
      COMMON /MAT/  AA(NMAX,NZ),index(nmax,NZ1),n                   
      DO I=1,N                                                          
        sum=0.d0
        DO J=1,INDEX(I,NZ1)                                              
          jj=index(i,j)                                                 
          jjj=jj+i-1                                              
          sum=sum+AA(I,J)*X(jj)                                        
        enddo                                                           
        v(i)=sum                                                        
      enddo                                                             
      END                                                               
      SUBROUTINE PRECOND(B)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(NMAX=29999,NZ=9,NZ1=NZ+1)
      DIMENSION B(NMAX)
      COMMON /MAT/  AA(NMAX,NZ),KA(NMAX,NZ1),N                   
      SMALL=1D-6
      do i=1,n
        b(i)=b(i)*SMALL
        do j=1,ka(i,10)
          aa(i,j)=aa(i,j)*SMALL
        enddo
      enddo
      END
