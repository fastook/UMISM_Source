      PARAMETER(MAX=29999)
C **************************************************************C
C                                                               C
C   PROGRAM:  MAKE1A                                            C
C                                                               C
C   DATE:  11 23 87                                             C
C   PROGRAMMER:  FASTOOK                                        C
C                                                               C
C   FUNCTION:                                                   C
C            GENERATES A RECTANGULAR DATA SET WITH              C
C               IXX ROWS                                        C
C               IYY COLUMNS                                     C
C            LOWER LEFT CORNER AT XLL,YLL                       C
C            UPPER RIGHT CORNER AT XUR,YUR                      C
C            AND OUTPUTS IT IN OLD DATA FORMAT                  C
C           PROMPTS FOR BED AND SURF AT EACH STEP               C
C **************************************************************C
      DIMENSION XX(MAX),YY(MAX),NLINE(MAX)
      DIMENSION BGEN(300,300)
      CHARACTER HED*80
      NDT=10
      INTER=1
      DT=100.
      RZERO=0.
      RFLUX=1.E6
      RFLUX1=1.
      RFLUX2=1.
      RFLUX3=1.
      RFLUX4=1.
      print *,'input flux across:'
      print 23,'bottom','right','top','left'
23    format(1x,4a15)
      print '(1x,4f15.6)',rflux1,rflux2,rflux3,rflux4
      read *,rflux1,rflux2,rflux3,rflux4
      ADOT=.1
      BED=100.0
      FRACT=0.0
      FLOWA=2.0
      SLDGB=.02
      PSURF=0.0
      HTICE=0.0
      KODE=0
      WRITE(*,*) 'INPUT TITLE'
      READ(*,2000) HED
2000  FORMAT(A80)
      WRITE(*,*) 'INPUT NUMBER OF COLUMNS AND ROWS'
      READ(*,*) IXX,IYY
      CALL BEDGEN(IXX,IYY,BGEN)
      WRITE(*,*) 'INPUT X,Y - LOWER LEFT CORNER'
      READ(*,*) XLL,YLL
      XLL=XLL*1000.
      YLL=YLL*1000.
      WRITE(*,*) 'INPUT X,Y - UPPER RIGHT CORNER'
      READ(*,*) XUR,YUR
      XUR=XUR*1000.
      YUR=YUR*1000.
      XD=(XUR-XLL)/FLOAT(IXX-1)
      YD=(YUR-YLL)/FLOAT(IYY-1)
      DO 100 I=1,IXX
      XX(I)=XLL+(I-1)*XD
100   CONTINUE
      DO 200 I=1,IYY
      YY(I)=YLL+(I-1)*YD
200   CONTINUE
      NUMNP=IXX*IYY
      NUMEL=(IXX-1)*(IYY-1)
      NUMGBC=2*(IXX-1)+2*(IYY-1)
1000  FORMAT (A80,/,7I6,F8.0)
c     WRITE(7,1000) HED,NUMNP,NUMEL,IXX,IYY,NUMGBC,NDT,INTER,DT
      WRITE(7,*) HED
      WRITE(7,*) NUMNP,NUMEL,IXX,IYY,NUMGBC,NDT,INTER,DT
      N=1
      DO 300 J=1,IYY
        DO 300 I=1,IXX
C         BED=(I-1)*200
C         WRITE(*,*) N
C         READ(*,*) BED,PSURF
          BED=BGEN(I,J)
          bed=100.
c          bed=500.-.0005*sqrt(xx(i)**2+yy(j)**2)
          adot=min(0.25,1.-.0000015*sqrt(xx(i)**2+yy(j)**2))
          adot=min(0.25,1.-.0000030*sqrt(xx(i)**2+yy(j)**2))
          dist=sqrt(xx(i)**2+yy(j)**2)
          adot=0.1*(1-dist/1e6)
      print *,adot
          IF(BED.LT.0.) THEN
            HTICE=0.
          PSURF=0.
          ELSE
            HTICE=BED
            PSURF=BED
          ENDIF
c          htice=5000.+bed-.001*sqrt(xx(i)**2+yy(j)**2)
          temp=-50
          htice=bed
C         IF(I.EQ.1) THEN
          IF(I.EQ.1 .OR. J.EQ.1 .OR. I.EQ.IXX .OR. J.EQ.IYY) THEN
c           WRITE(7,1001) N,KODE+1,XX(I),YY(J),HTICE,
            WRITE(7,*) N,KODE+1,XX(I),YY(J),HTICE,
     &                    ADOT,FRACT,PSURF,BED,FLOWA,SLDGB,
     &                    temp,8,1,0.,0.
          ELSE
c           WRITE(7,1001) N,KODE,XX(I),YY(J),HTICE,
            WRITE(7,*) N,KODE,XX(I),YY(J),HTICE,
     &                    ADOT,FRACT,PSURF,BED,FLOWA,SLDGB,
     &                    temp,8,1,0.,0.,1.9e5
          ENDIF
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3)
        N=N+1
300   CONTINUE
      N=1
      DO 400 I=1,IYY-1
        ISTART=1+(I-1)*IXX
        DO 390 J=1,IXX-1
          KX1=ISTART+J-1
          KX2=KX1+1
          KX3=KX2+IXX
          KX4=KX3-1
c         WRITE(7,1002) N,KX1,KX2,KX3,KX4
          WRITE(7,*) N,KX1,KX2,KX3,KX4,0.,0.
 1002 FORMAT(5I6,1PE17.10)
          N=N+1
390     CONTINUE
400   CONTINUE
      DO I=1,IXX-1
c       WRITE(7,1007) I,I+1,RFLUX1,1.
        WRITE(7,*) I,I+1,RFLUX1,1.
 1007 FORMAT(2I6,1PE13.6,i10)
      ENDDO
      DO I=1,IYY-1
c       WRITE(7,1007) I*IXX,(I+1)*IXX,RFLUX2,2.
        WRITE(7,*) I*IXX,(I+1)*IXX,RFLUX2,2.
      ENDDO
      DO  I=1,IXX-1
c       WRITE(7,1007) NUMNP-I+1,NUMNP-I,RFLUX3,3.
        WRITE(7,*) NUMNP-I+1,NUMNP-I,RFLUX3,3.
      ENDDO
      DO I=IYY-1,1,-1
c       WRITE(7,1007) (I+1)*IXX-(IXX-1),I*IXX-(IXX-1),RFLUX4,4.
        WRITE(7,*) (I+1)*IXX-(IXX-1),I*IXX-(IXX-1),RFLUX4,4.
      ENDDO
      NBOUND=2*IXX+2*(IYY-2)+1
      WRITE(9,*) NBOUND
      N=1
      DO 800 I=1,IXX
        NLINE(N)=I
        N=N+1
800   CONTINUE
      DO 810 I=2,IYY
        NLINE(N)=I*IXX
        N=N+1
810   CONTINUE
      DO 820 I=2,IXX
        NLINE(N)=NUMNP-I+1
        N=N+1
820   CONTINUE
      DO 830 I=1,IYY
        NLINE(N)=NUMNP-IXX+1-I*IXX
        N=N+1
830   CONTINUE
      WRITE(9,3000) (NLINE(N),N=1,NBOUND)
3000  FORMAT(16I6)
      DO 900 I=1,IYY
        DO 890 J=1,IXX
          NLINE(J)=J+(I-1)*IXX
890     CONTINUE
        WRITE(11,3000) IXX
        WRITE(11,3000) (NLINE(JJ),JJ=1,IXX)
900   CONTINUE
      STOP
      END
      FUNCTION XRAND(KX)
      IF(KX.GT.0) IX=KX
      IY=65539*IX
      IF(IY.LT.0) IY=IY+2147483647+1
      XRAND=.4656613E-9*FLOAT(IY)
      IX=IY
      RETURN
      END
      SUBROUTINE BEDGEN(IMAX,JMAX,B)
      DIMENSION B(300,300)
      KX=23
C     PRINT *,'INPUT RANDOM SEED'
C     READ(*,*) KX
      BASE=2000.
      HGTH=1000.
      B(1,1)=(XRAND(KX)-.5)*2000.
C     B(1,1)=1000.
      B(1,1)=BASE
C       DO 10 I=2,IMAX
C         B(I,1)=B(I-1,1)+(XRAND(0)-.5)*HGTH
C 10    CONTINUE
C       DO 20 J=2,JMAX
C         B(1,J)=(B(1,J-1)+B(2,J-1))*.5+(XRAND(0)-.5)*HGTH
C         B(JMAX,J)=(B(JMAX,J-1)+B(JMAX-1,J-1))*.5+(XRAND(0)-.5)*HGTH
C         DO 15 I=2,IMAX-1
C          B(I,J)=(B(I-1,J-1)+B(I,J-1)+B(I+1,J-1))/3.+(XRAND(0)-.5)*HGTH
C 15      CONTINUE
C 20    CONTINUE
      DO 10 J=2,JMAX
        B(1,J)=B(1,J-1)+(XRAND(0)-.5)*HGTH
10    CONTINUE
      DO 20 I=2,IMAX
        B(I,1)=(B(I-1,1)+B(I-1,2))*.5+(XRAND(0)-.5)*HGTH
        B(I,IMAX)=(B(I-1,IMAX)+B(I-1,IMAX-1))*.5+(XRAND(0)-.5)*HGTH
        DO 15 J=2,JMAX-1
         B(I,J)=(B(I-1,J-1)+B(I-1,J)+B(I-1,J+1))/3.+(XRAND(0)-.5)*HGTH
15      CONTINUE
20    CONTINUE
      RETURN
      END
