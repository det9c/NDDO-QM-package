      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(500),INDXR(500),INDXC(500)
      DO 11 J=1,N
 11                IPIV(J)=0
C *** MAIN LOOP OVER COLS TO BE REDUCED
      DO 22 I=1,N
      BIG=0
C *** OUTER LOOP OF PIVOT ELEMENT SEARCH
      DO 13 J=1,N
        IF(IPIV(J).NE.1)THEN
          DO 12 K=1,N
            IF(IPIV(K).EQ.0)THEN
              IF(ABS(A(J,K)).GE.BIG)THEN
                BIG=ABS(A(J,K))
                IROW=J
                ICOL=K
              ENDIF
            ELSE IF(IPIV(K).GT.1)THEN
c              WRITE(*,109)(IPIV(L),L=1,N),I,J,K
c 109                     FORMAT(' IPIV,I,J,K'/(20I4))
              STOP
            ENDIF
 12                          CONTINUE
        ENDIF
 13              CONTINUE
      IPIV(ICOL)=IPIV(ICOL)+1
      IF(IROW.NE.ICOL)THEN
        DO 14 L=1,N
          DUM=A(IROW,L)
          A(IROW,L)=A(ICOL,L)
          A(ICOL,L)=DUM
 14                    CONTINUE
        DO 15 L=1,M
          DUM=B(IROW,L)
          B(IROW,L)=B(ICOL,L)
          B(ICOL,L)=DUM
 15                    CONTINUE
      ENDIF
      INDXR(I)=IROW
      INDXC(I)=ICOL
      IF(A(ICOL,ICOL).EQ.0)THEN
c        WRITE(*,110)ICOL,ICOL
c 110         FORMAT(' SING MAT A(',I4,',',I4,')=0')
        STOP
      ENDIF
      PIVINV=1/A(ICOL,ICOL)
      A(ICOL,ICOL)=1
      DO 16 L=1,N
        A(ICOL,L)=A(ICOL,L)*PIVINV
 16              CONTINUE
      DO 17 L=1,M
        B(ICOL,L)=B(ICOL,L)*PIVINV
 17              CONTINUE
      DO 21 LL=1,N
        IF(LL.NE.ICOL)THEN
          DUM=A(LL,ICOL)
          A(LL,ICOL)=0
          DO 18 L=1,N
            A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
 18                          CONTINUE
          DO 19 L=1,M
            B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
 19                          CONTINUE
        ENDIF
 21              CONTINUE
 22                           CONTINUE
c      WRITE(6,112)(INDXR(I),I=1,N)
c      WRITE(6,113)(INDXC(I),I=1,N)
 112                               FORMAT(' INDXR',15I4)
 113                                       FORMAT(' INDXC',15I4)
      DO 24 J=1,N
        L=N+1-J
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
 23                          CONTINUE
        ENDIF
 24              CONTINUE
c      WRITE(6,119)A
c 119     FORMAT(' A='/(3E16.6))
      RETURN
      END
