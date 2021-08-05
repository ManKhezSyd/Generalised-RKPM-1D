    Module  RKPM1D_Intel
    
    USE MKL95_LAPACK
    USE MKL95_BLAS
    
    !=====================================================================
    Contains
    !=====================================================================
    ! Window Functions including:
    ! 1) WFun1D  : Weight Function in 1D
    ! 2) DWFun1D : First Derivation of Weight Function in 1D
    ! 3) DDWFun1D: Second Derivation of The Weight Function in 1D
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    SUBROUTINE DFAC( N, X)
    !=====================================================================
    
    !==============
    implicit none
    !==============
    
    
    INTEGER     N, J
    
    REAL(8)     X
    
    REAL(8)     FACT
    
    !********************************************************************
    
    FACT = 1._8
    
    DO J = N, 1, -1
        FACT = FACT * J
    END DO
    
    X = FACT
    
    !=====================================================================
    END SUBROUTINE  DFAC
    !=====================================================================
    
    !*********************************************************************
	!=====================================================================
    SUBROUTINE LEGZO(N,X,W)
    !=====================================================================
    !=========================================================
    !Purpose : Compute the zeros of Legendre polynomial Pn(x)
    !in the interval [-1,1], and the corresponding
    !weighting coefficients for Gauss-Legendre integration
    !Input :   n    --- Order of the Legendre polynomial
    !Output:   X(n) --- Zeros of the Legendre polynomial
    !          W(n) --- Corresponding weighting coefficients
    !=========================================================

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION X(N),W(N), X1(N), W1(N)
    N0=(N+1)/2
    DO 45 NR=1,N0
    Z=DCOS(3.141592653589793238_8 * (NR-0.25_8)/(1._8 * N))  
10  Z0=Z
    P=1._8
    DO 15 I=1,NR-1
15      P=P*(Z-X(I))
        F0 = 1._8
        IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0._8
            F1=Z
            DO 20 K=2,N
                PF=(2._8-1._8/K)*Z*F1-(1._8-1._8/K)*F0
                PD=K*(F1-Z*PF)/(1._8-Z*Z)
                F0=F1
20              F1=PF
                IF (Z.EQ.0.0) GO TO 40
                    FD=PF/P
                    Q=0._8
                    DO 35 I=1,NR-1
                        WP=1._8
                        DO 30 J=1,NR-1
                        IF (J.NE.I) WP=WP*(Z-X(J))
30                          CONTINUE
35                          Q=Q+WP
                            GD=(PD-Q*FD)/P
                            Z=Z-FD/GD
                            IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40                              X(NR)=-Z
                                X(N+1-NR)=+Z
                                W(NR)=2._8/((1._8-Z*Z)*PD*PD)
45                              W(N+1-NR)=W(NR)
                                RETURN
                                                              
    !=====================================================================
    END SUBROUTINE LEGZO
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine WFun1D( PrC, KiS, dmax, Ci, WT, Phee)
    !=====================================================================
    ! Variables:
    
    !==============
    implicit none
    !==============
    
    Integer                             WT
    ! Flag Type Integer: This Parameter determines the Type of The Weight
    ! Function

    Real(8)                             LC1, LC2, LC3, LC4
    ! Local Real Careers

    Real(8)                             Phee
    ! Weight Function Value

    Real(8)                             r
    ! Raduis: Distance Between Node and Particle

    Real(8)                             ri
    ! Raduis of Influence Domain of Particle I

    Real(8)                             dmax
    ! Scalling Parameter

    Real(8)                             Ci
    ! Characteristic Nodal Spacing Distance

    Real(8)                             XG1, PrC
    ! Particle Cordinate

    Real(8)                             XG2, KiS
    ! Node Cordinate

    Real(8)                             Z, Z1
    ! Z

    ! End of Variables Declaration
     
    !====================================================================
    ! Executive Part:

    XG1 = PrC
    XG2 = KiS

    r   = XG2 - XG1
    ri  = dmax * Ci

    Z   = DABS( r / ri)
    Z1  = r / ri
     
    !---------------------------------------------------------------------
    Phee = 0._8
    
    IF( Z >= 1._8) THEN
        
        Phee = 0._8
        
    ELSE
        
        IF( WT == 1) THEN

            IF( Z >= 0.5_8) THEN
                Phee = ( 4._8 / 3._8) - ( 4._8 * Z) + ( 4._8 * Z ** 2) - ( 4._8 / 3._8 * Z ** 3)
            ELSE
                Phee = ( 2._8 / 3._8) - ( 4._8 * Z ** 2) + ( 4._8 * Z ** 3)
            END IF

        ELSE IF( WT == 2) THEN

            IF( Z >= 3._8 / 5._8) THEN
                Phee = ( 625._8 / 384._8) - ( 625._8 / 96._8 * Z) + ( 625._8 / 64._8 * Z ** 2) - ( 625._8 / 96._8 * Z ** 3) + ( 625._8 / 384._8 * Z ** 4)
            ELSE IF( Z >= 1._8 / 5._8) THEN
                Phee = ( 55._8 / 96._8) + ( 25._8 / 48._8 * Z) - ( 125._8 / 16._8 * Z ** 2) + ( 625._8 / 48._8 * Z ** 3) - ( 625._8 / 96._8 * Z ** 4)
            ELSE
                Phee = ( 115._8 / 192._8) - ( 125._8 / 32._8 * Z ** 2) + ( 625._8 / 64._8 * Z ** 4)
            END IF

        ELSE IF( WT == 3) THEN
         
            Phee = ( 1._8 - Z1 ** 2) ** 3

        ELSE IF( WT == 4) THEN
        
            Phee = ( 1._8 - Z1 ** 2) ** 9

        ELSE IF( WT == 5) THEN
        
            Phee = ( 1._8 - Z1 ** 2) ** 13
		
        ELSE IF( WT == 6) THEN
        
            Phee = ( 1._8 - Z1 ** 2) ** 16
		
        END IF
        
    END IF

    Phee = 1._8 / ri * Phee

    ! End of Executive Part
    !=====================================================================
    END Subroutine WFun1D
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine DWFun1D( PrC, KiS, dmax, Ci, WT, PhPr)
    !=====================================================================
    ! Variables:
    
    !==============
    implicit none
    !==============
    
    Integer                             WT
    ! Flag Type Integer: This Parameter determines the Type of The Weight
    ! Function

    Real(8)                             LC1, LC2, LC3, LC4
    ! Local Real Careers

    Real(8)                             PhPr
    ! Weight Function Value

    Real(8)                             r
    ! Raduis: Distance Between Node and Particle

    Real(8)                             ri
    ! Raduis of Influence Domain of Particle I

    Real(8)                             dmax
    ! Scalling Parameter

    Real(8)                             Ci
    ! Characteristic Nodal Spacing Distance

    Real(8)                             XG1, PrC
    ! Particle Cordinate

    Real(8)                             XG2, KiS
    ! Node Cordinate

    Real(8)                             Z1, Z
    ! Z

    ! End of Variables Declaration
    !====================================================================
    ! Executive Part:

    XG1 = PrC
    XG2 = KiS

    r   = XG2 - XG1
    ri  = dmax * Ci

    Z1  = ( r / ri)
    Z   = DABS( r / ri)
     
    PhPr = 0._8
	 
    !====================================================================
    IF( Z >= 1._8) THEN
        
        PhPr = 0._8
        
    ELSE
    
        IF( WT == 1) THEN

            IF( Z >= 0.5_8) THEN
                IF( Z1 >= 0._8) THEN
		            PhPr = -4._8 + 8._8 * Z1 - 4._8 * Z1 ** 2
                ELSE
                    PhPr = 4._8 + 8._8 * Z1 + 4._8 * Z1 ** 2
                END IF
            ELSE
                IF( Z1 >= 0._8) THEN
                    PhPr = -8._8 * Z1 + 12._8 * Z1 ** 2
                ELSE
                    PhPr = -8._8 * Z1 - 12._8 * Z1 ** 2
                END IF
            END IF
	   
        ELSE IF( WT == 2) THEN

            IF( Z >= 3._8 / 5._8) THEN
                IF( Z1 >= 0._8) THEN
                    PhPr = -625._8 / 96._8 + 625._8 / 32._8 * Z1 - 625._8 / 32._8 * Z1 ** 2 + 625._8 / 96._8 * Z1 ** 3
                ELSE
                    PhPr = 625._8 / 96._8 + 625._8 / 32._8 * Z1 + 625._8 / 32._8 * Z1 ** 2 + 625._8 / 96._8 * Z1 ** 3
                END IF
            ELSE IF( Z >= 1._8 / 5._8) THEN
                IF( Z1 >= 0._8) THEN
                    PhPr = 25._8 / 48._8 - 125._8 / 8._8 * Z1 + 625._8 / 16._8 * Z1 ** 2 - 625._8 / 24._8 * Z1 ** 3
                ELSE
                    PhPr = -25._8 / 48._8 - 125._8 / 8._8 * Z1 - 625._8 / 16._8 * Z1 ** 2 - 625._8 / 24._8 * Z1 ** 3
                END IF
            ELSE
                IF( Z1 >= 0._8) THEN
                    PhPr = -125._8 / 16._8 * Z1 + 625._8 / 16._8 * Z1 ** 3
                ELSE
                    PhPr = -125._8 / 16._8 * Z1 + 625._8 / 16._8 * Z1 ** 3
                END IF
            END IF

        ELSE IF ( WT == 3) THEN

            PhPr = -6._8 * Z1 * ( Z1 ** 2 - 1._8) ** 2

        ELSE IF( WT == 4) THEN
    
            PhPr = -18._8 * Z1 * ( Z1 ** 2 - 1._8) ** 8

        ELSE IF( WT == 5) THEN

            PhPr = -26._8 * Z1 * ( Z1 ** 2 - 1._8) ** 12

	    ELSE IF( WT == 6) THEN

            PhPr = 32._8 * Z1 * ( Z1 ** 2 - 1._8) ** 15

        END IF
        
    END IF

    PhPr = 1._8 / ri ** 2 * PhPr

    ! End of Executive Part
    !=====================================================================
    END Subroutine DWFun1D
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine DDWFun1D( PrC, KiS, dmax, Ci, WT, PhZe)
    !=====================================================================
    ! Variables:
    
    !==============
    implicit none
    !==============
    
    Integer                             WT
    ! Flag Type Integer: This Parameter determines the Type of The Weight
    ! Function

    Real(8)                             LC1, LC2, LC3, LC4
    ! Local Real Careers

    Real(8)                             PhZe
    ! Weight Function Value

    Real(8)                             r
    ! Raduis: Distance Between Node and Particle

    Real(8)                             ri
    ! Raduis of Influence Domain of Particle I

    Real(8)                             dmax
    ! Scalling Parameter

    Real(8)                             Ci
    ! Characteristic Nodal Spacing Distance

    Real(8)                             XG1, PrC
    ! Particle Cordinate

    Real(8)                             XG2, KiS
    ! Node Cordinate

    Real(8)                             Z1, Z
    ! Z

    ! End of Variables Declaration
    !=====================================================================
    ! Executive Part:

    XG1 = PrC
    XG2 = KiS

    r   = XG2 - XG1
    ri  = dmax * Ci

    Z1  = ( r / ri)
    Z   = DABS( r / ri)
    
    PhZe = 0._8
	 
    !---------------------------------------------------------------------
    
    IF( Z >= 1._8) THEN
        PhZe = 0._8
    ELSE
        IF( WT == 1) THEN

            IF( Z >= 0.5_8) THEN
                PhZe = 8._8 - 8._8 * Z
            ELSE
                PhZe = -8._8 + 24._8 * Z
            END IF
	   
        ELSE IF( WT == 2) THEN
	   
            IF( Z >= 3._8 / 5._8) THEN
                PhZe = 625._8 / 32._8 - 625._8 / 16._8 * Z + 625._8 / 32._8 * Z ** 2
            ELSE IF( Z >= 1._8 / 5._8) THEN
                PhZe = -125._8 / 8._8 + 625._8 / 8._8 * Z - 625._8 / 8._8 * Z ** 2
            ELSE
                PhZe = -125._8 / 16._8 + 1875._8 / 16._8 * Z ** 2
            END IF

        ELSE IF( WT == 3) THEN
            
            PhZe = -6._8 * ( Z1 ** 2 - 1._8) * ( 5._8 * Z1 ** 2 - 1._8)

	    ELSE IF( WT == 4) THEN

            PhZe = -18._8 * ( ( Z1 ** 2 - 1._8) ** 7) * ( 17._8 * Z1 ** 2 - 1._8)
	
        ELSE IF( WT == 5) THEN
            
            PhZe = -26._8 * ( ( Z1 ** 2 - 1._8) ** 11) * ( 25._8 * Z1 ** 2 - 1._8)

	    ELSE IF( WT == 6) THEN

			PhZe = 32._8 * ( Z1 ** 2 - 1._8) ** 14 * ( 31._8 * Z1 ** 2 - 1._8) 

        END IF
        
    END IF

    PhZe = 1._8 / ri ** 3 * PhZe


    ! End of Executive Part
    !=====================================================================
    END Subroutine DDWFun1D
    !=====================================================================
    
    !=====================================================================
    !*********************************************************************
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine SAlphaBeta( Alpha, Beta, VL)
    !=====================================================================
    ! Variables:
    
    !==============
    implicit none
    !==============
    
    Integer                             Alpha, Beta
	! Alpha, Beta

	Real(8)                             VL
	! Flag Value
    
    !=====================================================================
    
    VL = 0._8
    IF( Beta <= Alpha) THEN
        VL = 1._8
    END IF
    
    !=====================================================================
    END Subroutine SAlphaBeta
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
	Subroutine mteta( KiS, teta, PrX, DtX, TCi, dmax, PNUM, WT, VL)
	!=====================================================================
	! Variables:
    
    !==============
    implicit none
    !==============
    
	Integer                             I, J
	! Counters

	Integer                             PNUM
	! Number of Particles

	Integer                             teta
	! Order of m

	Integer                             WT
	! Window Function Type

	Real(8)                             Ci
	! Raduis of Influence For each Point

	Real(8)                             LC1, LC2 , LC3 , LC4  , LC5 , LC6
	Real(8)                             LC7, LC8 , LC9 , LC10 , LC11, LC12
	! Local Careers

	Real(8)                             dmax
	! Scaling Parameter

	Real(8)                             Phee
	! Returned Phee Value: Weight Function Value

	Real(8)                             VL
	! Mteta FinaL Value

	Real(8)                             KiS, PrC
	! Node Cordinate, Particle Cordinate

	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate

	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz

	Real(8),    DIMENSION( PNUM)    ::  TCi
	! Raduises For Particles
	
	!********************************************************************
	
	LC1 = 0._8
	LC2 = 0._8
    
    VL = 1._8

	IF( teta >= 0) THEN
        
        DO I = 1, PNUM

            PrC = PrX( I)
            Ci = TCi( I)
            
            CALL WFun1D( PrC, KiS, dmax, Ci, WT, Phee)
            
            LC1 = LC1 + ((KiS - PrC) ** teta) * Phee * DtX( I)
            
        END DO 
        
        VL = LC1     
        
    END IF
    
    
    !=====================================================================
    END Subroutine mteta
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine Dmteta( KiS, teta, PrX, DtX, TCi, dmax, PNUM, WT, VL)
    !=====================================================================
    
    !==============
    implicit none
    !==============
    
    Integer                             I, J
	! Counters
	 
	Integer                             PNUM
	! Number of Particles
	 
	Integer                             teta
	! Order of m
	 
	Integer                             WT
	! Window Function Type
	 
	Real(8)                             Ci
	! Raduis of Influence For Each Point
	 
	Real(8)                             LC1, LC2 , LC3 , LC4 , LC5 , LC6
	Real(8)                             LC7, LC8 , LC9 , LC10, LC11, LC12
	! Local Careers
	 
	Real(8)                             dmax
	! Scaling Parameter
	 
	Real(8)                             Phee, PhPr
	! Returned Phee Value; Weight Function Value
	! Returned PhPr Value; DWeight Function Value
	 
	Real(8)                             VL
	! DMteta Final Value
	 
	Real(8)                             KiS, PrC
	! Node Cordinate, Particle Cordinate
	 
	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate
	 
	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz
	 
	Real(8),    DIMENSION( PNUM)    ::  TCi
	! Raduises For Particles
	
	!********************************************************************
	
	LC1 = 0._8
	LC2 = 0._8
    
    VL = 0._8
	
	IF( teta >= 0) THEN
	
	    DO I = 1, PNUM
	    
	        PrC = PrX( I)
	        Ci = TCi( I)
	        
	        CALL WFun1D( PrC, KiS, dmax, Ci, WT, Phee)
	        
	        CALL DWFun1D( PrC, KiS, dmax, Ci, WT, PhPr)
	        
	        LC3 = KiS - PrC
	        
	        IF( teta == 0) THEN
	        
	            LC1 = LC1 + ((KiS - PrC) ** teta) * PhPr * DtX( I)
	            
	        ELSE
	        
	            LC1 = LC1 + (teta * (KiS - PrC) ** (teta - 1) * Phee * DtX( I)) + (((KiS - PrC) ** teta) * PhPr * DtX( I))
	            
	        END IF
	        
	    END DO
	    
	    VL = LC1
	    
	END IF
	
	!=====================================================================
    END Subroutine Dmteta
	!=====================================================================
    
    !*********************************************************************
	!=====================================================================
	Subroutine DDmteta( KiS, teta, PrX, DtX, TCi, dmax, PNUM, WT, VL)
	!=====================================================================
    
    !==============
    implicit none
    !==============
	
	Integer                             I, J
	! Counters
	 
	Integer                             PNUM
	! Number of Particles
	 
	Integer                             teta
	! Order of m
	 
	Integer                             WT
	! Window Function Type
	 
	Real(8)                             Ci
	! Raduis of Influence For Each Point
	 
	Real(8)                             LC1, LC2 , LC3 , LC4 , LC5 , LC6
	Real(8)                             LC7, LC8 , LC9 , LC10, LC11, LC12
	! Local Careers
	 
	Real(8)                             dmax
	! Scaling Parameter
	 
	Real(8)                             Phee, PhPr, PhZe
	! Returned Phee Value; Weight Function Value
	! Returned PhPr Value; DWeight Function Value
	 
	Real(8)                             VL
	! DMteta Final Value
	 
	Real(8)                             KiS, PrC
	! Node Cordinate, Particle Cordinate
	 
	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate
	 
	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz
	 
	Real(8),    DIMENSION( PNUM)    ::  TCi
	! Raduises For Particles
	
	!********************************************************************
	
	LC1 = 0._8	
	LC2 = 0._8
    
    VL = 0._8
	
	IF( teta >= 0) THEN
    
        DO I = 1, PNUM
        
            PrC = PrX( I)
            Ci = TCi( I)
            
            CALL WFun1D( PrC, KiS, dmax, Ci, WT, Phee)
            
            CALL DWFun1D( PrC, KiS, dmax, Ci, WT, PhPr)
            
            CALL DDWFun1D( PrC, KiS, dmax, Ci, WT, PhZe)
            
            LC3 = KiS - PrC
            
            IF( teta == 0) THEN
            
                LC1 = LC1 + ((KiS - PrC) ** teta) * PhZe * DtX( I)
                
            ELSE IF( teta == 1) THEN
            
                LC1 = LC1 + 2._8 * teta * (KiS - PrC) ** (teta - 1) * PhPr * DtX( I) + ((KiS - PrC) ** teta) * PhZe * DtX( I)
                
            ELSE
            
                LC1 = LC1 + 2._8 * teta * (KiS - PrC) ** (teta - 1) * PhPr * DtX( I) + ((KiS - PrC) ** teta) * PhZe * DtX( I) + teta * ( teta - 1) * (KiS - PrC) ** ( teta - 2) * Phee * DtX( I)
                
            END IF
            
        END DO
        
        VL = LC1
        
    END IF
    
    !=====================================================================
    END Subroutine DDmteta
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine momment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, MA)
    !=====================================================================
    
    !==============
    implicit none
    !==============
    
    Integer                             i, t, j
	! Counters
	
	INTEGER                             N
	! N
	 
	Integer                             teta
	! m Order
	 
	Integer                             eta
	! PolyNomial Order
	 
	Integer                             ILC1, ILC2
	! Integer Local Career
	 
	Integer                             ORD
	! Order of Polynomial

	Integer                             PNUM
	! Number of Particles
	 
	Integer                             WT
	! Window Function Type
	 
	Real(8)                             dmax
	! Scaling Parameter
	 
	Real(8)                             KiS
	! Node Cordinate
	 
	Real(8)                             S1, S2
	! Zarayeb
	 
	Real(8)                             VL
	! Returned Value
	 
	Real(8)                             LC1, LC2, LC3, LC4, LC5, LC6
	! Local Careers
	 
	ReaL(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate
	 
	ReaL(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz
	 
	ReaL(8),    DIMENSION( PNUM)    ::  TCi
	! Total Infulence Raduises
	 
	Real(8),    DIMENSION( ORD + 1, ORD + 1)  ::  MA
	! Momment Matrix
	
	!********************************************************************
	
	LC1 = 0._8
	LC2 = 0._8
	LC3 = 0._8
	LC4 = 0._8
	LC5 = 0._8
	LC6 = 0._8
	
	DO t = 1, ORD + 1
	
	    DO i = 1, ORD + 1
	    
	        LC6 = 0._8
	        
	        DO j = 0, eta
	        
	            ILC1 = t - 1
	            ILC2 = j
	            CALL SAlphaBeta( ILC1, ILC2, S1)
	        
	            ILC1 = i - 1
	            ILC2 = j
	            CALL SAlphaBeta( ILC1, ILC2, S2)
	            
	            N = i - 1
	            CALL DFAC( N, LC1)
	            
	            N = t - 1
	            CALL DFAC( N, LC2) 
	        
                N = i - 1 - J
                CALL DFAC( N, LC3)
	        
                N = t - 1 - j
                CALL DFAC( N, LC4)
	        
	            teta = i + t - 2 * j - 2
	        
	            CALL mteta( KiS, teta, PrX, DtX, TCi, dmax, PNUM, WT, VL)       
	        
	            LC6 = LC6 + S1 * S2 * ( LC1 * LC2) / ( LC3 * LC4) * VL
	        
	        END DO
	        
	        MA( t, i) = LC6
	           
	    END DO

	END DO
	
	!========================================================================
    END Subroutine momment
	!========================================================================
    
    !************************************************************************
    !========================================================================
	Subroutine Dmomment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, DMA)
	!========================================================================
    
    !==============
    implicit none
    !==============
	
	Integer                             i, t, j
	! Counters
	
	INTEGER                             N
	! N

	Integer                             teta
	! m Order

	Integer                             eta
	! PolyNomial Order

	Integer                             ILC1, ILC2
	! Integer Local Career

	Integer                             ORD
	! Order of PolyNomial

	Integer                             PNUM
	! Number of Particles

	Integer                             WT
	! Window Function Type

	Real(8)                             dmax
	! Scaling Parameter

	Real(8)                             KiS
	! Node Cordinate

	Real(8)                             S1, S2
	! Zarayeb

	Real(8)                             VL
	! Returned Value

	Real(8)                             LC1, LC2, LC3, LC4, LC5, LC6
	! Local Careers

	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate

	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz

	Real(8),    DIMENSION( PNUM)    ::  TCi
	! Total Influence Raduises

	Real(8),    DIMENSION( ORD + 1, ORD + 1)  ::  DMA
	! DMomment Matrix
	
	!********************************************************************
	
	LC1 = 0._8
	LC2 = 0._8
	LC3 = 0._8
	LC4 = 0._8
	LC5 = 0._8
	LC6 = 0._8
	
	DO t = 1, ORD + 1
	
	    DO i = 1, ORD + 1
	    
	        LC6 = 0._8
	        
	        DO j = 0, eta
	        
	            ILC1 = t - 1
		        ILC2 = j
		        CALL SAlphaBeta( ILC1, ILC2, S1)
		        
		        ILC1 = i - 1
		        ILC2 = J
		        CALL SalphaBeta( ILC1, ILC2, S2)
		        
		        N = i - 1
		        CALL DFAC( N, LC1)
		        
		        N = t - 1
		        CALL DFAC( N, LC2)
		        
                N = i - 1 - j
                CALL DFAC( N, LC3)
		        
                N = t - 1 - j
                CALL DFAC( N, LC4)
		        
		        teta = i + t - 2 * j - 2

		        CALL Dmteta( KiS, teta, PrX, DtX, TCi, dmax, PNUM, WT, VL)

		        LC6 = LC6 + S1 * S2 * ( LC1 * LC2) / ( LC3 * LC4) * VL
		        
		    END DO

		    DMA( t, i) = LC6

	    END DO

	END DO
	
	!=====================================================================
    END Subroutine Dmomment
	!=====================================================================
    
    !*********************************************************************
	!=====================================================================
	Subroutine DDmomment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, DDMA)
	!=====================================================================
    
    !==============
    implicit none
    !==============
	
	Integer                             i, t, j
	! Counters
	
	INTEGER                             N
	! N

	Integer                             teta
	! m Order

	Integer                             eta
	! PolyNomial Order

	Integer                             ILC1, ILC2
	! Integer Local Career

	Integer                             ORD
	! Order of PolyNomial

	Integer                             PNUM
	! Number of Particles

	Integer                             WT
	! Window Function Type

	Real(8)                             dmax
	! Scaling Parameter

	Real(8)                             KiS
	! Node Cordinate

	Real(8)                             S1, S2
	! Zarayeb

	Real(8)                             VL
	! Returned Value

	Real(8)                             LC1, LC2, LC3, LC4, LC5, LC6
	! Local Careers

	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate

	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz

	Real(8),    DIMENSION( PNUM)    ::  TCi
	! Total Influence Raduises

	Real(8),    DIMENSION( ORD + 1, ORD + 1)  ::  DDMA
	! DMomment Matrix
	
	!********************************************************************
	
	LC1 = 0._8
	LC2 = 0._8
	LC3 = 0._8
	LC4 = 0._8
	LC5 = 0._8
	LC6 = 0._8
	
	DO t = 1, ORD + 1
	
	    DO i = 1, ORD + 1
	    
	        LC6 = 0._8
	        
	        DO j = 0, eta
	        
	            ILC1 = t - 1
		        ILC2 = j
		        CALL SAlphaBeta( ILC1, ILC2, S1)
		        
		        ILC1 = i - 1
		        ILC2 = J
		        CALL SalphaBeta( ILC1, ILC2, S2)
		        
		        N = i - 1
		        CALL DFAC( N, LC1)
		        
		        N = t - 1
		        CALL DFAC( N, LC2)
		        
                N = i - 1 - j
                CALL DFAC( N, LC3)
		        
                N = t - 1 - j
                CALL DFAC( N, LC4)
		        
		        teta = i + t - 2 * j - 2
		        
		        CALL DDmteta( KiS, teta, PrX, DtX, TCi, dmax, PNUM, WT, VL)

		        LC6 = LC6 + S1 * S2 * ( LC1 * LC2) / ( LC3 * LC4) * VL
		        
		    END DO
		    
		    DDMA( t, i) = LC6
		    
		END DO
		
    END DO
    
    !=====================================================================
    END Subroutine DDmomment
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine GenRS1D( PNUM, SPNM, KiS, PrX, DtX, TCi, dmax, eta, ORD, WT, ZAY)
    !=====================================================================
    
    !==============
    implicit none
    !==============
    
    Integer                             I, II, J, JJ, K, KK
	Integer                             L, P, R
	! Counters

	Integer                             EQN, IPATH
	! Number of Equations, Path of Subroutine

	Integer                             LA, LB, LC
	! Leading Dimensions

	Integer                             ORD, eta
	! Order of Polynomials and Order of Shapes

	Integer                             PNUM, SPNM
	! Number of Particles, Specified Particles Number

	Integer                             R1, R2, R3
	! Local Integer Careers

	Integer                             WT
	! Window Function Type

	Real(8)                             dmax, Cri
	! Dilation Parameters

	Real(8)                             LC1, LC2, LC3, LC4, LC5, LC6
	! Local Careers

	Real(8)                             Phee
	! Window Function Types

	Real(8)                             KiS, PrC
	! Node Cordinate, Particle Cordinate

	Real(8),    ALLOCATABLE         ::  Base( :, :), BaseLU( :, :)
	! Base Vector

	Real(8),    ALLOCATABLE         ::  Bets( :, :)
	! Beta Vector

	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz

	Real(8),    ALLOCATABLE         ::  LCV11( :, :)
	! Local Career Matrix

	Real(8),    DIMENSION( PNUM)    ::  TCi
	! Total Raduises For Particles

	Real(8),    ALLOCATABLE         ::  MA( :, :), MALU( :, :)
    ! momment Matrix

	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles COrdinate

	Real(8),    ALLOCATABLE         ::  PT( :, :)
	! Polynomial Matrix

	Real(8),    DIMENSION( 1, 1: 1 + eta)  ::  ZAY
	! Shapes
	
	!***********************************************************************
	
	ALLOCATE( Base( 1: ORD + 1, 1))
	Base = 0._8
	Base( 1, 1) = 1._8
	!--------------------------------------------------------------------
	ALLOCATE( MA( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( Bets( 1: ORD + 1, 1))
	ALLOCATE( PT( 1: eta + 1, 1: ORD + 1))

	CALL momment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, MA)
	
	
	!--------------------------------------------------------------------
	ALLOCATE( MALU( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( BaseLU( 1: ORD + 1, 1))
	
	MALU = MA
	BaseLU = Base
	
	CALL GESV( MALU, BaseLU)
	
	Bets = BaseLU
	!---------------------------------------------------------------------	
	PrC = PrX( SPNM)
    
    PT = 0._8
	
	DO R1 = 0, eta
	
	    IF( R1 == 0) THEN
	    
	        DO R2 = 0, ORD	        
	            PT( R1 + 1, R2 + 1) = ( KiS - PrC) ** R2	            
	        END DO
	        
	    ELSE IF( R1 == 1) THEN
	    
	        DO R2 = 0, ORD
	        
	            IF( R2 == 0) THEN
	            
	                PT( R1 + 1, R2 + 1) = 0._8
	                
	            ELSE
	            
	                PT( R1 + 1, R2 + 1) = -R2 * ( KiS - PrC) ** ( R2 - 1)
	                
	            END IF
	            
	        END DO
	        
	    ELSE IF( R1 == 2) THEN
	    
	        DO R2 = 0, ORD
	        
	            IF( R2 <= 1) THEN
	            
	                PT( R1 + 1, R2 + 1) = 0._8
	                
	            ELSE
	            
	                PT( R1 +1, R2 + 1) = R2 * ( R2 - 1) * ( KiS - PrC) ** ( R2 - 2)
	                
	            END IF
	            
	        END DO
	        
	    ELSE IF( R1 == 3) THEN
	    
	        DO R2 = 0, ORD
	        
	            IF( R2 <= 2) THEN
	            
	                PT( R1 + 1, R2 + 1) = 0._8
	                
	            ELSE
	            
	                PT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( R2 - 2) * ( KiS - PrC) ** ( R2 - 3)
	            
	            END IF
	            
	        END DO
	        
	    END IF
	    
	END DO
	!--------------------------------------------------------------------
	Cri = TCi( SPNM)

	Call WFun1D( PrC, KiS, dmax, Cri, WT, Phee)

	LC1 = Phee * DtX( SPNM)

	DO I = 0, eta

	    LC3 = 0._8

	    DO J = 1, ORD + 1 

            LC3 = LC3 + Bets( J, 1) * PT( I + 1, J)

        END DO
        
        ZAY( 1, I + 1) = LC1 * LC3

    END DO
    !--------------------------------------------------------------------
    
    DEALLOCATE( Base)
	DEALLOCATE( MA)
	DEALLOCATE( Bets)
	DEALLOCATE( PT)
	DEALLOCATE( MALU)
	DEALLOCATE( BaseLU)
    
    !=====================================================================
    END Subroutine GenRS1D
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine DGenRS1D( PNUM, SPNM, KiS, PrX, DtX, TCi, dmax, eta, ORD, WT, DZAY)
    !=====================================================================
    
    !==============
    implicit none
    !==============
    
    Integer                             I, II, J, JJ, K, KK
	Integer                             L, P, R
	! Counters

	Integer                             EQN, IPATH
	! Number of Equations, Path of Subroutine

	Integer                             LA, LB, LC
	! Leading Dimension

	Integer                             ORD, eta
	! Order of PolyNomilas and Order of Shapes

	Integer                             PNUM, SPNM
	! Number of Particles, Specified Particles Number

	Integer                             R1, R2, R3
	! Local Integer Career

	Integer                             WT
	! Window Function Type

	Real(8)                             ALF, Bet
	! Alpha and Beta

	Real(8)                             dmax, Cri
	! Dilation Parameter

	Real(8)                             LC1, LC2, LC3, LC4, LC5, LC6
	! Local Careers

	Real(8)                             Phee, PhPr
	! Values

	Real(8)                             KiS, PrC
	! Node Cordinate, Particle Cordinate

	Real(8),    ALLOCATABLE         ::  Base( :, :), BaseLU( :, :)
	! Base Vector

	Real(8),    ALLOCATABLE         ::  Bets( :, :), Kets( :, :)
	! Beta Vector

	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz

	Real(8),    ALLOCATABLE         ::  LCV11( :, :), LCV12( :, :), LCV13( :, :)
	Real(8),    ALLOCATABLE         ::  LCV1( :, :), LCV1LU( :, :)
	! Local Career Matrix

	Real(8),    DIMENSION( PNUM)    ::  TCi
	! TotaL Raduises for Particles

	Real(8),    ALLOCATABLE         ::  MA( :, :), DMA( :, :), MALU( :, :)
	! momment Matrix, Dmomment Matrix

	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate

	Real(8),    ALLOCATABLE         ::  PT( :, :), DPT( :, :)
	! Polynomial Matrix

	Real(8),    DIMENSION( 1, 1: 1 + eta)  ::  DZAY
	! Shapes
	
	!********************************************************************
	
	ALLOCATE( Base( 1: ORD + 1, 1))
	Base = 0._8
	Base( 1, 1) = 1._8
	!--------------------------------------------------------------------
	ALLOCATE(  MA( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( DMA( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( Bets( 1: ORD + 1, 1))
	ALLOCATE( Kets( 1: ORD + 1, 1))
	ALLOCATE(  PT( 1: eta + 1, 1: ORD + 1))
	ALLOCATE( DPT( 1: eta + 1, 1: ORD + 1))
	ALLOCATE( LCV1( 1: ORD + 1, 1))

	CALL momment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, MA)

	CALL Dmomment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, DMA)
	!--------------------------------------------------------------------

	!--------------------------------------------------------------------
	ALLOCATE( MALU( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( BaseLU( 1: ORD + 1, 1))
	
	MALU = MA
	BaseLU = Base
	
	CALL GESV( MALU, BaseLU)
	
	Bets = BaseLU
	!---------------------------------------------------------------------
	
	!---------------------------------------------------------------------
	CALL GEMM( DMA, Bets, LCV1, 'N', 'N')
	LCV1 = -LCV1
	! MP Matrix Multiplied By Bets Vector
	!---------------------------------------------------------------------
	
	!---------------------------------------------------------------------
	ALLOCATE( LCV1LU( 1: ORD + 1, 1))
	
	MALU = MA
	LCV1LU = LCV1
	
	CALL GESV( MALU, LCV1LU)
	
	Kets = LCV1LU
	! Kets Vector Computed
	!--------------------------------------------------------------------
	
	!--------------------------------------------------------------------
	PrC = PrX( SPNM)

	DO R1 = 0, eta
	
	    IF( R1 == 0) THEN
	    
	        DO R2 = 0, ORD
	        
	            PT( R1 + 1, R2 + 1) = ( KiS - PrC) ** R2
	            
	            IF( R2 == 0 ) THEN
	                DPT( R1 + 1, R2 + 1) = 0._8
	            ELSE
	                DPT( R1 + 1, R2 + 1) = R2 * ( KiS - PrC) ** ( R2 - 1)
	            END IF
	            
            END DO
            
        ELSE IF( R1 == 1) THEN
        
            DO R2 = 0, ORD
            
                IF( R2 == 0) THEN
                    PT( R1 + 1, R2 + 1) = 0._8
                ELSE
                    PT( R1 + 1, R2 + 1) = -R2 * ( KiS - PrC) ** ( R2 - 1)
                END IF
                
                IF( R2 == 0 .OR. R2 == 1) THEN
                    DPT( R1 + 1, R2 + 1) = 0._8
                ELSE
                    DPT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( KiS - PrC) ** ( R2 - 2)
                END IF
                
            END DO
            
        ELSE IF( R1 == 2) THEN
        
            DO R2 = 0, ORD
            
                IF( R2 <= 1) THEN
                    PT( R1 + 1, R2 + 1) = 0._8
                ELSE
                    PT( R1 + 1, R2 + 1) = R2 * ( R2 - 1) * ( KiS - PrC) ** ( R2 - 2)
                END IF
                
                IF( R2 <= 2) THEN
                    DPT( R1 + 1, R2 + 1) = 0._8
                ELSE
                    DPT( R1 + 1, R2 + 1) = R2 * ( R2 - 1) * ( R2 - 2) * ( KiS - PrC) ** ( R2 - 3)
                END IF
                
            END DO
            
        ELSE IF( R1 == 3) THEN
        
            DO R2 = 0, ORD
            
                IF( R2 <= 2) THEN
                    PT( R1 + 1, R2 + 1) = 0._8
                ELSE
                    PT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( R2 - 2) * ( KiS - PrC) ** ( R2 - 3)
                END IF
                
                IF( R2 <= 3) THEN
                    DPT( R1 + 1, R2 + 1) = 0._8
                ELSE
                    DPT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( R2 - 2) * ( R2 - 3) * ( KiS - PrC) ** ( R2 - 4)
                END IF
                
            END DO
            
        END IF
        
    END DO
    !--------------------------------------------------------------------
    
    !------------------------------------------------------------------
	Cri = TCi( SPNM)

	CALL WFun1D( PrC, KiS, dmax, CrI, WT, Phee)

	CALL DWFun1D( PrC, KiS, dmax, Cri, WT, PhPr)

	LC3 = 0._8
	LC5 = 0._8

	DO I = 0, eta

	    LC3 = 0._8
	    LC5 = 0._8

	    DO J = 1, ORD + 1

            LC2 = Kets( J, 1) * PT( I + 1, J) + Bets( J, 1) * DPT( I + 1, J)
            LC3 = LC3 + LC2

            LC4 = Bets( J, 1) * PT( I + 1, J)
            LC5 = LC5 + LC4

        END DO
        
        DZAY( 1, I + 1) = ( Phee * LC3 + PhPr * LC5) * DtX( SPNM)

    END DO
    !--------------------------------------------------------------------
    
    !--------------------------------------------------------------------
    
    DEALLOCATE( Base)
	DEALLOCATE(  MA)
	DEALLOCATE( DMA)
	DEALLOCATE( Bets)
	DEALLOCATE( Kets)
	DEALLOCATE(  PT)
	DEALLOCATE( DPT)
	DEALLOCATE( LCV1)
	DEALLOCATE( MALU)
	DEALLOCATE( BaseLU)
	DEALLOCATE( LCV1LU)
    
    !=====================================================================
    END Subroutine DGenRS1D
    !=====================================================================
    
    !*********************************************************************
    !=====================================================================
    Subroutine DDGenRS1D( PNUM, SPNM, KiS, PrX, DtX, TCi, dmax, eta, ORD, WT, DDZAY)
    !=====================================================================
    
    !==============
    implicit none
    !==============
    
    Integer                             I, II, J, JJ, K, KK
	Integer                             L, P, R
	! Counters

	Integer                             EQN, IPATH
	! Number of Equations, Path of Subroutine

	Integer                             LA, LB, LC
	! Leading Dimension

	Integer                             ORD, eta
	! Order of PolyNomilas and Order of Shapes

	Integer                             PNUM, SPNM
	! Number of Particles, Specified Particles Number

	Integer                             R1, R2, R3
	! Local Integer Career

	Integer                             WT
	! Window Function Type

	Real(8)                             ALF, Bet
	! Alpha and Beta

	Real(8)                             dmax, Cri
	! Dilation Parameter

	Real(8)                             LC1, LC2, LC3, LC4, LC5, LC6
	Real(8)                             LC7, LC8, LC9, LC10
	! Local Careers

	Real(8)                             Phee, PhPr, PhZe
	! Values

	Real(8)                             KiS, PrC
	! Node Cordinate, Particle Cordinate

	Real(8),    ALLOCATABLE         ::  Base( :, :), BaseLU( :, :)
	! Base Vector

	Real(8),    ALLOCATABLE         ::  Bets( :, :), Kets( :, :), Xets( :, :)
	! Beta Vector

	Real(8),    DIMENSION( PNUM)    ::  DtX
	! Efraz

	Real(8),    ALLOCATABLE         ::  LCV11( :, :), LCV12( :, :), LCV13( :, :)
	Real(8),    ALLOCATABLE         ::  LCV1( :, :), LCV2( :, :), LCV3( :, :), LCV1LU( :, :)
	! Local Career Matrix

	Real(8),    DIMENSION( PNUM)    ::  TCi
	! TotaL Raduises for Particles

	Real(8),    ALLOCATABLE         ::  MA( :, :), DMA( :, :), DDMA( :, :), MALU( :, :)
	! momment Matrix, Dmomment Matrix

	Real(8),    DIMENSION( PNUM)    ::  PrX
	! Particles Cordinate

	Real(8),    ALLOCATABLE         ::  PT( :, :), DPT( :, :), DDPT( :, :)
	! Polynomial Matrix

	Real(8),    DIMENSION( 1, 1: 1 + eta)  ::  DDZAY
	! Shapes
	
	!********************************************************************
	
	ALLOCATE( Base( 1: ORD + 1, 1))
	Base = 0._8
	Base( 1, 1) = 1._8
	!--------------------------------------------------------------------
	ALLOCATE(   MA( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE(  DMA( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( DDMA( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( Bets( 1: ORD + 1, 1))
	ALLOCATE( Kets( 1: ORD + 1, 1))
	ALLOCATE( Xets( 1: ORD + 1, 1))
	ALLOCATE(   PT( 1: eta + 1, 1: ORD + 1))
	ALLOCATE(  DPT( 1: eta + 1, 1: ORD + 1))
	ALLOCATE( DDPT( 1: eta + 1, 1: ORD + 1))
	ALLOCATE( LCV1( 1: ORD + 1, 1))
	ALLOCATE( LCV2( 1: ORD + 1, 1))
	ALLOCATE( LCV3( 1: ORD + 1, 1))

    CALL momment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, MA)

	CALL Dmomment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, DMA)

	CALL DDmomment( KiS, PrX, DtX, PNUM, ORD, eta, TCi, dmax, WT, DDMA)
	!--------------------------------------------------------------------
	
	!--------------------------------------------------------------------
	ALLOCATE( MALU( 1: ORD + 1, 1: ORD + 1))
	ALLOCATE( BaseLU( 1: ORD + 1, 1))
	
	MALU = MA
	BaseLU = Base
	
	CALL GESV( MALU, BaseLU)
	
	Bets = BaseLU
	!---------------------------------------------------------------------
	
	!---------------------------------------------------------------------
	CALL GEMM( DMA, Bets, LCV1, 'N', 'N')
	LCV1 = -LCV1
	! MP Matrix Multiplied By Bets Vector
	!---------------------------------------------------------------------
	
	!---------------------------------------------------------------------
	ALLOCATE( LCV1LU( 1: ORD + 1, 1))
	
	MALU = MA
	LCV1LU = LCV1
	
	CALL GESV( MALU, LCV1LU)
	
	Kets = LCV1LU
	DEALLOCATE( LCV1LU)
	! Kets Vector Computed
	!--------------------------------------------------------------------
	
	!--------------------------------------------------------------------
	CALL GEMM( DDMA, Bets, LCV2, 'N', 'N')
	LCV2 = -LCV2
	! DDMA Matrix Multiplied By Bets Vector
	!--------------------------------------------------------------------
	
	!--------------------------------------------------------------------
	CALL GEMM( DMA, Kets, LCV3, 'N', 'N')
	LCV3 = -2._8 * LCV3
	! DMA Matrix Multiplied By Kets Vector
	!--------------------------------------------------------------------
	
	!--------------------------------------------------------------------
	LCV1 = LCV2 + LCV3
	!--------------------------------------------------------------------
	
	!--------------------------------------------------------------------
	ALLOCATE( LCV1LU( 1: ORD + 1, 1))
	
	MALU = MA
	LCV1LU = LCV1
	
	CALL GESV( MALU, LCV1LU)
	
	Xets = LCV1LU
	! Xets Vector Computed
	!--------------------------------------------------------------------
	
	!--------------------------------------------------------------------
	PrC = PrX( SPNM)
	
	DO R1 = 0, eta
	
	    IF( R1 == 0) THEN
	    
	        DO R2 = 0, ORD
	        
	            PT( R1 + 1, R2 + 1) = ( KiS - PrC) ** R2
	            
	            IF( R2 == 0 ) THEN
		            DPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DPT( R1 + 1, R2 + 1) = R2 * ( KiS - PrC) ** ( R2 - 1)
		        END IF
		        
		        IF( R2 <= 1) THEN
		            DDPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DDPT( R1 + 1, R2 + 1) = R2 * ( R2 - 1) * ( KiS - PrC) ** ( R2 - 2)
		        END IF
		        
		    END DO
		    
		ELSE IF( R1 == 1) THEN

	        DO R2 = 0, ORD

		        IF( R2 == 0) THEN
		            PT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            PT( R1 + 1, R2 + 1) = -R2 * ( KiS - PrC) ** ( R2 - 1)
                END IF

		        IF( R2 == 0 .OR. R2 == 1) THEN
		            DPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DPT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( KiS - PrC) ** ( R2 - 2)
		        END IF

		        IF( R2 <= 2) THEN
		            DDPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DDPT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( R2 - 2) * ( KiS - PrC) ** ( R2 - 3)
		        END IF

		    END DO
		    
		ELSE IF( R1 == 2) THEN

	        DO R2 = 0, ORD

		        IF( R2 <= 1) THEN
                    PT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            PT( R1 + 1, R2 + 1) = R2 * ( R2 - 1) * ( KiS - PrC) ** ( R2 - 2)
		        END IF

		        IF( R2 <= 2) THEN
		            DPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DPT( R1 + 1, R2 + 1) = R2 * ( R2 - 1) * ( R2 - 2) * ( KiS - PrC) ** ( R2 - 3)
		        END IF

		        IF( R2 <= 3) THEN
		            DDPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DDPT( R1 + 1, R2 + 1) = R2 * ( R2 - 1) * ( R2 - 2) * ( R2 - 3) * ( KiS - PrC) ** ( R2 - 4)
		        END IF

		    END DO
		    
		ELSE IF( R1 == 3) THEN

	        DO R2 = 0, ORD

		        IF( R2 <= 2) THEN
		            PT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            PT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( R2 - 2) * ( KiS - PrC) ** ( R2 - 3)
		        END IF

		        IF( R2 <= 3) THEN
		            DPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DPT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( R2 - 2) * ( R2 - 3) * ( KiS - PrC) ** ( R2 - 4)
		        END IF

		        IF( R2 <= 4) THEN
		            DDPT( R1 + 1, R2 + 1) = 0._8
		        ELSE
		            DDPT( R1 + 1, R2 + 1) = -R2 * ( R2 - 1) * ( R2 - 2) * ( R2 - 3) * ( R2 - 4) * ( KiS - PrC) ** ( R2 - 5)
		        END IF

	        END DO

	    END IF

	END DO
	!--------------------------------------------------------------------

	!--------------------------------------------------------------------
	Cri = TCi( SPNM)

	CALL WFun1D( PrC, KiS, dmax, CrI, WT, Phee)

	CALL DWFun1D( PrC, KiS, dmax, Cri, WT, PhPr)

	CALL DDWFun1D( PrC, KiS, dmax, Cri, WT, PhZe)

	LC3 = 0._8
	LC5 = 0._8
	LC7 = 0._8

	DO I = 0, eta

	    LC3 = 0._8
	    LC5 = 0._8
	    LC7 = 0._8

	    DO J = 1, ORD + 1

		    LC2 = Xets( J, 1) * PT( I + 1, J) + 2._8 * Kets( J, 1) * DPT( I + 1, J) + Bets( J, 1) * DDPT( I + 1, J)
		    LC3 = LC3 + LC2

		    LC4 = Kets( J, 1) * PT( I + 1, J) + Bets( J, 1) * DPT( I + 1, J)
		    LC5 = LC5 + LC4

		    LC6 = Bets( J, 1) * PT( I + 1, J)
		    LC7 = LC7 + LC6

	    END DO

	    DDZAY( 1, I + 1) = ( Phee * LC3 + 2._8 * PhPr * LC5 + PhZe * LC7) * DtX( SPNM)

	END DO

	!--------------------------------------------------------------------
	
	DEALLOCATE( Base)
	DEALLOCATE(   MA)
	DEALLOCATE(  DMA)
	DEALLOCATE( DDMA)
	DEALLOCATE( Bets)
	DEALLOCATE( Kets)
	DEALLOCATE( Xets)
	DEALLOCATE(   PT)
	DEALLOCATE(  DPT)
	DEALLOCATE( DDPT)
	DEALLOCATE( LCV1)
	DEALLOCATE( LCV2)
	DEALLOCATE( LCV3)
    DEALLOCATE( MALU)
	DEALLOCATE( BaseLU)
	DEALLOCATE( LCV1LU)
	
	!********************************************************************
	END Subroutine DDGenRS1D
	!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	!55555555555555555555555555555555555555555555555555555555555555555555
    

	 END Module RKPM1D_Intel