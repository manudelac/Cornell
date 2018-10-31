	!Manuel de la Cruz González DNI:70909708H, Compañero de grupo: Javier García Marcos
	PROGRAM ESPECTRO DEL BOTONOMIO

        INTEGER,PARAMETER :: n=40 ! adimensional	                    
        DOUBLE PRECISION,PARAMETER :: r1=0.05d0 ! fm
	INTEGER,PARAMETER :: nmax=40 ! adimensional
	DOUBLE PRECISION,PARAMETER :: rmax=8.d0 ! fm
	DOUBLE PRECISION,PARAMETER :: pi=3.141592d0 ! adimensioanl
	DOUBLE PRECISION,PARAMETER :: hc=0.197327d0 ! GeVfm
	DOUBLE PRECISION,PARAMETER :: mb=4.733d0 ! GeV/c²
	DOUBLE PRECISION,PARAMETER :: sigma=0.207d0 ! GeV²
	DOUBLE PRECISION,PARAMETER :: alfas=0.366d0 ! adimensional
	INTEGER,PARAMETER :: nr=3 !adim
	DOUBLE PRECISION rn(n),etan(n),nnl(n,0:n),betan(n,n)
	DOUBLE PRECISION ig(n,n),igt(n,n),igv1(n,n)
        DOUBLE PRECISION igv2(n,n),s(n,n),ss(n,n),t(n,n) 
	DOUBLE PRECISION ga,vn(n,n),ham(n,n)                                          
	DOUBLE PRECISION alphar(n),alphai(n),beta(n)
	DOUBLE PRECISION z(n,n),av(n,n),work(n)		       
	DOUBLE PRECISION aut(n),minim1,minim2,minim3
	DOUBLE PRECISION av1(n),av2(n),av3(n)
	DOUBLE PRECISION suma1,suma2,suma3,r,limsup,liminf 
	DOUBLE PRECISION fi,norma(nr),difmas,difmas1                 
	INTEGER i,j,k,l,info,aa,bb,cc,imin(nr)
	CHARACTER V

	liminf=0.d0 !fm
	limsup=1.5d0 !fm

c APERTURA DE FICHEROS PARA ESCRIBIR LOS VALORES
c ______________________________________________

c	OPEN(10,FILE='nnl0.dat')   ! fichero para la matriz N(n,l) con l=0
c       OPEN(11,FILE='nnl1.dat')   ! fichero para la matriz N(n,l) con l=1
c       OPEN(12,FILE='ig.dat')     ! fichero para la matriz ig(n,n) con l=0,1
c	OPEN(13,FILE='S0.dat')     ! fichero para la matriz S(n,n) con l=0
c  	OPEN(14,FILE='S1.dat')     ! fichero para la matriz S(n,n) con l=1
c       OPEN(15,FILE='igt.dat')    ! fichero para la matriz igt(n,n) con l=0,1
c	OPEN(16,FILE='T0.dat')     ! fichero para la matriz T(n,n) con l=0
c  	OPEN(17,FILE='T1.dat')     ! fichero para la matriz T(n,n) con l=1
c       OPEN(18,FILE='igv1.dat')   ! fichero para la matriz igv1(n,n) con l=0,1
c       OPEN(19,FILE='igv2.dat')   ! fichero para la matriz igv2(n,n) con l=0,1
c	OPEN(20,FILE='V0.dat')     ! fichero para la matriz V(n,n) con l=0
c  	OPEN(21,FILE='V1.dat')     ! fichero para la matriz V(n,n) con l=1
c	OPEN(22,FILE='H0.dat')     ! fichero para la matriz H(n,n) con l=0             
c	OPEN(23,FILE='H1.dat')	   ! fichero para la matriz H(n,n) con l=1		

c	OPEN(24,FILE='aut0.dat')   ! fichero para escribir los autovalores con l=0			
c	OPEN(25,FILE='aut1.dat')   ! fichero para escribir los autovalores con l=1		
c	OPEN(26,FILE='av0.dat')	   ! fichero para escribir los autovectores con l=0		
c	OPEN(27,FILE='av1.dat')	   ! fichero para escribir los autovectores con l=1		

c	OPEN(28,FILE='3av0.dat')    ! fichero para escribir los autovectores de mínima energía con l=0			
c	OPEN(29,FILE='3av1.dat')    ! fichero para escribir los autovectores de mínima energía con l=1			
c	OPEN(30,FILE='WaveFn1l0.dat')		!Ficheros para escribir las funciones de onda para la representación
c	OPEN(31,FILE='WaveFn2l0.dat')		
c	OPEN(32,FILE='WaveFn3l0.dat')
c	OPEN(33,FILE='WaveFn1l1.dat')
c	OPEN(34,FILE='WaveFn2l1.dat')
c	OPEN(35,FILE='WaveFn3l1.dat')
	OPEN(36,FILE='wavefun.dat')		!Fichero de salida de las funciones de onda
	


c PROGRAMA
c ________

	DO j=0,1			!Este bucle j recorre todo el programa dando valor l=0,1

	 DO i=1,n			!Funciones para calcular Nnl
c rn
	  rn(i)=0.d0

	  rn(i)=r1*((rmax/r1)**((i-1.d0)/
     &    (nmax-1.d0)))                               ! cálculo del vector r(n) necesario para calcular eta(n)
c Etan
          etan(i)=0.d0

	  etan(i)=1.d0/((rn(i))**(2.d0))               ! cálculo de eta(n) necesario para hallar los valores de la norma 							        N(n,l)
c Nnl
c================================================
	  nnl(i,j)=0.d0                                         ! inicialización de N(n,l) a CERO

	 nnl(i,j)=dsqrt(((2.d0**(j+2.d0))*	!Cálculo de Nnl
     &   ((2.d0*etan(i))**(j+1.5d0)))/
     &   ((dsqrt(pi))*(2.d0*j+1.d0)))
	
c	   IF (j==0) THEN                    ! bucle IF para escribir los valores de N(n,l) según sea l=0 ó l=1
c	    WRITE(10,*)(nnl(i,j)),i,j
c	     ELSE
c	     WRITE(11,*)(nnl(i,j)),i,j
c	   END IF
	 END DO

c CÁLCULO DE BETA(n,n)
c ____________________

	 DO k=1,n
	  DO l=1,n

	   betan(k,l)=0.d0              ! inicialización de beta(n,n) a CERO

	   betan(k,l)=etan(k)+etan(l)

          END DO
	 END DO

c CALCULO DE LA INTEGRAL COMÚN DE S(n,n;l=0,1), T(n,n;l=0,1) y V(n,n;l=0,1); I(2l,beta)
c _____________________________________________________________________________________

	 CALL GammaFun(((2.d0*j+3.d0)/2.d0),GA)     ! llamada a la subrutina que calcula valores de la función GAMMA
 	 DO k=1,n
	  DO l=1,n

	   ig(k,l)=0.d0                                                      ! inicialización de la integral a CERO

	   ig(k,l)=0.5d0*((betan(k,l))**((-1.d0)*
     &     ((2.d0*j+3.d0)/2.d0)))*GA                                         ! cálculo de la integral
c	        WRITE(12,*)ig(k,l),k,l,j                                         ! escritura en fichero para comprobar
	  END DO
         END DO

c CÁLCULO DE LA MATRIZ DE SOLAPAMIENTO S(n,n) PARA l=0,1
c ______________________________________________________

	 DO k=1,n
	  DO l=1,n

	   s(k,l)=0.d0                         ! inicialización de la matriz a CERO
	   ss(k,l)=0.d0
	   s(k,l)=nnl(k,j)*(nnl(l,j))*ig(k,l)
	   ss(k,l)=s(k,l)			!CLONO S PARA USARLA POSTERIORMENTE YA QUE DDGEV LA MODIFICA
	   					    
c	   IF (j==0) THEN 			    ! bucle IF para que escriba los valores de la matriz según sea l=0 ó 								l=1
c	   WRITE(13,*)(s(k,l)),k,l,j                ! escritura en fichero para S(n,n;l=0)
c	    ELSE
c	     WRITE(14,*)(s(k,l)),k,l,j              ! escritura en fichero para S(n,n;l=1)
c	   END IF
	  END DO
         END DO

c CALCULO DE LA INTEGRAL PARA T(n,n;l=0,1); I(2l+2,beta)
c ______________________________________________________

	CALL GammaFun(((2.d0*j+5.d0)/2.d0),GA)              ! llamada a la subrutina que calcula valores de la función GAMMA
 	 DO k=1,n
	  DO l=1,n

	   igt(k,l)=0.d0                                     ! inicialización de la integral a CERO

	   igt(k,l)=0.5d0*((betan(k,l))**((-1.d0)*
     &     ((2.d0*j+5.d0)/2.d0)))*GA     		     ! cálculo de la integral
c	   WRITE(15,*)igt(k,l),k,l,j                         ! escritura en fichero para comprobar
	  END DO
	 END DO

c CÁLCULO DE LA MATRIZ ENERGÍA CINÉTICA T(n,n) PARA l=0,1
c _______________________________________________________

	DO k=1,n
	 DO l=1,n

	  t(k,l)=0.d0                                             ! cálculo de la matriz T(n,n)

          t(k,l)=(-1.d0)*(((hc)**2.d0)/mb)                          ! inicialización de la matriz a CERO
     &    *(nnl(k,j))*(nnl(l,j))     
     &    *(4.d0*(etan(l)**2.d0)*
     &    igt(k,l)-2.d0*etan(l)
     &    *(2.d0*j+3.d0)*ig(k,l))
c	   IF (j==0) THEN                                         ! bucle IF para que escriba valores de la matriz según 										sea l=0 ó l=1
c	    WRITE(16,*)(t(k,l)),k,l,j                             ! escritura en fichero para T(n,n;l=0)
c	     ELSE
c	      WRITE(17,*)(t(k,l)),k,l,j                           ! escritura en fichero para T(n,n;l=1)
c 	   END IF
	 END DO
	END DO

c CÁLCULO DE LA 1ª INTEGRAL PARA V(n,n;l=0,1); I(2l-1,beta)
c ______________________________________________________

	CALL GammaFun(((2.d0*j+2.d0)/2.d0),GA)              ! llamada a la subrutina que calcula valores de la función GAMMA
 	 DO k=1,n
	  DO l=1,n

	   igv1(k,l)=0.d0                                                       ! inicialización de la integral a CERO

	   igv1(k,l)=0.5d0
     &     *((betan(k,l))**((-1.d0)*
     &     ((2.d0*j+2.d0)/2.d0)))*GA     ! cálculo de la integral
c	   WRITE(18,*)igv1(k,l),k,l,j                                           ! escritura en fichero para comprobar
	  END DO
	 END DO

c CÁLCULO DE LA 2ª INTEGRAL PARA V(n,n;l=0,1); I(2l+1,beta)
c _________________________________________________________

	CALL GammaFun(((2.d0*j+4.d0)/2.d0),GA)              ! llamada a la subrutina que calcula valores de la función GAMMA
 	 DO k=1,n
	  DO l=1,n

	   igv2(k,l)=0.d0                             ! inicialización de la integral a CERO
                                                    
	   igv2(k,l)=0.5d0*
     &     ((betan(k,l))**((-1.d0)*
     &     ((2.d0*j+4.d0)/2.d0)))*GA     ! cálculo de la integral
c	   WRITE(19,*)igv2(k,l),k,l,j                                           ! escritura en fichero para comprobar
	  END DO
	 END DO

c CÁLCULO DE LA MATRIZ POTENCIAL V(n,n) PARA l=0,1
c _______________________________________________________

	DO k=1,n
	 DO l=1,n

	  vn(k,l)=0.d0                                             ! cálculo de la matriz V(n,n)
	                                                           ! inicialización de la matriz a CERO

          vn(k,l)=(nnl(k,j))*(nnl(l,j))*
     &    (((-4.d0/3.d0)*(alfas)*(hc)*
     &    igv1(k,l))+(((sigma)/(hc))
     &    *igv2(k,l))+(2.d0*(mb)*ig(k,l)))

c	   IF (j==0) THEN                                         ! bucle IF para que escriba valores de la matriz según 									   sea l=0 ó l=1
c	    WRITE(20,*)(vn(k,l)),k,l,j                             ! escritura en fichero para V(n,n;l=0)
c	     ELSE
c	      WRITE(21,*)(vn(k,l)),k,l,j                           ! escritura en fichero para V(n,n;l=1)
c 	   END IF
	 END DO
	END DO


c CÁLCULO DE LA MATRIZ HAMILTONIANA H(n,n) PARA l=0,1
c__________________________________________________________

	DO k=1,n
	 DO l=1,n

	  ham(k,l)=0.d0                                             ! inicialización de la matriz a CERO

	  ham(k,l)=t(k,l)+vn(k,l)
c	  IF (j==0) THEN                                            ! bucle IF para que escriba valores de la matriz según sea l=0 ó l=1
c	    WRITE(22,*)(ham(k,l)),k,l,j                             ! escritura en fichero para H(n,n;l=0)
c	     ELSE
c	      WRITE(23,*)(ham(k,l)),k,l,j                           ! escritura en fichero para H(n,n;l=1)
c	  END IF
         END DO
        END DO




c SUBRUTINA DGGEV QUE RESUELVE EL PROBLEMA DE AUTOVALORES PARA l=0,1
c=========================================================================================
	CALL DGGEV('N', 'V', n, ham, n, s, 
     &             n, ALPHAR, ALPHAI,
     &             BETA, z, n, av, n, WORK, 16*n, INFO)

c Salida en pantalla de INFO
	 IF (j==0) THEN  
	  WRITE(*,*)              
	  WRITE(*,*)"Variable info de la subrutina DGGEV para l=0"
	  WRITE(*,*)INFO
	  WRITE(*,*)                        
	 ELSE
	  WRITE(*,*)"Variable info de la subrutina DGGEV para l=1"
	  WRITE(*,*)INFO
	  WRITE(*,*)  
	 END IF

	DO k=1,n

	 aut(k)=0.d0

	 aut(k)=ALPHAR(k)/BETA(k)              !Calculamos los AUTOVALORES según la definición de DGGEV

c	 IF (j==0) THEN                                   ! bucle IF escribe valores del vector según sea l=0 ó l=1 
c	  WRITE(24,*)aut(k),k,j                           ! escritura en fichero para aut(n;l=0)
c	   ELSE
c	  WRITE(25,*)aut(k),k,j                           ! escritura en fichero para aut(n;l=1)
c	 END IF
c        DO l=1,n
c	  IF (j==0) THEN                                   ! bucle IF para que escriba autovectores según sea l=0 ó l=1
c	   WRITE(26,*)av(k,l),k,l,j                           ! escritura en fichero para av(n,n;l=0)
c 	    ELSE
c	   WRITE(27,*)av(k,l),k,l,j                           ! escritura en fichero para av(n,n;l=1)
c	  END IF
c        END DO         
	END DO
                 


c BUCLE PARA ORDENAR LOS 3 AUTOVALORES MÁS BAJOS PARA l=0,1
c__________________________________________________________

	minim1=dabs(1.d0+aut(1))	!Creamos una variable para entrar en el bucle que compare el valor...
	minim2=dabs(1.d0+aut(1))	!...de los autovalores y escoja los 3 menores
	minim3=dabs(1.d0+aut(1))


 	DO k=1,n
	 IF (minim1.GT.aut(k)) THEN
          minim1=aut(k)
          aa=k				!Las variables aa,bb y cc nos guardan el índice de los autovalores más bajos...
          aut(k)=aut(1)  		!... para saber donde están almacenados en memoria
          aut(1)=minim1
         END IF
	END DO

	DO k=2,n
	 IF (minim2.GT.aut(k)) THEN
          minim2=aut(k)
          bb=k

          aut(k)=aut(2)  
          aut(2)=minim2
         END IF
	END DO

	DO k=3,n
	 IF (minim3.GT.aut(k)) THEN
          minim3=aut(k)
          cc=k

          aut(k)=aut(3)  
          aut(3)=minim3

         END IF
        END DO
					
c ESCRIBIMOS EN PANTALLA LAS ENERGÍAS MÁS BAJAS      
	 IF (j==0) THEN                
	  WRITE(*,*)"Energias mas bajas n=1,2,3 (respectivamente) para ond
     &a S (l=0)"
	  WRITE(*,"(F15.10)")aut(1),aut(2),aut(3) 
	  WRITE(*,*)                         
	 ELSE
	  WRITE(*,*)"Energias mas bajas n=1,2,3 (respectivamente) para ond
     &a P (l=1)"
	  WRITE(*,"(F15.10)")aut(1),aut(2),aut(3)
	  WRITE(*,*)  
	 END IF
	                        	

c	WRITE(*,*)aa,bb,cc


c BUCLES PARA ESCOGER Y NORMALIZAR LOS AUTOVECTORES DE LAS TRES ENERGÍAS MAS BAJAS PARA l=0,1
c____________________________________________________________________________________________

	DO k=1,n                !escojo autovectores con la posicion aa bb cc, la posición la hemos sacado de las 
	 av1(k)=av(k,aa)	!autoenergías mínimas
	 av2(k)=av(k,bb)
	 av3(k)=av(k,cc)
	END DO

c       DO k=1,n		!Compruebo escribiendo en archivo
c	 IF (j==0) THEN                
c	  WRITE(28,*)av1(k),av2(k),av3(k),k,j                           
c	 ELSE
c	  WRITE(29,*)av1(k),av2(k),av3(k),k,j
c	 END IF
c	END DO                          	

c	NORMALIZACIÓN
c==========================================================================

c Sumatorio previo a la normalización

	suma1=0.d0		!Normalizo cada vector propio de menor valor propio con la matriz S
	 DO k=1,n
	  DO l=1,n
	   suma1=suma1+av1(l)*ss(k,l)*av1(k) !Almaceno en memoria la variable suma para dividir los elementos...
	  END DO			     !... por su raiz cuadrada y que así queden normalizados
	 END DO

	suma2=0.d0
	 DO k=1,n
	  DO l=1,n
	   suma2=suma2+av2(l)*ss(k,l)*av2(k)
	  END DO
	 END DO

	suma3=0.d0
	 DO k=1,n
	  DO l=1,n
	   suma3=suma3+av3(l)*ss(k,l)*av3(k)
	  END DO
	 END DO

c NORMALIZACIÓN			

	DO k=1,n
	 av1(k)=av1(k)/(dsqrt(dabs(suma1)))
	 av2(k)=av2(k)/(dsqrt(dabs(suma2)))
	 av3(k)=av3(k)/(dsqrt(dabs(suma3)))
	END DO

	norma(1)=dsqrt(dabs(suma1))	!Guardo en memoria la norma de cada vector para utilizarla posteriormente...
	norma(2)=dsqrt(dabs(suma2))	!... al calcular las funciones de onda
	norma(3)=dsqrt(dabs(suma3))

c Compruebo normalización escribiendo en archivo, es correcta pues obtenemos 1

c	suma1=0.d0
c	 DO k=1,n
c	  DO l=1,n
c	   suma1=suma1+av1(l)*ss(k,l)*av1(k)
c	  END DO
c	 END DO
c
c	suma2=0.d0
c	 DO k=1,n
c	  DO l=1,n
c	   suma2=suma2+av2(l)*ss(k,l)*av2(k)
c	  END DO
c	 END DO
c
c	suma3=0.d0
c	 DO k=1,n
c	  DO l=1,n
c	   suma3=suma3+av3(l)*ss(k,l)*av3(k)
c	  END DO
c	 END DO
c
c	WRITE(*,*)suma1,suma2,suma3,j


	imin(1)=aa	!guardamos el indice para utilizarlo posteriormente
	imin(2)=bb
	imin(3)=cc


c CÁLCULO DE LAS FUNCIONES DE ONDA
c======================================================================

	DO i=1,nr		!nr=3. Consideramos las tres energías más bajas

	r=0.d0 
	

	DO l=0,99
	  r=((limsup-liminf)/99.d0)*l 	!Límite superior es 1,5 y el inferior 0

	    fi=0.d0
            DO k=1,n		!Phi=Nnl*(r^l)*Exp(-etan*(r^2))*Ĉnl



	      fi=fi+(Nnl(k,j))*(r**j)*		
     &        (dexp(-1.d0*etan(k)*(r**2.d0)))*
     &        ((av(k,imin(i)))/norma(i))
  
	    END DO
	
	    WRITE(36,'(2I4,2E24.12)')i,j,r,fi  !escribo las funciones en archivo con el formato requerido

c Con este case escribimos las funciones en diferentes archivos según su "n" y su "l" para representarlas en GNUPLOT

c          IF (j==0) THEN  
c		SELECT CASE (i)
c			CASE (1)                                       							   
c	    		WRITE(30,'(2I4,2E24.12)')i,j,r,fi 
c			CASE (2)                                       							   
c	    		WRITE(31,'(2I4,2E24.12)')i,j,r,fi  
c			CASE (3)                                       							   
c	    		WRITE(32,'(2I4,2E24.12)')i,j,r,fi  
c		END SELECT			                       
c	  ELSE
c		SELECT CASE (i)
c			CASE (1)                                       							   
c	    		WRITE(33,'(2I4,2E24.12)')i,j,r,fi  
c			CASE (2)                                       							   
c	    		WRITE(34,'(2I4,2E24.12)')i,j,r,fi 
c			CASE (3)                                       							   
c	    		WRITE(35,'(2I4,2E24.12)')i,j,r,fi   
c		END SELECT                        
c 	  END IF	
	END DO

	END DO
        

c DIFERENCIA DE MASA PARA nr=1, l=0
c=====================================================================


	IF (j==0) THEN	!Introducimos el bucle if por si deseamos calcular la dif de masas con otra l

	difmas=0.d0
	i=1

	DO k=1,n
	 DO l=1,n

	  difmas=difmas+av1(k)*			!Sumatorio de la diferencia de masas
     &    av1(l)*Nnl(k,j)*Nnl(l,j)
c     &    /(norma(i)**2.d0)			!Podemos normalizar aquí los vectores (otra opción)

	 END DO

	END DO

	 difmas=((8.d0*alfas*(hc**3.d0))/	!Expresión final
     &   (9.d0*(mb**2.d0)))*difmas

c Salida en pantalla de diferencia de masa n=1,l=0
	
	WRITE(*,*)"Separacion de masa entre diferentes proyecciones de spi
     &n para n=1,l=0:"
	WRITE(*,"(F15.10)")difmas
	WRITE(*,*)
	END IF



        END DO	!End do de l=0,1
	END PROGRAM

c Aquí finaliza el programa. Hemos incluido en el programa la subrutina que calcula la función gamma por comodidad.


c SUBRUTINA QUE CALCULA VALORES DE LA FUNCIÓN GAMMA	
c _________________________________________________

       SUBROUTINE GammaFun(X,GA)
C
C       ==================================================

C       Input :  x  --- Argument of Gamma(x)
C                       ( x is not equal to 0,-1,-2,...)
C       Output:  GA --- Gamma(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0,
     &          -0.420026350340952D-1,
     &          0.1665386113822915D0,
     &          -.421977345555443D-1,
     &          -.96219715278770D-2, 
     &           .72189432466630D-2,
     &          -.11651675918591D-2, 
     &          -.2152416741149D-3,
     &          .1280502823882D-3, 
     &          -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
