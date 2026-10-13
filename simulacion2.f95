module simulacion2
    implicit none
    
    ! Decaro las variables globales del sistema
    real*8 :: coordenadas(4), step_size, tiempo, tiempoFinal
    real*8 :: coordenadaEntera(4), mediaCoordenada(4), coordTierra(2), coordMarte(2), coordLuna(2), coordLT(2)
    real*8 :: error(4), errorMaximo, s(5), hMax, sMax, hMaxLimit
    integer*8 :: i, unit_resultado, j, k, z
    real*8 :: altura, r, vOrbCohete, distanciaCohete, alturaTierra
    !Valores de impulso
    real*8 :: tiempoImpulsos(1), impulso(1), direccionImpulso(1, 2), energia, energiaInicial, distanciaMarte(2), modDistanciaMarte
    integer :: boolImpulso(1)

    ! Aquí defino el valor de las constantes como parámetros
    real*8, parameter :: pi = 4.0d0 * atan(1.0d0)  ! Valor de pi
    real*8, parameter :: G = 6.67430e-11           ! Constante gravitacional en Nm^2/kg^2
    real*8, parameter :: Ms = 1.989e30             ! Masa del Sol en kg
    real*8, parameter :: Mt = 5.972e24             ! Masa de la Tierra en kg
    real*8, parameter :: Ml = 7.34767309e22        ! Masa de la Luna en kg
    real*8, parameter :: Mm = 6.4171e23            ! Masa de Marte en kg
    real*8, parameter :: m = 2000.0                ! Masa del objeto en kg
    real*8, parameter :: L = 1.496e11              ! Distancia en metros del Sol a la Tierra

    real*8, parameter :: Mtotal = Ms + Mt + Ml + Mm  ! Masa total del sistema

    real*8, parameter :: tCaract = sqrt(L**3 / (G * Mtotal))  ! Tiempo característico del sistema

    real*8, parameter :: MtCV = Mt / Mtotal  ! Masa de la Tierra en términos de la masa total
    real*8, parameter :: MlCV = Ml / Mtotal  ! Masa de la Luna en términos de la masa
    real*8, parameter :: MmCV = Mm / Mtotal  ! Masa de Marte en términos de la masa total
    real*8, parameter :: MsCV = Ms / Mtotal  ! Masa del Sol en términos de la masa total

    ! Radios reescalados
    real*8, parameter :: radioTierra = 6.371e6 / L ! Radio reescalado de la Tierra
    real*8, parameter :: radioMarte = 3.3895e6 / L ! Radio reescalado de Marte
    real*8, parameter :: radioLuna = 1.7374e6 / L  ! Radio reescalado de la Luna

    ! Velocidades angulares
    real*8, parameter :: wt = tCaract*2.0 * pi / 31557600.0  ! Velocidad angular de la Tierra alrededor del Sol
    real*8, parameter :: wm = tCaract*2.0 * pi / 59353000.0  ! Velocidad angular de Marte alrededor del Sol
    real*8, parameter :: wl = tCaract*2.0 * pi / 2360592.0   ! Velocidad angular de la Luna alrededor de la Tierra

    ! Distancias orbitales reescaladas
    real*8, parameter :: dOrbTierra = 1.0             ! Distancia orbital de la Tierra al Sol
    real*8, parameter :: dOrbMarte = 2.279e11 / L     ! Distancia orbital de Marte al Sol
    real*8, parameter :: dOrbLuna = 3.844e8 / L       ! Distancia orbital de la Luna a la Tierra

    ! Desfases iniciales
    real*8, parameter :: desfIniTierra = 0.0          ! Desfase inicial de la Tierra
    real*8, parameter :: desfIniMarte = 55.0d0*pi/180.0d0   ! Desfase inicial de Marte
    real*8, parameter :: desfIniLuna = 0.0            ! Desfase inicial de la Luna
    

    public :: ejecutarSimulacion
    contains

        function calcularSimulacion(t, v) result(res)
            real*8 :: res(6), velocidadRel(2), vOrbSol, vOrbTierra, wc
            real*8 :: t, v, theta, distMinima, distMarte, energia, distSol, distTierra, velocidad
            real*8 :: velocidadTierra(2), velocidadRelativa(2), magnitudVelocidadRelativa, tiempoMarte
            integer*8 :: numEntero, unit_resultado
            real*8 :: boolLlegada

            !unit_resultado = 10
            !open(unit=unit_resultado, file='coordenadas.txt', status='replace')

            !open(unit=unit_resultado, file='coordenadas.txt', status='replace')

            tiempoMarte = 0.0d0
            boolLlegada = 0.0
            z = 0
            tiempoImpulsos = (/t/)
            impulso = (/v /)
            do k = 1, size(tiempoImpulsos)
                boolImpulso(k) = 0 ! Inicializo el booleano de impulso
            end do
        
                !----------------------------------------------------------------------------------
            ! VALORES INICIALES DEL COHETE'
            altura = dOrbTierra*L + radioTierra*L + 250.0e3
            alturaTierra = altura - dOrbTierra*L
            r = altura/L
            !vOrbCohete = sqrt(G * Mt / altura) / L + wt * dOrbTierra * cos(wt * 0.0 + desfIniTierra)
            !vOrbCohete = sqrt(G*Ms /altura)*tCaract/L + sqrt(G*Mt/alturaTierra)*tCaract/(L)
            vOrbSol = sqrt(G * Ms / (dOrbTierra * L)) * tCaract / L

            vOrbTierra = sqrt(G*Mt/alturaTierra)*tCaract/(L)

            wc = vOrbTierra/(alturaTierra/L)

            coordTierra = calcularTierra(tiempoImpulsos(1))

            
            
            coordenadas(1) = coordTierra(1) + alturaTierra*cos(wc*tiempoImpulsos(1))/L   ! Posición inicial en el eje x
            coordenadas(2) = coordTierra(2) + alturaTierra*sin(wc*tiempoImpulsos(1))/L      ! Posición inicial en el eje y

            !coordenadas(3) = dOrbTierra * wt * sin(wt * tiempoImpulsos(1) + desfIniTierra) + &
            !    vOrbTierra * cos(wc * tiempoImpulsos(1)) ! Velocidad inicial en el eje y
            !coordenadas(4) = -dOrbTierra * wt * cos(wt * tiempoImpulsos(1) + desfIniTierra) + &
            !    vOrbTierra * sin(wc * tiempoImpulsos(1)) ! Velocidad inicial en el eje x

            coordenadas(3) = -dOrbTIerra*wt*sin(wt* tiempoImpulsos(1) + desfIniTierra) - &
                alturaTierra*wc*sin(wc * tiempoImpulsos(1))/L
            
            coordenadas(4) = dOrbTierra*wt*cos(wt*tiempoImpulsos(1) + desfIniTierra) + &
                            alturaTierra*wc*cos(wc * tiempoImpulsos(1))/L
        
            !print *, coordenadas(3)
            !print *, coordenadas(4)
            !  Defino el bucle principal
        
            tiempo = tiempoImpulsos(1) ! Tiempo inicial
            tiempoFinal = 365*24*3600.0d0/tCaract ! 12 meses en segundos
            step_size = 60.0d0/tCaract! El tiempo esta normalizado con la masa
            errorMaximo = step_size/100.0d0 ! Error máximo permitido
            j = 0
            !hMaxLimit = 3*60/m ! Límite máximo del paso de tiempo
        
            s(1) = 1.0d-8

            
            modDistanciaMarte = lMarte(coordenadas, tiempo)


            distMinima = modDistanciaMarte
        
            energiaInicial = calcularEnergia(coordenadas, tiempo)
            do while (tiempo < tiempoFinal)
                !En el bucle aplicaremos el método de la h adaptativa
            
            
                coordenadaEntera = rungeKutta(coordenadas, tiempo, step_size)
                mediaCoordenada = rungeKutta(coordenadas, tiempo, step_size / 2.0d0)
                mediaCoordenada = rungeKutta(mediaCoordenada, tiempo + step_size / 2.0d0, step_size / 2.0d0)
            
                error = abs(coordenadaEntera - mediaCoordenada)*16.0d0/15.0d0
            
                do i = 2, 5
                    s(i) = (error(i - 1) / errorMaximo)**0.2d0
                end do
            
                sMax = maxVal(s)
            
                hMax = step_size/sMax
            
                if (sMax > 2.d0) then
                    !print *, t, " ", v, " ", z

                    step_size = hMax
                    z = z + 1
                
                else if (sMax < 2.0d0) then
                
                    coordenadas = mediaCoordenada
                
                    ! Escribo las coordenadas en un fichero
                    coordTierra = calcularTierra(tiempo+step_size)
                    coordMarte = calcularMarte(tiempo+step_size)
                    coordLuna = calcularLuna(tiempo+step_size)
                    coordLT = calcularLunaTierra(tiempo+step_size)
                
                    !numEntero = 10**9
                    !if(mod(j, numEntero) == 0) then ! Cada 100 iteraciones
                    !    print *,tiempo*tCaract/(3600.0d0), &
                    !        coordenadas(1), coordenadas(2), &
                    !        coordTierra(1), coordTierra(2), &
                    !        coordMarte(1), coordMarte(2), &
                    !        coordLuna(1), coordLuna(2), &
                    !        coordLT(1), coordLT(2)
                    !end if

                    !if(mod(j, 1) == 0) then ! Cada 100 iteraciones
                    !    print *, tiempo*tCaract/(3600.0d0), &
                    !        coordenadas(1), coordenadas(2), &
                    !        coordTierra(1), coordTierra(2), &
                    !        coordMarte(1), coordMarte(2), &
                    !        coordLuna(1), coordLuna(2), &
                    !        coordLT(1), coordLT(2)
                    !end if
                !
                    do k = 1, size(tiempoImpulsos)
                        if (tiempo >= tiempoImpulsos(k) .AND. boolImpulso(k) == 0) then
                            ! Calcular la velocidad de la Tierra en el instante actual
        
                            velocidadRel = velRelTierra(tiempo, coordenadas)
                            velocidad = sqrt(velocidadRel(1)**2 + velocidadRel(2)**2)
                            direccionImpulso(1, 1) = velocidadRel(1) / velocidad
                            direccionImpulso(1, 2) = velocidadRel(2) / velocidad
        
                            ! Aplicar el impulso en la dirección calculada
                            coordenadas(3) = coordenadas(3) + impulso(k) * direccionImpulso(1, 1)
                            coordenadas(4) = coordenadas(4) + impulso(k) * direccionImpulso(1, 2)
                            boolImpulso(k) = 1
                            !print *, "Impulso aplicado en el tiempo:", tiempo*tCaract/(60.0d0), "min"
                        end if
                    enddo

                    modDistanciaMarte = lMarte(coordenadas, tiempo)

                    if(lMarte(coordenadas, tiempo) < radioMarte) then
                        !print *, 'Cohete ha llegado a Marte (colisión)'
                        boolLlegada = -1.0d0
                        exit
                    else if(lMarte(coordenadas, tiempo) < radioMarte + (1.0e6 / L)) then
                        !print *, 'Cohete está en órbita cerca de Marte'
                        boolLlegada = 1.5
                        exit
                    end if

                    if(modDistanciaMarte < distMinima) then
                        distMinima = modDistanciaMarte
                        distSol = sqrt(coordenadas(1)**2 + coordenadas(2)**2)
                        distTierra = sqrt((coordenadas(1) - coordTierra(1))**2 + (coordenadas(2) - coordTierra(2))**2)
                    endif

                    



                    if(lLuna(coordenadas, tiempo) < radioLuna) then
                        !print *, 'Cohete ha llegado a la Luna'
                        exit
                    end if

                    
                    if (lTierra(coordenadas, tiempo) < radioTierra) then
                        !print *, 'Cohete ha llegado a la Tierra'
                        exit
                    end if


                    if(tiempo >= 1.2*30*24*3600/tCaract .AND. lTierra(coordenadas, tiempo) < 1.3*dOrbLuna) then
                        exit
                    end if

                    if(sqrt(coordenadas(1)**2 + coordenadas(2)**2) < 0.9*dOrbTierra) then
                        !print *, 'Cohete ha llegado al Sol'
                        exit
                    end if
                
                    

                    tiempo = tiempo + step_size
                    z = 0
                    j = j + 1

                    !if (mod(j, 10**6) == 0) then
                    !    print *, tiempo*tCaract/3600.0d0, " ", t, " " , v
                    !end if

                
                    if (step_size < hMax) then
                        step_size = step_size * 2.0d0
                        !step_size = min(2 * step_size, hMaxLimit)
                    end if
                
                end if


            end do
            !close(unit_resultado)

            res(1) = distMinima
            res(2) = v
            res(3) = tiempo
            res(4) = distSol
            res(5) = distTierra
            res(6) = boolLlegada

            !close(unit_resultado)

        end function calcularSimulacion

        function representarSimulacion(t, v, numeroCohete, generacion) result(res)
            real*8 :: res(5), velocidadRel(2), vOrbSol, vOrbTierra, wc
            real*8 :: t, v, theta, distMinima, distMarte, energia, distSol, distTierra, velocidad
            real*8 :: velocidadTierra(2), velocidadRelativa(2), magnitudVelocidadRelativa
            integer*8 :: numEntero, unit_resultado
            integer :: numeroCohete, generacion
            character(len=50) :: coheteStr, generacionStr, fileName

            unit_resultado = 10
            !open(unit=unit_resultado, file='coordenadas.txt', status='replace')

            

            write(coheteStr, '(I0)') numeroCohete
            write(generacionStr, '(I0)') generacion

            ! Concatenar las cadenas para formar el nombre del archivo
            fileName = 'coordenadas/coordenadas_N' // trim(adjustl(coheteStr)) // '_gen' // trim(adjustl(generacionStr)) // '.txt'

            ! Abrir el archivo con el nombre generado
            open(unit=unit_resultado, file=fileName, status='replace')

            z = 0
            tiempoImpulsos = (/t/)
            impulso = (/v /)
            do k = 1, size(tiempoImpulsos)
                boolImpulso(k) = 0 ! Inicializo el booleano de impulso
            end do
        
                !----------------------------------------------------------------------------------
            ! VALORES INICIALES DEL COHETE'
            altura = dOrbTierra*L + radioTierra*L + 250.0e3
            alturaTierra = altura - dOrbTierra*L
            r = altura/L
            !vOrbCohete = sqrt(G * Mt / altura) / L + wt * dOrbTierra * cos(wt * 0.0 + desfIniTierra)
            !vOrbCohete = sqrt(G*Ms /altura)*tCaract/L + sqrt(G*Mt/alturaTierra)*tCaract/(L)
            vOrbSol = sqrt(G * Ms / (dOrbTierra * L)) * tCaract / L

            vOrbTierra = sqrt(G*Mt/alturaTierra)*tCaract/(L)

            wc = vOrbTierra/(alturaTierra/L)

            coordTierra = calcularTierra(tiempoImpulsos(1))

            
            
            coordenadas(1) = coordTierra(1) + alturaTierra*cos(wc*tiempoImpulsos(1))/L   ! Posición inicial en el eje x
            coordenadas(2) = coordTierra(2) + alturaTierra*sin(wc*tiempoImpulsos(1))/L      ! Posición inicial en el eje y

            !coordenadas(3) = dOrbTierra * wt * sin(wt * tiempoImpulsos(1) + desfIniTierra) + &
            !    vOrbTierra * cos(wc * tiempoImpulsos(1)) ! Velocidad inicial en el eje y
            !coordenadas(4) = -dOrbTierra * wt * cos(wt * tiempoImpulsos(1) + desfIniTierra) + &
            !    vOrbTierra * sin(wc * tiempoImpulsos(1)) ! Velocidad inicial en el eje x

            coordenadas(3) = -dOrbTIerra*wt*sin(wt* tiempoImpulsos(1) + desfIniTierra) - &
                alturaTierra*wc*sin(wc * tiempoImpulsos(1))/L
            
            coordenadas(4) = dOrbTierra*wt*cos(wt*tiempoImpulsos(1) + desfIniTierra) + &
                            alturaTierra*wc*cos(wc * tiempoImpulsos(1))/L
        
            !print *, coordenadas(3)
            !print *, coordenadas(4)
            !  Defino el bucle principal
        
            tiempo = tiempoImpulsos(1) ! Tiempo inicial
            tiempoFinal = 365*24*3600.0d0/tCaract ! 12 meses en segundos
            step_size = 60.0d0/tCaract! El tiempo esta normalizado con la masa
            errorMaximo = step_size/100.0d0 ! Error máximo permitido
            j = 0
            !hMaxLimit = 3*60/m ! Límite máximo del paso de tiempo
        
            s(1) = 1.0d-8

            
            modDistanciaMarte = lMarte(coordenadas, tiempo)


            distMinima = modDistanciaMarte
        
            energiaInicial = calcularEnergia(coordenadas, tiempo)
            do while (tiempo < tiempoFinal)
                !En el bucle aplicaremos el método de la h adaptativa
            
            
                coordenadaEntera = rungeKutta(coordenadas, tiempo, step_size)
                mediaCoordenada = rungeKutta(coordenadas, tiempo, step_size / 2.0d0)
                mediaCoordenada = rungeKutta(mediaCoordenada, tiempo + step_size / 2.0d0, step_size / 2.0d0)
            
                error = abs(coordenadaEntera - mediaCoordenada)*16.0d0/15.0d0
            
                do i = 2, 5
                    s(i) = (error(i - 1) / errorMaximo)**0.2d0
                end do
            
                sMax = maxVal(s)
            
                hMax = step_size/sMax
            
                if (sMax > 2.d0) then
                    !print *, t, " ", v, " ", z

                    step_size = hMax
                    z = z + 1
                
                else if (sMax < 2.0d0) then
                
                    coordenadas = mediaCoordenada
                
                    ! Escribo las coordenadas en un fichero
                    coordTierra = calcularTierra(tiempo+step_size)
                    coordMarte = calcularMarte(tiempo+step_size)
                    coordLuna = calcularLuna(tiempo+step_size)
                    coordLT = calcularLunaTierra(tiempo+step_size)
                
                    !numEntero = 10**9
                    !if(mod(j, numEntero) == 0) then ! Cada 100 iteraciones
                    !    print *,tiempo*tCaract/(3600.0d0), &
                    !        coordenadas(1), coordenadas(2), &
                    !        coordTierra(1), coordTierra(2), &
                    !        coordMarte(1), coordMarte(2), &
                    !        coordLuna(1), coordLuna(2), &
                    !        coordLT(1), coordLT(2)
                    !end if

                    if(mod(j, 1) == 0) then ! Cada 100 iteraciones
                        write(unit_resultado,*), tiempo*tCaract/(3600.0d0), &
                            coordenadas(1), coordenadas(2), &
                            coordTierra(1), coordTierra(2), &
                            coordMarte(1), coordMarte(2), &
                            coordLuna(1), coordLuna(2), &
                            coordLT(1), coordLT(2)
                    end if
                !
                    do k = 1, size(tiempoImpulsos)
                        if (tiempo >= tiempoImpulsos(k) .AND. boolImpulso(k) == 0) then
                            ! Calcular la velocidad de la Tierra en el instante actual
        
                            velocidadRel = velRelTierra(tiempo, coordenadas)
                            velocidad = sqrt(velocidadRel(1)**2 + velocidadRel(2)**2)
                            direccionImpulso(1, 1) = velocidadRel(1) / velocidad
                            direccionImpulso(1, 2) = velocidadRel(2) / velocidad
        
                            ! Aplicar el impulso en la dirección calculada
                            coordenadas(3) = coordenadas(3) + impulso(k) * direccionImpulso(1, 1)
                            coordenadas(4) = coordenadas(4) + impulso(k) * direccionImpulso(1, 2)
                            boolImpulso(k) = 1
                            !print *, "Impulso aplicado en el tiempo:", tiempo*tCaract/(60.0d0), "min"
                        end if
                    enddo

                    modDistanciaMarte = lMarte(coordenadas, tiempo)

                    if(modDistanciaMarte < distMinima) then
                        distMinima = modDistanciaMarte
                        distSol = sqrt(coordenadas(1)**2 + coordenadas(2)**2)
                        distTierra = sqrt((coordenadas(1) - coordTierra(1))**2 + (coordenadas(2) - coordTierra(2))**2)
                    endif

                    if(modDistanciaMarte < 200.0e3/L) then
                        !print *, 'Cohete ha llegado a marte'
                        exit
                    end if

                    
                    if (lTierra(coordenadas, tiempo) < radioTierra) then
                        !print *, 'Cohete ha llegado a la Tierra'
                        exit
                    end if


                    if(tiempo >= 1.2*30*24*3600/tCaract .AND. lTierra(coordenadas, tiempo) < 1.3*dOrbLuna) then
                        exit
                    end if

                    if(sqrt(coordenadas(1)**2 + coordenadas(2)**2) < 0.9*dOrbTierra) then
                        !print *, 'Cohete ha llegado al Sol'
                        exit
                    end if
                
                    

                    tiempo = tiempo + step_size
                    z = 0
                    j = j + 1

                    !if (mod(j, 10**6) == 0) then
                    !    print *, tiempo*tCaract/3600.0d0, " ", t, " " , v
                    !end if

                
                    if (step_size < hMax) then
                        step_size = step_size * 2.0d0
                        !step_size = min(2 * step_size, hMaxLimit)
                    end if
                
                end if


            end do
            close(unit_resultado)

            res(1) = distMinima
            res(2) = v
            res(3) = tiempo
            res(4) = distSol
            res(5) = distTierra

            !close(unit_resultado)

        end function representarSimulacion

        subroutine ejecutarSimulacion(t, v, resultado)
            real*8 :: t, v, resultado(6)

            ! Llamar a la función calcularSimulacion
            resultado = calcularSimulacion(t, v)
        end subroutine ejecutarSimulacion

        function representarOrbita(t, v, numeroCohete, generacion, dMin) result(res)
            real*8 :: res(7), velocidadRel(2), vOrbSol, vOrbTierra, wc, dMin, tiempoMarte, tiempoLlegada
            real*8 :: t, v, theta, distMinima, distMarte, energia, distSol, distTierra, velocidad
            real*8 :: velocidadTierra(2), velocidadRelativa(2), magnitudVelocidadRelativa, velocidadImpulso(2)
            integer*8 :: numEntero, unit_resultado
            integer :: numeroCohete, generacion, f
            character(len=50) :: coheteStr, generacionStr, fileName
            real*8 :: velocidadRelativaMarte(2), energiaCineticaRelativa, energiaPotencialRelativa, energiaRelativaMarte

            unit_resultado = 10
            !open(unit=unit_resultado, file='coordenadas.txt', status='replace')

            
            f = 0
            write(coheteStr, '(I0)') numeroCohete
            write(generacionStr, '(I0)') generacion

            ! Concatenar las cadenas para formar el nombre del archivo
            fileName = 'orbitaCoordenada/coordenadas' // trim(adjustl(coheteStr)) // '_gen' // trim(adjustl(generacionStr)) // '.txt'

            ! Abrir el archivo con el nombre generado
            open(unit=unit_resultado, file=fileName, status='replace')

            z = 0
            tiempoImpulsos = (/t/)
            impulso = (/v /)
            do k = 1, size(tiempoImpulsos)
                boolImpulso(k) = 0 ! Inicializo el booleano de impulso
            end do
        
                !----------------------------------------------------------------------------------
            ! VALORES INICIALES DEL COHETE'
            altura = dOrbTierra*L + radioTierra*L + 250.0e3
            alturaTierra = altura - dOrbTierra*L
            r = altura/L
            !vOrbCohete = sqrt(G * Mt / altura) / L + wt * dOrbTierra * cos(wt * 0.0 + desfIniTierra)
            !vOrbCohete = sqrt(G*Ms /altura)*tCaract/L + sqrt(G*Mt/alturaTierra)*tCaract/(L)
            vOrbSol = sqrt(G * Ms / (dOrbTierra * L)) * tCaract / L

            vOrbTierra = sqrt(G*Mt/alturaTierra)*tCaract/(L)

            wc = vOrbTierra/(alturaTierra/L)

            coordTierra = calcularTierra(tiempoImpulsos(1))

            
            
            coordenadas(1) = coordTierra(1) + alturaTierra*cos(wc*tiempoImpulsos(1))/L   ! Posición inicial en el eje x
            coordenadas(2) = coordTierra(2) + alturaTierra*sin(wc*tiempoImpulsos(1))/L      ! Posición inicial en el eje y

            !coordenadas(3) = dOrbTierra * wt * sin(wt * tiempoImpulsos(1) + desfIniTierra) + &
            !    vOrbTierra * cos(wc * tiempoImpulsos(1)) ! Velocidad inicial en el eje y
            !coordenadas(4) = -dOrbTierra * wt * cos(wt * tiempoImpulsos(1) + desfIniTierra) + &
            !    vOrbTierra * sin(wc * tiempoImpulsos(1)) ! Velocidad inicial en el eje x

            coordenadas(3) = -dOrbTIerra*wt*sin(wt* tiempoImpulsos(1) + desfIniTierra) - &
                alturaTierra*wc*sin(wc * tiempoImpulsos(1))/L
            
            coordenadas(4) = dOrbTierra*wt*cos(wt*tiempoImpulsos(1) + desfIniTierra) + &
                            alturaTierra*wc*cos(wc * tiempoImpulsos(1))/L
        
            !print *, coordenadas(3)
            !print *, coordenadas(4)
            !  Defino el bucle principal
        
            tiempo = tiempoImpulsos(1) ! Tiempo inicial
            tiempoFinal = 365*0.7*24*3600.0d0/tCaract ! 12 meses en segundos
            step_size = 60.0d0/tCaract! El tiempo esta normalizado con la masa
            errorMaximo = step_size/100.0d0 ! Error máximo permitido
            j = 0
            !hMaxLimit = 3*60/m ! Límite máximo del paso de tiempo
        
            s(1) = 1.0d-8

            
            modDistanciaMarte = lMarte(coordenadas, tiempo)


            distMinima = modDistanciaMarte
        
            energiaInicial = calcularEnergia(coordenadas, tiempo)
            z=0
            do while (tiempo < tiempoFinal)
                !En el bucle aplicaremos el método de la h adaptativa
            
            
                coordenadaEntera = rungeKutta(coordenadas, tiempo, step_size)
                mediaCoordenada = rungeKutta(coordenadas, tiempo, step_size / 2.0d0)
                mediaCoordenada = rungeKutta(mediaCoordenada, tiempo + step_size / 2.0d0, step_size / 2.0d0)
            
                error = abs(coordenadaEntera - mediaCoordenada)*16.0d0/15.0d0
            
                do i = 2, 5
                    s(i) = (error(i - 1) / errorMaximo)**0.2d0
                end do
            
                sMax = maxVal(s)
            
                hMax = step_size/sMax
            
                if (sMax > 2.d0) then
                    !print *, t, " ", v, " ", z

                    step_size = hMax
                    z = z + 1
                
                else if (sMax < 2.0d0) then
                
                    coordenadas = mediaCoordenada
                
                    ! Escribo las coordenadas en un fichero
                    coordTierra = calcularTierra(tiempo+step_size)
                    coordMarte = calcularMarte(tiempo+step_size)
                    coordLuna = calcularLuna(tiempo+step_size)
                    coordLT = calcularLunaTierra(tiempo+step_size)
                
                    !numEntero = 10**9
                    !if(mod(j, numEntero) == 0) then ! Cada 100 iteraciones
                    !    print *,tiempo*tCaract/(3600.0d0), &
                    !        coordenadas(1), coordenadas(2), &
                    !        coordTierra(1), coordTierra(2), &
                    !        coordMarte(1), coordMarte(2), &
                    !        coordLuna(1), coordLuna(2), &
                    !        coordLT(1), coordLT(2)
                    !end if

                    if(mod(j, 1) == 0) then ! Cada 100 iteraciones
                        write(unit_resultado,*), tiempo*tCaract/(3600.0d0), &
                            coordenadas(1), coordenadas(2), &
                            coordTierra(1), coordTierra(2), &
                            coordMarte(1), coordMarte(2), &
                            coordLuna(1), coordLuna(2), &
                            coordLT(1), coordLT(2)
                    end if
                !
                    do k = 1, size(tiempoImpulsos)
                        if (tiempo >= tiempoImpulsos(k) .AND. boolImpulso(k) == 0) then
                            ! Calcular la velocidad de la Tierra en el instante actual
        
                            velocidadRel = velRelTierra(tiempo, coordenadas)
                            velocidad = sqrt(velocidadRel(1)**2 + velocidadRel(2)**2)
                            direccionImpulso(1, 1) = velocidadRel(1) / velocidad
                            direccionImpulso(1, 2) = velocidadRel(2) / velocidad
        
                            ! Aplicar el impulso en la dirección calculada
                            coordenadas(3) = coordenadas(3) + impulso(k) * direccionImpulso(1, 1)
                            coordenadas(4) = coordenadas(4) + impulso(k) * direccionImpulso(1, 2)
                            boolImpulso(k) = 1
                            !print *, "Impulso aplicado en el tiempo:", tiempo*tCaract/(60.0d0), "min"
                        end if
                    enddo

                    modDistanciaMarte = lMarte(coordenadas, tiempo)

                    if(modDistanciaMarte < distMinima) then
                        distMinima = modDistanciaMarte
                        distSol = sqrt(coordenadas(1)**2 + coordenadas(2)**2)
                        distTierra = sqrt((coordenadas(1) - coordTierra(1))**2 + (coordenadas(2) - coordTierra(2))**2)
                    endif

                    if(modDistanciaMarte < dMin*1.2 .AND. f == 0) then
                        ! Vamos a calcular la orbita para orbitar marte
                        velocidadImpulso = calcularVelocidadOrbita(coordenadas, tiempo)

                        coordenadas(3) = velocidadImpulso(1) + coordenadas(3)
                        coordenadas(4) = velocidadImpulso(2) + coordenadas(4)

                        tiempoFinal = tiempo + 24.0d0 * 3600.0d0 / tCaract

                        f = 1 ! Marcamos que el impulso ya se realizó
                    end if
                    
                    if (lTierra(coordenadas, tiempo) < radioTierra) then
                        !print *, 'Cohete ha llegado a la Tierra'
                        exit
                    end if


                    if(tiempo >= 1.2*30*24*3600/tCaract .AND. lTierra(coordenadas, tiempo) < 1.3*dOrbLuna) then
                        exit
                    end if

                    if(sqrt(coordenadas(1)**2 + coordenadas(2)**2) < 0.9*dOrbTierra) then
                        !print *, 'Cohete ha llegado al Sol'
                        exit
                    end if
                
                    

                    tiempo = tiempo + step_size
                    z = 0
                    j = j + 1

                    !if (mod(j, 10**6) == 0) then
                    !    print *, tiempo*tCaract/3600.0d0, " ", t, " " , v
                    !end if

                
                    if (step_size < hMax) then
                        step_size = step_size * 2.0d0
                        !step_size = min(2 * step_size, hMaxLimit)
                    end if
                
                end if


            end do
            close(unit_resultado)

            ! Calcular la energía cinética relativa a Marte
            velocidadRelativaMarte = velRelMarte(tiempo, coordenadas)
            energiaCineticaRelativa = 0.5d0 * (velocidadRelativaMarte(1)**2 + velocidadRelativaMarte(2)**2)

            ! Calcular la energía potencial gravitatoria relativa a Marte
            energiaPotencialRelativa = -MmCV / (modDistanciaMarte)

            ! Calcular la energía total relativa a Marte
            energiaRelativaMarte = energiaCineticaRelativa + energiaPotencialRelativa

            

            res(1) = distMinima
            res(2) = v
            res(3) = tiempo
            res(4) = distSol
            res(5) = distTierra
            res(6) = impulso(1)
            res(7) = energiaRelativaMarte

            !close(unit_resultado)

        end function representarOrbita

        function calcularVelocidadOrbita(coordenadas, t) result(res)
            real*8 :: coordenadas(4), t, res(2)
            real*8 :: coordMarte(2), r(2), distMarte, tangente(2)
            real*8 :: velRelMarte(2), momentAngMarte, signMarte
            real*8 :: velocidadOrbMarte, velRelDeseada(2), velMarte(2), velDeseadaInercial(2)
        
            ! Coordenadas de Marte
            coordMarte = calcularMarte(t)
        
            ! Posición relativa respecto a Marte
            r = coordenadas(1:2) - coordMarte
            distMarte = sqrt(r(1)**2 + r(2)**2)
        
            ! Velocidad relativa al centro de Marte
            velMarte(1) = -dOrbMarte * sin(wm * t + desfIniMarte) * wm
            velMarte(2) =  dOrbMarte * cos(wm * t + desfIniMarte) * wm
            velRelMarte = coordenadas(3:4) - velMarte
        
            ! Momento angular respecto a Marte (r × v)
            momentAngMarte = r(1) * velRelMarte(2) - r(2) * velRelMarte(1)
            signMarte = momentAngMarte / abs(momentAngMarte)
        
            ! Vector normal (radial) y signo
            tangente(1) = -r(2) / distMarte
            tangente(2) = r(1) / distMarte
            tangente = tangente * signMarte
        
            ! Velocidad orbital deseada
            velocidadOrbMarte = sqrt(G*Mm/(distMarte*L)) * tCaract / L
            velRelDeseada = velocidadOrbMarte * tangente
        
            ! Velocidad deseada en el sistema inercial
            velDeseadaInercial = velMarte + velRelDeseada
        
            ! Impulso necesario
            res = velDeseadaInercial - coordenadas(3:4)
        end function calcularVelocidadOrbita
        
        


        ! Funciones para localizar los planetas
        function calcularTierra(t) result(res)
            real*8 :: t, res(2)

            res(1) = dOrbTierra * cos(wt * t + desfIniTierra)
            res(2) = dOrbTierra * sin(wt * t + desfIniTierra)
        end function calcularTierra


        function calcularMarte(t) result(res)
            real*8 :: t, res(2)

            res(1) = dOrbMarte * cos(wm * t + desfIniMarte)
            res(2) = dOrbMarte * sin(wm * t + desfIniMarte)
        end function calcularMarte

        function calcularLunaTierra(t) result(res)
            real*8 :: t, res(2)

            res(1) = dOrbLuna * cos(wl * t + desfIniLuna)
            res(2) = dOrbLuna * sin(wl * t + desfIniLuna)
        end function calcularLunaTierra

        
        function calcularLuna(t) result(res)
            real*8 :: t, res(2), relativo(2), tierra(2)

            tierra = calcularTierra(t)
            relativo = calcularLunaTierra(t)

            res = tierra + relativo
        end function calcularLuna

        ! Funciones del hamiltoniano, cada una es la derivada

        function funcionX(x, y, px, py, t) result(res)
            real*8 :: x, y, px, py, t, res

            res = px
            
        end function funcionX

        function funcionY(x, y, px, py, t) result(res)
            real*8 :: x, y, px, py, t, res

            res = py
            
        end function funcionY

        function funcionPx(x, y, px, py, t) result(res)
            real*8 :: x, y, px, py, t, res
            real*8 :: dTierra(2), dLuna(2), dMarte(2)

            dTierra = calcularTierra(t)
            dLuna = calcularLuna(t)
            dMarte = calcularMarte(t)


            res = MlCV * (-x + dLuna(1)) / ((x - dLuna(1))**2 + (y - dLuna(2))**2)**1.5 &
                + MmCV * (-x + dMarte(1)) / ((x - dMarte(1))**2 + (y - dMarte(2))**2)**1.5 &
                - MsCV * x / (x**2 + y**2)**1.5 &
                + MtCV * (-x + dTierra(1)) / ((x - dTierra(1))**2 + (y - dTierra(2))**2)**1.5
        end function funcionPx

        function funcionPy(x, y, px, py, t) result(res)
            real*8 :: x, y, px, py, t, res
            real*8 :: dTierra(2), dLuna(2), dMarte(2)

            dTierra = calcularTierra(t)
            dLuna = calcularLuna(t)
            dMarte = calcularMarte(t)

            res = MlCV * (-y + dLuna(2)) / ((x - dLuna(1))**2 + (y - dLuna(2))**2)**1.5 &
                + MmCV * (-y + dMarte(2)) / ((x - dMarte(1))**2 + (y - dMarte(2))**2)**1.5 &
                - MsCV * y / (x**2 + y**2)**1.5 &
                + MtCV * (-y + dTierra(2)) / ((x - dTierra(1))**2 + (y - dTierra(2))**2)**1.5
        end function funcionPy


        ! Función para calcular runge-Kutta
        function rungeKutta(coordenadas, t, h) result(res)
            real*8 :: coordenadas(4), t, res(4), h
            real*8 :: ncoordenadas(4)
            real*8 :: k1(4), k2(4), k3(4), k4(4)
            
            ! k1
            k1(1) = h * funcionX(coordenadas(1), coordenadas(2), coordenadas(3), coordenadas(4), t)
            k1(2) = h * funcionY(coordenadas(1), coordenadas(2), coordenadas(3), coordenadas(4), t)
            k1(3) = h * funcionPx(coordenadas(1), coordenadas(2), coordenadas(3), coordenadas(4), t)
            k1(4) = h * funcionPy(coordenadas(1), coordenadas(2), coordenadas(3), coordenadas(4), t)

            ! k2
            ncoordenadas = coordenadas + 0.5 * k1
            k2(1) = h * funcionX(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)
            k2(2) = h * funcionY(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)
            k2(3) = h * funcionPx(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)
            k2(4) = h * funcionPy(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)

            ! k3
            ncoordenadas = coordenadas + 0.5 * k2
            k3(1) = h * funcionX(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)
            k3(2) = h * funcionY(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)
            k3(3) = h * funcionPx(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)
            k3(4) = h * funcionPy(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + 0.5 * h)

            ! k4
            ncoordenadas = coordenadas + k3
            k4(1) = h * funcionX(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + h)
            k4(2) = h * funcionY(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + h)
            k4(3) = h * funcionPx(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + h)
            k4(4) = h * funcionPy(ncoordenadas(1), ncoordenadas(2), ncoordenadas(3), ncoordenadas(4), t + h)

            ! Calculamos el resultado
            res = coordenadas + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

        end function rungeKutta

        function calcularFlujo(t, coordenadas) result(flujo)
            real*8 :: t, flujo
            real*8 :: coordenadas(2), coordTierra(2), coordMarte(2), coordLuna(2)
            real*8 :: aporteT, aporteM, aporteL
            real*8 :: derivadaT(2), derivadaM(2), derivadaL(2)
            real*8 :: distanciaTierra(2), distanciaMarte(2), distanciaLuna(2)
            real*8 :: deltaT, deltaM, deltaL

            coordTierra = calcularTierra(t)
            coordMarte = calcularMarte(t)
            coordLuna = calcularLuna(t)

            ! Derivadas de las posiciones orbitales
            derivadaT(1) = (-1.0d0)*dOrbTierra*sin(wt*t)*wt
            derivadaT(2) = dOrbTierra*cos(wt*t)*wt

            derivadaM(1) = (-1.0d0)*dOrbMarte*sin(wm*t)*wm
            derivadaM(2) = dOrbMarte*cos(wm*t)*wm

            derivadaL(1) = (-1.0d0)*dOrbLuna*sin(wl*t)*wl
            derivadaL(2) = dOrbLuna*cos(wl*t)*wl

            ! Calculamos las distancias
            distanciaTierra(1) = coordenadas(1) - coordTierra(1)
            distanciaTierra(2) = coordenadas(2) - coordTierra(2)

            distanciaMarte(1) = coordenadas(1) - coordMarte(1)
            distanciaMarte(2) = coordenadas(2) - coordMarte(2)

            distanciaLuna(1) = coordenadas(1) - coordLuna(1)
            distanciaLuna(2) = coordenadas(2) - coordLuna(2)

            ! Calculamos el aporte de cada planeta
            deltaT = MtCV
            deltaM = MmCV
            deltaL = MlCV

            aporteT = dot_product(derivadaT, distanciaTierra) * deltaT / (sqrt(distanciaTierra(1)**2 + distanciaTierra(2)**2))**3
            aporteM = dot_product(derivadaM, distanciaMarte) * deltaM / (sqrt(distanciaMarte(1)**2 + distanciaMarte(2)**2))**3
            aporteL = dot_product(derivadaL, distanciaLuna) * deltaL / (sqrt(distanciaLuna(1)**2 + distanciaLuna(2)**2))**3

            ! Calculamos el flujo total
            flujo = aporteT + aporteM + aporteL

        end function calcularFlujo

        function calcularEnergia(coordenadas, t) result(energia)
            real*8 :: coordenadas(4), t, energia
            real*8 :: coordTierra(2), coordMarte(2), coordLuna(2)
            real*8 :: flujo

            ! Calculamos las posiciones de los cuerpos celestes
            coordTierra = calcularTierra(t)
            coordMarte = calcularMarte(t)
            coordLuna = calcularLuna(t)

            ! Calculamos la energía potencial gravitacional y cinética
            energia = 0.5d0 * (coordenadas(3)**2 + coordenadas(4)**2) &
                - MsCV / sqrt(coordenadas(1)**2 + coordenadas(2)**2) &
                - MtCV / sqrt((coordenadas(1) - coordTierra(1))**2 + (coordenadas(2) - coordTierra(2))**2) &
                - MmCV / sqrt((coordenadas(1) - coordMarte(1))**2 + (coordenadas(2) - coordMarte(2))**2) &
                - MlCV / sqrt((coordenadas(1) - coordLuna(1))**2 + (coordenadas(2) - coordLuna(2))**2)

            ! Calculamos el flujo
            flujo = calcularFlujo(t, coordenadas)

            ! Ajustamos la energía con el flujo
            energia = energia - flujo
        end function calcularEnergia

        function lMarte(coordenadas, t) result(res)
            real*8 :: coordenadas(4), t, vDistanciaMarte(2), coordMarte(2), res

            coordMarte = calcularMarte(t)

            vDistanciaMarte(1) = coordenadas(1) - coordMarte(1)
            vDistanciaMarte(2) = coordenadas(2) - coordMarte(2)

            res = sqrt(vDistanciaMarte(1)**2 + vDistanciaMarte(2)**2)

        end function lMarte

        function lTierra(coordenadas, t) result(res)
            real*8 :: coordenadas(4), t, vDistanciaTierra(2), coordTierra(2), res

            coordTierra = calcularTierra(t)

            vDistanciaTierra(1) = coordenadas(1) - coordTierra(1)
            vDistanciaTierra(2) = coordenadas(2) - coordTierra(2)

            res = sqrt(vDistanciaTierra(1)**2 + vDistanciaTierra(2)**2)

        end function lTierra

        function lLuna(coordenadas, t) result(res)
            real*8 :: coordenadas(4), t, vDistanciaLuna(2), coordLuna(2), res

            coordLuna = calcularLuna(t)

            vDistanciaLuna(1) = coordenadas(1) - coordLuna(1)
            vDistanciaLuna(2) = coordenadas(2) - coordLuna(2)

            res = sqrt(vDistanciaLuna(1)**2 + vDistanciaLuna(2)**2)

        end function lLuna

        function velRelTierra(t, coordenadas)
            real*8 :: t, calcVelocidadTierra(2), velRelTierra(2), coordenadas(4)

            calcVelocidadTierra(1) = -dOrbTierra * sin(wt * t + desfIniTierra) * wt
            calcVelocidadTierra(2) = dOrbTierra * cos(wt * t + desfIniTierra) * wt

            ! Calcular la velocidad relativa del cohete respecto a la Tierra
            velRelTierra(1) = coordenadas(3) - calcVelocidadTierra(1)
            velRelTierra(2) = coordenadas(4) - calcVelocidadTierra(2)

        end function velRelTierra

        function velRelMarte(t, coordenadas) 
            real*8 :: t, velocidadMarte(2), velRelMarte(2), coordenadas(4)

            velocidadMarte(1) = -dOrbMarte * sin(wm * t + desfIniMarte) * wm
            velocidadMarte(2) = dOrbMarte * cos(wm * t + desfIniMarte) * wm

            ! Calcular la velocidad relativa del cohete respecto a Marte
            velRelMarte(1) = coordenadas(3) - velocidadMarte(1)
            velRelMarte(2) = coordenadas(4) - velocidadMarte(2)

        end function velRelMarte


    
end module simulacion2