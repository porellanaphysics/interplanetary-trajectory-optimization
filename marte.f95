program marte
    implicit none



    ! Decaro las variables globales del sistema
    real*8 :: m, G, Ms, Mt, Ml, Mm, wt, wm, wl, L
    real*8 :: MtCV, MlCV, MmCV, MsCV ! Masas tras el cambio de variable
    real*8 :: dOrbTierra, dOrbMarte, dOrbLuna, desfIniTierra, desfIniMarte, desfIniLuna
    real*8 :: tCaract, Mtotal
    real*8 :: radioTierra, radioMarte, radioLuna
    real*8 :: coordenadas(4), step_size, tiempo, tiempoFinal
    real*8 :: coordenadaEntera(4), mediaCoordenada(4), coordTierra(2), coordMarte(2), coordLuna(2), coordLT(2)
    real*8 :: velocidadTierra(2), velocidadRelativa(2), magnitudVelocidadRelativa(2), velocidad
    real*8 :: error(4), errorMaximo, s(5), hMax, sMax, hMaxLimit
    integer*8 :: i, unit_resultado, j, k
    real*8 :: pi, altura, r, vOrbCohete, distanciaCohete, alturaTierra
    !Valores de impulso
    real*8 :: tiempoImpulsos(1), impulso(1), direccionImpulso(1, 2), energia, energiaInicial, theta
    real*8 :: velocidadRel(2)
    integer :: boolImpulso(1)
    real*8 :: semiEjeMayor, vOrbitaMarte, vTransferencia
    real*8 :: energiaOrbital, energiaTransferencia, velDeseadaInercial(2), velRelDeseada(2)
    real*8 :: flujo, aporteT, aporteM, aporteL


    !Aquí defino el valor de las constantes
    ! Inicializamos las constantes
    pi = 4.0d0 * atan(1.0d0)  ! Valor de pi
    G = 6.67430e-11           ! Constante gravitacional en Nm^2/kg^2
    Ms = 1.989e30             ! Masa del Sol en kg
    Mt = 5.972e24             ! Masa de la Tierra en kg
    Ml = 7.34767309e22        ! Masa de la Luna en kg
    Mm = 6.4171e23            ! Masa de Marte en kg
    m = 2000.0                ! Masa del objeto en kg
    L = 1.496e11              ! Distancia en metros del Sol a la Tierra

    mTotal = Ms + Mt + Ml + Mm  ! Masa total del sistema

    tCaract = sqrt(L**3 / (G * Mtotal))  ! Tiempo característico del sistema
    print *, "Tiempo característico:", tCaract/(3600.0d0*24.0d0), "días"

    MtCV = Mt / mTotal  ! Masa de la Tierra en términos de la masa total
    MlCV = Ml / mTotal  ! Masa de la Luna en términos de la masa
    MmCV = Mm / mTotal  ! Masa de Marte en términos de la masa total
    MsCV = Ms / mTotal  ! Masa del Sol en términos de la masa total

    ! Radios reescalados
    radioTierra = 6.371e6 / L ! Radio reescalado de la Tierra
    radioMarte = 3.3895e6 / L ! Radio reescalado de Marte
    radioLuna = 1.7374e6 / L  ! Radio reescalado de la Luna

    ! Velocidades angulares
    wt = tCaract*2.0 * pi / 31557600.0  ! Velocidad angular de la Tierra alrededor del Sol
    wm = tCaract*2.0 * pi / 59353000.0  ! Velocidad angular de Marte alrededor del Sol
    wl = tCaract*2.0 * pi / 2360592.0   ! Velocidad angular de la Luna alrededor de la Tierra

    ! Distancias orbitales reescaladas
    dOrbTierra = 1.0             ! Distancia orbital de la Tierra al Sol
    dOrbMarte = 2.279e11 / L     ! Distancia orbital de Marte al Sol
    dOrbLuna = 3.844e8 / L       ! Distancia orbital de la Luna a la Tierra

    ! Desfases iniciales
    desfIniTierra = 0.0          ! Desfase inicial de la Tierra
    desfIniMarte = pi/2.0d0           ! Desfase inicial de Marte
    desfIniLuna = 0.0            ! Desfase inicial de la Luna

    !----------------------------------------------------------------------------------
    ! VALORES INICIALES DEL COHETE
    altura = dOrbTierra*L + radioTierra*L + 250.0e3
    alturaTierra = altura - dOrbTierra*L
    r = altura/L
    !vOrbCohete = sqrt(G * Mt / altura) / L + wt * dOrbTierra * cos(wt * 0.0 + desfIniTierra)
    vOrbCohete = sqrt(G*Ms /altura)*tCaract/L + sqrt(G*Mt/alturaTierra)*tCaract/(L) 
    

    coordenadas(1) = r ! Posición inicial en el eje x
    coordenadas(2) = 0.0        ! Posición inicial en el eje y

    coordenadas(3) = 0.0! Velocidad inicial en el eje y
    coordenadas(4) = vOrbCohete ! Velocidad inicial en el eje x

    print *, sqrt(2*G*Mt/alturaTierra)*tCaract/L - sqrt(G*Mt/alturaTierra)*tCaract/L, "m/s"

    ! Cálculo de la velocidad necesaria para una transferencia de Hohmann hacia Marte

    ! Semi-eje mayor de la órbita de transferencia
    semiEjeMayor = (dOrbTierra + dOrbMarte) / 2.0

    ! Velocidad orbital de Marte en su órbita
    vOrbitaMarte = sqrt(G * Ms / (dOrbMarte * L)) * tCaract / L

    ! Energía orbital inicial en la órbita de la Tierra
    energiaOrbital = -G * Ms / (2.0 * dOrbTierra * L)

    ! Energía orbital en la órbita de transferencia
    energiaTransferencia = -G * Ms / (2.0 * semiEjeMayor * L)

    ! Velocidad necesaria para la transferencia de Hohmann
    vTransferencia = sqrt(2.0 * (energiaTransferencia + G * Ms / (dOrbTierra * L))) * tCaract / L &
                     - sqrt(G * Mt / (radioTierra * L)) * tCaract / L &
                     - wt * dOrbTierra

    print *, "Velocidad necesaria para la transferencia de Hohmann:", vTransferencia, "unidades normalizadas"
    print *, "Velocidad orbital de Marte:", vOrbitaMarte, "unidades normalizadas"

    ! Valores de impulso -----------------------------------------
    tiempoImpulsos = (/0.20828468252499205/)
    impulso = (/0.0/)
    !reccionImpulso = reshape((/cos(1.56347473), sin(1.56347473)/), (/1, 2/)) ! Dirección del impulso (x, y) basado en el ángulo
    do k = 1, size(tiempoImpulsos)
        boolImpulso(k) = 0 ! Inicializo el booleano de impulso
    end do

    print *, "Velocidad orbital del cohete:", impulso(1)


!-------------------------------------------
    ! Fichero de resultados
    unit_resultado = 10
    open(unit=unit_resultado, file='coordenadas.txt', status='replace')


    !  Defino el bucle principal

    tiempo = 0.0d0
    tiempoFinal = 15*24*3600.0d0/tCaract ! 12 meses en segundos
    step_size = 60.0d0/tCaract! El tiempo esta normalizado con la masa
    errorMaximo = step_size/100.0d0 ! Error máximo permitido
    j = 0
    !hMaxLimit = 3*60/m ! Límite máximo del paso de tiempo

    s(1) = 1.0d-8



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
            
            step_size = hMax

        else if (sMax < 2.0d0) then

            coordenadas = mediaCoordenada

            ! Escribo las coordenadas en un fichero
            coordTierra = calcularTierra(tiempo+step_size)
            coordMarte = calcularMarte(tiempo+step_size)
            coordLuna = calcularLuna(tiempo+step_size)
            coordLT = calcularLunaTierra(tiempo+step_size)

            if(mod(j, 1) == 0) then ! Cada 100 iteraciones
                write(unit_resultado, *) tiempo*tCaract/(3600.0d0), &
                    coordenadas(1), coordenadas(2), &
                    coordTierra(1), coordTierra(2), &
                    coordMarte(1), coordMarte(2), &
                    coordLuna(1), coordLuna(2), &
                    coordLT(1), coordLT(2)
            end if

            !do k = 1, size(tiempoImpulsos)
            !    if (tiempo >= tiempoImpulsos(k) .AND. boolImpulso(k) == 0) then
            !        ! Calcular la velocidad de la Tierra en el instante actual
            !        velocidadTierra(1) = -dOrbTierra * sin(wt * tiempo + desfIniTierra) * wt
            !        velocidadTierra(2) = dOrbTierra * cos(wt * tiempo + desfIniTierra) * wt
!
            !        ! Calcular la velocidad relativa del cohete respecto a la Tierra
            !        
            !        velocidadRelativa(1) = coordenadas(2) - velocidadTierra(1)
            !        velocidadRelativa(2) = coordenadas(3) - velocidadTierra(2)
!
            !        ! Normalizar la velocidad relativa para obtener la dirección del impulso
            !    
            !        magnitudVelocidadRelativa = sqrt(velocidadRelativa(1)**2 + velocidadRelativa(2)**2)
            !        direccionImpulso(1, 1) = velocidadRelativa(1) / magnitudVelocidadRelativa(1)
            !        direccionImpulso(1, 2) = velocidadRelativa(2) / magnitudVelocidadRelativa(2)
!
            !        ! Aplicar el impulso en la dirección calculada
            !        coordenadas(3) = coordenadas(3) + impulso(k) * direccionImpulso(1, 2)
            !        coordenadas(4) = coordenadas(4) + impulso(k) * direccionImpulso(1, 1)
            !        boolImpulso(k) = 1
            !        print *, "Impulso aplicado en el tiempo:", tiempo*tCaract/(3600.0d0), "s"
            !    end if
            !enddo

            !do k = 1, size(tiempoImpulsos)
            !    if (tiempo >= tiempoImpulsos(k) .AND. boolImpulso(k) == 0) then
            !        ! Calcular la velocidad de la Tierra en el instante actual
            !        velocidad = sqrt(coordenadas(3)**2 + coordenadas(4)**2)
            !        direccionImpulso(1, 1) = coordenadas(3) / velocidad
            !        direccionImpulso(1, 2) = coordenadas(4) / velocidad
!
            !        ! Aplicar el impulso en la dirección calculada
            !        coordenadas(3) = coordenadas(3) + impulso(k) * direccionImpulso(1, 1)
            !        coordenadas(4) = coordenadas(4) + impulso(k) * direccionImpulso(1, 2)
            !        boolImpulso(k) = 1
            !        print *, "Impulso aplicado en el tiempo:", tiempo*tCaract/(60.0d0), "min"
            !    end if
            !enddo

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
                    print *, "Impulso aplicado en el tiempo:", tiempo*tCaract/(60.0d0), "min"
                end if
            enddo


            


            tiempo = tiempo + step_size
            !if(mod(j, 100) == 0) then
            !    print*, tiempo*tCaract, " ", step_size
            !end if
            j = j + 1
            
            
            !if (mod(j, 100) == 0) then
            !    print *, "Tiempo:", tiempo * m, "s"
            !    print*, step_size * m, "s"
            !    print*, "------------------"
            !end if
            

            if (step_size < hMax) then
                step_size = step_size * 2.0d0
                !step_size = min(2 * step_size, hMaxLimit)
            end if
        
        end if
        

    end do

    print *, j


    close(unit_resultado)
    
    contains
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

end program marte