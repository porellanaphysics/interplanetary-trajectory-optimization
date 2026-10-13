program optimizacion
    use omp_lib
    use simulacion2
    use randomnumber
    implicit none

    integer*8, parameter :: N = 2000, gen = 2000, numPadres = 3
    real*8, parameter :: probMutacion = 0.20
    real*8, allocatable :: poblacion(:,:), poblacionNueva(:,:), fitnessPoblacion(:)
    integer, allocatable :: indices(:)
    integer :: p, b
    real*8 :: random, resultado(5)
    integer :: c, num1, num2
    real*8, allocatable :: topFitness(:)
    integer*8 :: numEntero
    real*8 :: tempFitness
    integer :: tempIndex
    integer :: bestIndex
    real*8 :: bestFitness

    real*8 :: vMin, vMax, tMax

    call dran_ini(67657)

    ! Inicializamos el entorno para la paralelización
    !$OMP PARALLEL
    !$OMP SINGLE
    print *, "Número de hilos disponibles: ", omp_get_max_threads()
    !$OMP END SINGLE
    !$OMP END PARALLEL

    allocate(poblacion(N, N), poblacionNueva(N, N), fitnessPoblacion(N), indices(numPadres), topFitness(numPadres))

    ! Ponemos los límites de las variables para que escoja de manera aleatoria
    vMin = 4100*tCaract/L
    vMax = 5000*tCaract/L

    tMax = 1*30*24*3600/tCaract

    ! Ahora generamos la poblacion con valores aleatorios

    !! ESTO SE PARALELIZA
    !$OMP PARALLEL DO PRIVATE(random) SHARED(poblacion, vMin, vMax, tMax)
    do p = 1, N
        random = dran_u()
        !print*, random
        poblacion(p, 1) = vMin + (vMax - vMin)*random
        !print *, "Velocidad inicial: ", poblacion(p, 1)*L/tCaract
        random = dran_u()
        poblacion(p, 2) = tMax * dran_u()
    end do 
    !$OMP END PARALLEL DO

    
    ! Ahora iniciamos el blucle de generacion:

    do p = 1, gen
        !! ESTO SE PARALELIZA
        print*, "Calculando ", p
        !$OMP PARALLEL DO PRIVATE(resultado) SHARED(poblacion, fitnessPoblacion)
        do b = 1, N
            
            !print *, "C ", b, p, " | Velocidad: ", poblacion(b, 1), " Tiempo: ", poblacion(b, 2)
            resultado = calcularSimulacion(poblacion(b, 2), poblacion(b, 1))
            !print *, "F ", b, p, " | Velocidad: ", poblacion(b, 1), " Tiempo: ", poblacion(b, 2)
            
            if(resultado(1) <= dOrbLuna/10) then
                print *, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                print *, "¡Hemos llegado!"
                print *, "Distancia: ", resultado(1)
                print *, "Velocidad: ", poblacion(b, 1)
                print *, "Tiempo: ", poblacion(b, 2)
                print *, "V_Real: ", poblacion(b, 1)*L/tCaract
                print *, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
            end if

            ! PONGO LAS EXCEPCIONES
            if (resultado(5) < dOrbLuna*1.5) then ! esto es que no sale de la orbita terrestre
                fitnessPoblacion(b) = 1.0e-6
                !print *, "no ha salido de la tierra: ", b
            else if(resultado(4) < dOrbTierra*0.9) then ! limita que el cohete vaya al sistema solar interior
                fitnessPoblacion(b) = 1.0e-6
                !print *, "se ha ido al sol: ", b
            
            else if(resultado(4) > dOrbMarte*1.2) then ! limita que el cohete vaya al sistema solar interior
                fitnessPoblacion(b) = 1.0e-6

            else if(resultado(3) < 2*30*24*3600/tCaract) then !Limita que el cohete no lo alcance en menos de 3 meses
                fitnessPoblacion(b) = 1.0e-6
                !print *, "el tiempo es muy lento: ", b
            else 
                fitnessPoblacion(b) = fitness(resultado(1))
            end if

        end do
        !$OMP END PARALLEL DO

        print *, "Eligiendo ", p

        !RESULTADOS FINALES CUANDO SE LLEGA AL FINAL
        

        ! Seleccionamos los 5 valores más altos de fitnessPoblacion
        print *, "Seleccionando los mejores individuos de la generación ", p

        !topFitness = -1.0e10  ! Inicializamos con un valor muy bajo
        ! Inicializamos los valores de topFitness e indices
        topFitness = -1.0e10
        indices = -1

        ! Usamos un algoritmo simple para encontrar los mejores fitness
        do c = 1, numPadres
            

            bestFitness = -1.0e10
            bestIndex = -1

            do b = 1, N
            if (fitnessPoblacion(b) > bestFitness .and. all(indices /= b)) then
                bestFitness = fitnessPoblacion(b)
                bestIndex = b
            end if
            end do

            topFitness(c) = bestFitness
            indices(c) = bestIndex
        end do

        ! Ordenamos las listas de topFitness e indices en función de los valores de topFitness
        do c = 1, numPadres-1
            do b = c+1, numPadres
            if (topFitness(c) < topFitness(b)) then
                ! Intercambiamos los valores de topFitness
                tempFitness = topFitness(c)
                topFitness(c) = topFitness(b)
                topFitness(b) = tempFitness

                ! Intercambiamos los valores de indices
                
                tempIndex = indices(c)
                indices(c) = indices(b)
                indices(b) = tempIndex
            end if
            end do
        end do


        if(p == gen) then
            ! Encontramos el mejor individuo de la última generación
            topFitness(1) = -1.0e10
            indices(1) = -1

            do c = 1, N
                if (fitnessPoblacion(c) > topFitness(1)) then
                    topFitness(1) = fitnessPoblacion(c)
                    indices(1) = c
                end if
            end do

            ! Mostramos los resultados del mejor individuo
            print *, "Mejor individuo encontrado:"
            print *, "Fitness: ", topFitness(1)
            print *, "Velocidad: ", poblacion(indices(1), 1)
            print *, "Tiempo: ", poblacion(indices(1), 2)
            print*,  "V_Real: ",  poblacion(indices(1), 1)*L/tCaract

            ! Salimos del programa
            exit
        end if

        ! Imprimimos solo el mejor individuo seleccionado
        print *, "Mejor individuo seleccionado en la generación ", p, ":"
        print *, "  Fitness: ", topFitness(1)
        print *, "  Velocidad: ", poblacion(indices(1), 1)
        print *, "  Tiempo: ", poblacion(indices(1), 2)
        print *, "  V_Real: ", poblacion(indices(1), 1)*L/tCaract
        print *, "-------------------------------------------------------------"

        

        ! Ahora los índices contienen las posiciones de los numPadre valores más altos
        ! y top fitness sus valores

        print *, "Generando Individuos ", p

        !Ahora generamos una nueva poblacion a partir de estos
        !$OMP PARALLEL DO PRIVATE(num1, num2, random) SHARED(poblacionNueva, poblacion, indices, vMax, vMin)
        do b = 1, N-numPadres
            num1 = i_dran(numPadres)
            num2 = i_dran(numPadres)

            poblacionNueva(b, 1) = 0.5*(poblacion(indices(num1), 1) + poblacion(indices(num2), 1))
            poblacionNueva(b, 2) = 0.5*(poblacion(indices(num1), 2) + poblacion(indices(num2), 2))

            ! Aplicamos la mutacion
            random = dran_u()
            numEntero = 10
            if(random < probMutacion) then
            poblacionNueva(b,1) = poblacionNueva(b,1) + (vMax + vMin)*0.1*dran_u()*((-1)**i_dran(numEntero)) ! Puede variar un 10% de la velocidad
            if (poblacionNueva(b,1) > vMax) poblacionNueva(b,1) = vMax
            if (poblacionNueva(b,1) < vMin) poblacionNueva(b,1) = vMin
            endif

            random = dran_u()
            if(random < probMutacion) then 
            poblacionNueva(b,2) = poblacionNueva(b,2) + 24*3600*dran_u()/tCaract*((-1)**i_dran(numEntero)) ! Puede variar 1 dias
            if (poblacionNueva(b,2) > tMax) poblacionNueva(b,2) = tMax
            if (poblacionNueva(b,2) < 0.0) poblacionNueva(b,2) = 0.0
            endif

        end do
        !$OMP END PARALLEL DO

        ! Rellenamos la población nueva con los padres seleccionados
        do b = N-numPadres+1, N
            poblacionNueva(b, 1) = poblacion(indices(b-(N-numPadres)), 1)
            poblacionNueva(b, 2) = poblacion(indices(b-(N-numPadres)), 2)
        end do

        ! Actualizamos la población para la siguiente generación
        poblacion = poblacionNueva

        ! Imprimimos la generación actual
        !print *, "Generación ", p, ":"
        !do b = 1, N
        !    print *, "Individuo ", b, ": Velocidad = ", poblacion(b, 1), ", Tiempo = ", poblacion(b, 2)
        !    print *, "------"
        !end do
        !print *, "-------------------------------------------------------------"


    end do 

    ! Liberamos la memoria asignada
    deallocate(poblacion, poblacionNueva, fitnessPoblacion, indices, topFitness)

    ! Abrimos un archivo para escribir los resultados
    !open(unit=10, file="resultados.txt", status="unknown", action="write")
!
    !! Escribimos los resultados del mejor individuo
    !write(10, *) "Mejor individuo encontrado:"
    !write(10, *) "Fitness: ", topFitness(1)
    !write(10, *) "Velocidad: ", poblacion(indices(1), 1)
    !write(10, *) "Tiempo: ", poblacion(indices(1), 2)
    !write(10, *) "V_Real: ", poblacion(indices(1), 1)*L/tCaract
!
    !! Cerramos el archivo
    !close(10)


    contains
        function fitness(distancia) result(res)
            real*8 :: res, distancia

            res = 1/(1 + distancia)
        end function fitness
    
end program optimizacion