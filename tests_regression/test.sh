#!/bin/bash

spack env activate wass

startdir=$(pwd)
cooldown=60

cd ..
mkdir -p build
cd build

cmake ..
cmake --build . -j4

#experiments (to be done both in direct solve and fma; varying mpi processes)
#-- plain scalar phi, v1
#-- plain scalar phi, v2; different function
#-- basic 2-component phi, (v1, v2)

nthreads=1
ncores_total=4

#modes="direct fma"
modes="direct fma"

echo "reference problems (np=1 fixed)"
#generate the baseline - this should be skipped if ref* files are already present
#logs are saved as ref_* while all other isntances are reg_*
for mode in ${modes}
do
    echo "solver mode ${mode}"
    for func in 1 2
    do
        cp ${startdir}/ref_parameters_bem_3_f${func}_${mode}.prm ./parameters_bem_3.prm

        #references should be for np=1; np>1 are regressions
        nprocs=1
        nthreads=1
        
        if [ ! -f "${startdir}/ref_f${func}_${mode}.log" ]
        then
            echo "reference problem ${func} - nprocs ${nprocs} - nthreads ${nthreads} - solver ${mode}"
            export OMP_NUM_THREADS=${nthreads}
            perf stat --detailed mpirun --np ${nprocs} --map-by node:PE=1 --bind-to none bem_fma_3d ${nthreads} > log 2>&1
            
            result=${?}

            if (( result == 0 ))
            then
                for fn in result_scalar_results result_vector_results #scalar_error vector_error
                do
                    cp ./${fn}.vtu ${startdir}/ref_${fn}_f${func}_${mode}.vtu
                    
                    echo "reg_np${nprocs}_${fn}_f${func}_${mode}.vtu"
                    meshio-info ${startdir}/ref_${fn}_f${func}_${mode}.vtu
                done
                
                #the first is for simple reading, the one appending is for python mega pandas dataframing
                cp log ${startdir}/ref_f${func}_${mode}.log
                cat log >> ${startdir}/ref_f${func}_${mode}_log.log
            else
                echo "error during reference problem ${func} - nprocs ${nprocs} - nthreads ${nthreads} - solver ${mode}"
                
                #the first is for simple reading, the one appending is for python mega pandas dataframing
                cp log ${startdir}/ref_f${func}_${mode}.err
                cat log >> ${startdir}/ref_f${func}_${mode}_log.err
            fi
            
            echo
            sleep ${cooldown}
        fi
    done
done

#exit 0

nthreads=1

#solve the simple problem instances with multiple mpi procs
#logs are saved as reg_*
for mode in ${modes}
do
    echo "regression problems of mode ${mode}"
    for func in 1 2
    do
        cp ${startdir}/ref_parameters_bem_3_f${func}_${mode}.prm ./parameters_bem_3.prm

        #references should be for np=1; np>1 are regressions
        for nprocs in 2 4
        do 
            echo "regression problem ${func} - nprocs ${nprocs} - nthreads ${nthreads} - solver ${mode}"
            export OMP_NUM_THREADS=${nthreads}
            if ((nprocs > 1))
            then
                perf stat --detailed mpirun --np ${nprocs} --map-by node:PE=${nthreads} bem_fma_3d ${nthreads} > log 2>&1
            else
                perf stat --detailed mpirun --np ${nprocs} --map-by node:PE=${nthreads} --bind-to none bem_fma_3d ${nthreads} > log 2>&1
            fi

            result=${?}

            if (( result == 0 ))
            then
                for fn in result_scalar_results result_vector_results #scalar_error vector_error
                do
                    cp ./${fn}.vtu ${startdir}/reg_np${nprocs}_nt${nthreads}_${fn}_f${func}_${mode}.vtu
                    
                    echo "reg_np${nprocs}_nt${nthreads}_${fn}_f${func}_${mode}.vtu"
                    meshio-info ${startdir}/reg_np${nprocs}_nt${nthreads}_${fn}_f${func}_${mode}.vtu
                done
                
                #the first is for simple reading, the one appending is for python mega pandas dataframing
                cp log ${startdir}/reg_np${nprocs}_nt${nthreads}_f${func}_${mode}.log
                cat log >> ${startdir}/reg_np${nprocs}_nt${nthreads}_f${func}_${mode}_log.log
            else
                echo "error during regression problem ${func} - nprocs ${nprocs} - nthreads ${nthreads} - solver ${mode}"
                
                #the first is for simple reading, the one appending is for python mega pandas dataframing
                cp log ${startdir}/reg_np${nprocs}_nt${nthreads}_f${func}_${mode}.err
                cat log >> ${startdir}/reg_np${nprocs}_nt${nthreads}_f${func}_${mode}_log.err
            fi
            echo
            sleep ${cooldown}
        done
    done
done

#exit 0

#solve complex problems - each instance shall generate two sets of result files, with f1 and f2 suffix to be compared with references
for mode in ${modes}
do
    cp ${startdir}/ref_parameters_bem_3_complex_${mode}.prm ./parameters_bem_3.prm

    for nprocs in 1 2 4
    do 
        echo "complex problem - nprocs ${nprocs} - nthreads ${nthreads} - solver ${mode}"
        export OMP_NUM_THREADS=${nthreads}
        if ((nprocs > 1))
        then
            perf stat --detailed mpirun --np ${nprocs} --map-by node:PE=${nthreads} bem_fma_3d ${nthreads} > log 2>&1
        else
            perf stat --detailed mpirun --np ${nprocs} --map-by node:PE=${nthreads} --bind-to none bem_fma_3d ${nthreads} > log 2>&1
        fi

        result=${?}

        if (( result == 0 ))
        then
            #first component has no suffix on filename
            for fn in result_scalar_results result_vector_results #scalar_error vector_error
            do
                cp ./${fn}.vtu ${startdir}/reg_np${nprocs}_nt${nthreads}_${fn}_complex_f1_${mode}.vtu
                
                echo "reg_np${nprocs}_${fn}_complex_f1.vtu"
                meshio-info ${startdir}/reg_np${nprocs}_nt${nthreads}_${fn}_complex_f1_${mode}.vtu
            done
            #do the remaining comps
            for func in 2
            do
                for fn in result_${func}_scalar_results result_${func}_vector_results #scalar_${func}_error vector_${func}_error
                do
                    fn_redux="${fn/_${func}/}"
                    cp ./${fn}.vtu ${startdir}/reg_np${nprocs}_nt${nthreads}_${fn_redux}_complex_f${func}_${mode}.vtu
                    
                    echo "reg_np${nprocs}_nt${nthreads}_${fn_redux}_complex_f${func}.vtu"
                    meshio-info ${startdir}/reg_np${nprocs}_nt${nthreads}_${fn_redux}_complex_f${func}_${mode}.vtu
                done
            done
            
            #the first is for simple reading, the one appending is for python mega pandas dataframing
            cp log ${startdir}/reg_np${nprocs}_nt${nthreads}_complex_${mode}.log
            cat log >> ${startdir}/reg_np${nprocs}_nt${nthreads}_complex_${mode}_log.log
        else
            echo "error during complex problem - nprocs ${nprocs} - nthreads ${nthreads} - solver ${mode}"
            
            #the first is for simple reading, the one appending is for python mega pandas dataframing
            cp log ${startdir}/reg_np${nprocs}_nt${nthreads}_complex_${mode}.err
            cat log >> ${startdir}/reg_np${nprocs}_nt${nthreads}_complex_${mode}_log.err
        fi
        echo
        sleep ${cooldown}
    done
done

#exit 0

cd ${startdir}
./test_compares.sh
