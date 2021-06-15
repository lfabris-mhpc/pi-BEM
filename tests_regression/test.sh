#!/bin/bash

spack env activate wass

startdir=$(pwd)

cd ..
mkdir -p build
cd build

cmake ..
cmake --build . -j4

#experiments (to be done both in direct solve and fma; varying mpi processes)
#-- plain scalar phi, v1
#-- plain scalar phi, v2; different function
#-- basic 2-component phi, (v1, v2)

tbb_nthreads=1

echo "reference problems (np=1 fixed)"
#generate the baseline - this should be skipped if ref* files are already present
#logs are saved as ref_* while all other isntances are reg_*
for mode in direct fma
do
    echo "solver mode ${mode}"
    for func in 1 2
    do
        cp ${startdir}/ref_parameters_bem_3_f${func}_${mode}.prm ./parameters_bem_3.prm

        #references should be for np=1; np>1 are regressions
        for np in 1
        do
            if [ ! -f "${startdir}/ref_f${func}_${mode}_log.log" ]
            then
                echo "reference problem ${func} - np ${np} - ${mode}"
                export OMP_NUM_THREADS=${tbb_nthreads}
                perf stat --detailed mpirun --np ${np} --map-by node:PE=1 --bind-to none bem_fma_3d ${tbb_nthreads} > log 2>&1
                
                result=${?}
    
                if (( result == 0 ))
                then
                    for fn in result_scalar_results result_vector_results #scalar_error vector_error
                    do
                        cp ./${fn}.vtu ${startdir}/ref_${fn}_f${func}_${mode}.vtu
                        
                        echo "reg_np${np}_${fn}_f${func}_${mode}.vtu"
                        meshio-info ${startdir}/ref_${fn}_f${func}_${mode}.vtu
                    done
                    
                    #awk, grep or cp from log to the startdir for the info we're interested in
                    cp log ${startdir}/ref_f${func}_${mode}_log.log
                else
                    echo "error during reference problem ${func} - np ${np} - ${mode}"
                fi
                echo
            fi
        done
    done
done

#exit 0

#solve the simple problem instances with multiple mpi procs
#logs are saved as reg_*
for mode in direct fma
do
    echo "regression problems of mode ${mode}"
    for func in 1 2
    do
        cp ${startdir}/ref_parameters_bem_3_f${func}_${mode}.prm ./parameters_bem_3.prm

        #references should be for np=1; np>1 are regressions
        for np in 2 4
        do 
            echo "regression problem ${func} - np ${np} - ${mode}"
            export OMP_NUM_THREADS=${tbb_nthreads}
            perf stat --detailed mpirun --np ${np} --map-by node:PE=1 --bind-to none bem_fma_3d ${tbb_nthreads} > log 2>&1

            result=${?}

            if (( result == 0 ))
            then
                for fn in result_scalar_results result_vector_results #scalar_error vector_error
                do
                    cp ./${fn}.vtu ${startdir}/reg_np${np}_${fn}_f${func}_${mode}.vtu
                    
                    echo "reg_np${np}_${fn}_f${func}_${mode}.vtu"
                    meshio-info ${startdir}/reg_np${np}_${fn}_f${func}_${mode}.vtu
                done
                
                #awk, grep or cp from log to the startdir for the info we're interested in
                cp log ${startdir}/reg_np${np}_f${func}_${mode}_log.log
            else
                echo "error during regression problem ${func} - np ${np} - ${mode}"
            fi
            echo
        done
    done
done

#exit 0

#solve complex problems - each instance shall generate two sets of result files, with f1 and f2 suffix to be compared with references
for mode in direct fma
do
    cp ${startdir}/ref_parameters_bem_3_complex_${mode}.prm ./parameters_bem_3.prm

    for np in 1 2 4
    do 
        echo "complex problem - np ${np} - ${mode}"
        export OMP_NUM_THREADS=${tbb_nthreads}
	    perf stat --detailed mpirun --np ${np} --map-by node:PE=1 --bind-to none bem_fma_3d ${tbb_nthreads} > log 2>&1

        result=${?}

        if (( result == 0 ))
        then
            #first component has no suffix on filename
            for fn in result_scalar_results result_vector_results #scalar_error vector_error
            do
                cp ./${fn}.vtu ${startdir}/reg_np${np}_${fn}_complex_f1_${mode}.vtu
                
                echo "reg_np${np}_${fn}_complex_f1.vtu"
                meshio-info ${startdir}/reg_np${np}_${fn}_complex_f1_${mode}.vtu
            done
            #do the remaining comps
            for func in 2
            do
                for fn in result_${func}_scalar_results result_${func}_vector_results #scalar_${func}_error vector_${func}_error
                do
                    fn_redux="${fn/_${func}/}"
                    cp ./${fn}.vtu ${startdir}/reg_np${np}_${fn_redux}_complex_f${func}_${mode}.vtu
                    
                    echo "reg_np${np}_${fn_redux}_complex_f${func}.vtu"
                    meshio-info ${startdir}/reg_np${np}_${fn_redux}_complex_f${func}_${mode}.vtu
                done
            done
            
            #awk, grep or cp from log to the startdir for the info we're interested in
            cp log ${startdir}/reg_np${np}_complex_${mode}_log.log
        else
            echo "error during complex problem - np ${np} - ${mode}"
        fi
        echo
    done
done

#exit 0

cd ${startdir}
#do the comparison - loop on np
#python mesh_compare.py ref_result_scalar_results_f1_direct.vtu ref_result_scalar_results_f1_fma.vtu 0.005
#python mesh_compare.py ref_result_scalar_results_f1_direct.vtu reg_np1_result_scalar_results_complex_f1_direct.vtu 0.005
#python mesh_compare.py ref_result_scalar_results_f2_direct.vtu reg_np1_result_scalar_results_complex_f2_direct.vtu 0.005
#python mesh_compare.py ref_result_scalar_results_f1_direct.vtu reg_np1_result_scalar_results_complex_f1_fma.vtu 0.005
#python mesh_compare.py ref_result_scalar_results_f2_direct.vtu reg_np1_result_scalar_results_complex_f2_fma.vtu 0.005
#python mesh_compare.py ref_result_scalar_results_f1_fma.vtu reg_np1_result_scalar_results_complex_f1_direct.vtu 0.005
#python mesh_compare.py ref_result_scalar_results_f2_fma.vtu reg_np1_result_scalar_results_complex_f2_direct.vtu 0.005

atol=0.005
echo "regressions instances"
for mode in direct fma
do
    for func in 1 2
    do
        for np in 2 4
        do
            for type in result_scalar_results result_vector_results #scalar_error vector_error
            do
                python mesh_compare.py ref_${type}_f${func}_${mode}.vtu reg_np${np}_${type}_f${func}_${mode}.vtu ${atol}
            done
        done
    done
done

echo "complex instances"
for mode in direct fma
do
    for func in 1 2
    do
        for np in 1 2 4
        do
            for type in result_scalar_results result_vector_results #scalar_error vector_error
            do
                python mesh_compare.py ref_${type}_f${func}_${mode}.vtu reg_np${np}_${type}_complex_f${func}_${mode}.vtu ${atol}
            done
        done
    done
done