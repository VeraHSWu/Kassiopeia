option (KGeoBag_ENABLE_APP "Build KGeoBag applications" ON)
if (KGeoBag_ENABLE_APP)

    find_package(OpenMP)  # optional

    # names
    set( APPLICATION_NAMES
        KGeoBag
        RandomExtensionTester
        MetricsExtensionTester
        EulerAngles
        GeometryDump
        MergeTriangles
    )

    # libraries
    set( APPLICATION_LIBRARIES
        KGeoBagBindings
        KGeoBagCore
        KGeoBagRandom
    )

    if( KGeoBag_USE_VTK )
        list( APPEND APPLICATION_NAMES
            GeometryViewer
            MeshViewer
            AxialMeshViewer
            PointTester
            NormalTester
            DistanceTester
            OutsideTester
            RandomPointTester
            MeshIntersectionTester
        )
        list( APPEND APPLICATION_LIBRARIES
            KGeoBagVtkVisualization
        )
    endif( KGeoBag_USE_VTK )

    if( KGeoBag_USE_ROOT )
        list( APPEND APPLICATION_NAMES
            GeometryViewerROOT
        )
        list( APPEND APPLICATION_LIBRARIES
            KGeoBagROOTVisualization
        )
    endif( KGeoBag_USE_ROOT )

    # target
    foreach( NAME ${APPLICATION_NAMES} )
        add_executable( ${NAME} Source/${NAME}.cc )
    endforeach()

    # install
    foreach( NAME ${APPLICATION_NAMES} )
        target_link_libraries( ${NAME} ${APPLICATION_LIBRARIES} )
        if(OpenMP_FOUND)
            target_compile_definitions(${NAME} PRIVATE USE_OPENMP)
            target_link_libraries(${NAME} OpenMP::OpenMP_CXX)
        endif(OpenMP_FOUND)
        kasper_install_executables( ${NAME} )
    endforeach()

endif(KGeoBag_ENABLE_APP)
