# header files
set( WRITERS_HEADER_BASENAMES
    KSWritersMessage.h
    KSWriteASCII.h
)
if( Kassiopeia_USE_VTK )
    set( WRITERS_HEADER_BASENAMES
        ${WRITER_HEADER_BASENAMES}
        KSWriteVTK.h
    )
endif( Kassiopeia_USE_VTK )
if( Kassiopeia_USE_ROOT )
    set( WRITERS_HEADER_BASENAMES
        ${WRITER_HEADER_BASENAMES}
        KSWriteROOT.h
        KSWriteROOTCondition.h
        KSWriteROOTConditionOutput.h
        KSWriteROOTConditionPeriodic.h
        KSWriteROOTConditionTerminator.h
        KSWriteROOTConditionStep.h
    )
endif( Kassiopeia_USE_ROOT )
set( WRITERS_HEADER_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/Include
)
foreach( BASENAME ${WRITERS_HEADER_BASENAMES} )
    list( APPEND WRITERS_HEADER_FILENAMES ${WRITERS_HEADER_PATH}/${BASENAME} )
endforeach( BASENAME )

# source files
set( WRITERS_SOURCE_BASENAMES
    KSWritersMessage.cxx
    KSWriteASCII.cxx
)
if( Kassiopeia_USE_VTK )
    set( WRITERS_SOURCE_BASENAMES
        ${WRITERS_SOURCE_BASENAMES}
        KSWriteVTK.cxx
    )
endif( Kassiopeia_USE_VTK )
if( Kassiopeia_USE_ROOT )
    set( WRITERS_SOURCE_BASENAMES
        ${WRITERS_SOURCE_BASENAMES}
        KSWriteROOT.cxx
        KSWriteROOTCondition.cxx
    )
endif( Kassiopeia_USE_ROOT )
set( WRITERS_SOURCE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/Source
)
foreach( BASENAME ${WRITERS_SOURCE_BASENAMES} )
    list( APPEND WRITERS_SOURCE_FILENAMES ${WRITERS_SOURCE_PATH}/${BASENAME} )
endforeach( BASENAME )

# data files
set( WRITERS_DATA_BASENAMES
)
set( WRITERS_DATA_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/Data
)
foreach( BASENAME ${WRITERS_DATA_BASENAMES} )
    list( APPEND WRITERS_DATA_FILENAMES ${WRITERS_DATA_PATH}/${BASENAME} )
endforeach( BASENAME )

# library
add_library( KassiopeiaWriters SHARED
    ${WRITERS_SOURCE_FILENAMES} ${WRITERS_HEADER_FILENAMES} )
target_include_directories(KassiopeiaWriters
    PUBLIC $<BUILD_INTERFACE:${WRITERS_HEADER_PATH}> $<INSTALL_INTERFACE:include> )
target_link_libraries( KassiopeiaWriters
    PUBLIC
        KassiopeiaOperators
        ROOT::Tree
        ROOT::Physics
)

if( Kassiopeia_USE_VTK )
    target_link_libraries( KassiopeiaWriters PUBLIC KommonVtk )
    if( VTK_MAJOR_VERSION GREATER 5 )
        target_compile_definitions( KassiopeiaWriters PRIVATE VTK6 )
    endif()
endif( Kassiopeia_USE_VTK )

# install
kasper_install_headers( ${WRITERS_HEADER_FILENAMES} )
kasper_install_data( ${WRITERS_DATA_FILENAMES} )
kasper_install_libraries( KassiopeiaWriters )
