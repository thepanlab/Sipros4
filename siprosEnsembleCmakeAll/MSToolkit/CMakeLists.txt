aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/expat-2.2.9 expatSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/MSToolkit MSToolkitSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/mzIMLTools mzIMLToolsSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/mzParser mzParserSRC)
## add sqlite
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/sqlite-3.32.1 sqliteSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/zlib-1.2.11 zlibSRC)

## add sqlite
# add_library(MStoolkit ${expatSRC} ${MSToolkitSRC} ${mzIMLToolsSRC} ${mzParserSRC} ${sqliteSRC} ${zlibSRC})
add_library(MStoolkit ${expatSRC} ${MSToolkitSRC} ${mzIMLToolsSRC} ${mzParserSRC} ${zlibSRC})
## for .so lib
# set_target_properties(MStoolkit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
## for .a lib
# set_target_properties(MStoolkit PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_compile_definitions(MStoolkit PUBLIC GCC HAVE_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)
## add sqlite
# target_compile_definitions(MStoolkit PUBLIC GCC HAVE_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW)

