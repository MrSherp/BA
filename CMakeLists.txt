add_executable(main main.cpp)
target_link_libraries(main LINK_PUBLIC ${Boost_LIBRARIES} LINK_PUBLIC generallib)

IF ( USE_OPENMP )
target_compile_definitions(main PUBLIC USE_OPENMP)
set_target_properties(main PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
ENDIF ( )
