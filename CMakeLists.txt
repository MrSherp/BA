# add the executable
add_executable(main main.cpp)
target_link_libraries(main aol)
target_link_libraries (main ${Boost_LIBRARIES})
