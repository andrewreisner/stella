add_subdirectory(single-block)

add_executable(driver-sb ${CMAKE_CURRENT_SOURCE_DIR}/single-block/driver.c)
target_link_libraries(driver-sb single-block stella)

add_executable(debug-sb ${CMAKE_CURRENT_SOURCE_DIR}/single-block/debug.c)
target_link_libraries(debug-sb single-block stella)
