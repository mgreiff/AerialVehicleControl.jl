cmake_minimum_required (VERSION 3.7)

project(control)

include_directories(include)

file(GLOB SOURCES "src/*.c")

file(GLOB HEADERS "include/*.h")

message("source files: " ${SOURCES})

add_library(${PROJECT_NAME}_lib STATIC ${SOURCES} ${HEADERS})
target_link_libraries (${PROJECT_NAME}_lib -lm)

install (TARGETS ${PROJECT_NAME}_lib DESTINATION lib)
install (FILES ${HEADERS} ${CL} DESTINATION include/${PROJECT_NAME})