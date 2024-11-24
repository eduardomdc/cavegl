cmake_minimum_required(VERSION 3.9)
project (
    learn-opengl
    LANGUAGES CXX
)

set(OpenGL_GL_PREFERENCE LEGACY)
set(GLAD_INCLUDE_DIRS /usr/include/glad)
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(OpenMP)

add_subdirectory(glad)

include_directories(
    ${OPENGL_INCLUDE_DIRS} 
    ${GLFW_INCLUDE_DIRS}
    #${GLAD_INCLUDE_DIRS}
)

set(MY_SOURCE_FILES
    ../src/main.cpp
    ../src/external/image_loader.cpp
    ../src/camera.cpp
    ../src/app.cpp
    ../src/cave.cpp
    ../src/marchingcubes.cpp
    ../src/perlin.cpp
    ../src/external/simplex.cpp
    ../src/light.cpp
    ../src/boids.cpp
)

link_directories(${GLAD_INCLUDE_DIRS})

add_executable(main ${MY_SOURCE_FILES})

target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX ${OPENGL_LIBRARIES} glfw glad)
