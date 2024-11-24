#include <cstddef>
#include <ctime>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "glm/ext/matrix_transform.hpp"
#include "glm/ext/vector_float3.hpp"
#include "glm/geometric.hpp"
#include "glm/trigonometric.hpp"
#include "external/stb_image.h"
#include "shader.hpp"
#include "math.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "app.hpp"
#include <random>
#include <time.h>

App* app;

void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

int main(int argc, char** argv){
    int nthreads = atoi(argv[1]);
    unsigned long int seed = time(nullptr);
    srand(seed);
    app = new App(nthreads);
    app->init();
    glfwTerminate();
    return 0;
    glfwSetCursorPosCallback(app->window, mouse_callback);
    glfwSetFramebufferSizeCallback(app->window, framebuffer_size_callback);
    glfwSetKeyCallback(app->window, key_callback);

    //render loop
    while(!glfwWindowShouldClose(app->window)){
        // input
        app->input();
        // rendering commands
        app->render();
        // update physics
        app->update();
    }
    // deallocate resources
    glDeleteVertexArrays(1, &app->VAO);
    glDeleteBuffers(1, &app->VBO);
    glDeleteBuffers(1, &app->EBO);
    glDeleteProgram(app->shader->id);
    delete app->cave;
    delete app->boidsys;
    //close glfw
    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height){
    glViewport(0,0,width,height);
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), float(width) / float(height), 0.1f, 1000.0f);
    app->shader->use();
	int projectionLoc = glGetUniformLocation(app->shader->id, "projection");
	glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
    app->fishader->use();
    projectionLoc = glGetUniformLocation(app->fishader->id, "projection");
	glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
    
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos){
    float xoffset = xpos - app->lastX;
    float yoffset = app->lastY - ypos;
	
    app->lastX = xpos;
    app->lastY = ypos;
    const float sensitivity = 0.05f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;
    if (app->camera.fps){
        app->camera.yaw += xoffset;
        app->camera.pitch += yoffset;
    }

    if (app->camera.pitch > 89) app->camera.pitch = 89;
    if (app->camera.pitch < -89) app->camera.pitch = -89;

    glm::vec3 direction;
    
    direction.x = cos(glm::radians(app->camera.yaw )) * cos(glm::radians(app->camera.pitch));
    direction.y = sin(glm::radians(app->camera.pitch));
    direction.z = sin(glm::radians(app->camera.yaw)) * cos(glm::radians(app->camera.pitch));
    
    //direction = -app->camera.pos;
    if (app->camera.fps){
        app->camera.front = glm::normalize(direction);
    } else {
        glm::vec3 right = glm::normalize(glm::cross(app->camera.front, app->camera.up));
        app->camera.pos += right*xoffset;
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS){
        app->camera.fps = !app->camera.fps;
    }
}