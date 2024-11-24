#include "app.hpp"
#include <GLFW/glfw3.h>
#include <cstdio>
#include <iostream>
#include "boids.hpp"
#include "glm/ext/matrix_float4x4.hpp"
#include "glm/ext/matrix_transform.hpp"
#include "external/stb_image.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/geometric.hpp"
#include "glm/trigonometric.hpp"
#include "cave.hpp"
#include <random>
#include <time.h>

App::App(int nthreads){
	camera = Camera();
    cave = new Cave(800, 50, 800, nthreads);
    boidsys = new BoidSystem(300, 1, cave);
    camera.pos = {75, 100, 75};
    camera.front = -glm::normalize(camera.pos-glm::vec3(cave->sizex/2.0f, 0.0f, cave->sizez/2.0f));
    texture1 = 0;
    fishtex = 0;
}

int load_texture(const char* file, unsigned int format, unsigned int index) {
    int width, height, nrChannels;
    unsigned int texture;
    
    unsigned char *data = stbi_load(file, &width, &height, &nrChannels, 0);
    if (!data) {
        std::cout<<"Failed to load texture data"<<std::endl;
        return -1;
    }
    glGenTextures(1, &texture); // texture array of size 1
    glActiveTexture(GL_TEXTURE0+index);
    glBindTexture(GL_TEXTURE_2D, texture);
    // image data to gl texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, format, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);
    // free image data
    stbi_image_free(data);
    return texture;
}

int App::init(){
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    // glfw window creation
    // --------------------
    window = glfwCreateWindow(800, 600, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    //enable Z buffering
    glEnable(GL_DEPTH_TEST);

	// texture wrapping mirrored
    // s,t,r == x,y,z
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    // texture filtering methods
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    // ----  CAVE GRAPHICS
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(
        GL_ARRAY_BUFFER,
        cave->data.size() * sizeof(float), 
        static_cast<void*>(cave->data.data()), 
        GL_STATIC_DRAW
    );

    /*
        vertex data attributes
        |   position    |   normal  |   texmap  |
        |   x, y, z     |   x, y, z |   x, y    |
    */
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture mapping attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6*sizeof(float)));
    glEnableVertexAttribArray(2);
    // safety unbinding
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    // LOADING TEXTURE ------------
    texture1 = load_texture("../textures/rock.jpg", GL_RGB, 0);
    fishtex = load_texture("../textures/clownfish.png", GL_RGBA, 1);
    // LOAD Caustics textures
    for (int i=1; i<=32; i++){
        char filename[100];
        sprintf(filename, "../textures/caustics/caustics_%03d.bmp", i);
        unsigned int name = load_texture(filename, GL_RGB, i+1);
        causticstex.push_back(name);
    }
    // set shader samplers for textures
    // CAVE SHADER- --------
    shader = new Shader("../src/shaders/v.vs", "../src/shaders/f.fs");
    shader->use();
    glUniform1i(glGetUniformLocation(shader->id, "texture1"), 0);
	glUniform3fv(glGetUniformLocation(shader->id, "objectColor"), 1, glm::value_ptr(cave->rockcolor));
    glUniform3fv(glGetUniformLocation(shader->id, "ambientLightColor"), 1, glm::value_ptr(cave->ambientLight->intensity*cave->ambientLight->color));
    glUniform3fv(glGetUniformLocation(shader->id, "lightPos"), 1, glm::value_ptr(cave->light->pos));
	glUniform3fv(glGetUniformLocation(shader->id, "lightColor"), 1, glm::value_ptr(cave->light->color));
    glUniform3fv(glGetUniformLocation(shader->id, "lightDir"), 1, glm::value_ptr(cave->light->dir));
    // model matrix
    glm::mat4 model = glm::mat4(1.0f);

	int modelLoc = glGetUniformLocation(shader->id, "model");
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model)); 

    //projection matrix
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.1f, 1000.0f);
	int projectionLoc = glGetUniformLocation(shader->id, "projection");
	glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
    // FISH SHADER --- --------
    fishader = new Shader("../src/shaders/fish.vs", "../src/shaders/fish.fs");
    fishader->use();
    glUniform1i(glGetUniformLocation(fishader->id, "texture1"), 1);
	glUniform3fv(glGetUniformLocation(fishader->id, "objectColor"), 1, glm::value_ptr(glm::vec3(1.0f)));
    glUniform3fv(glGetUniformLocation(fishader->id, "ambientLightColor"), 1, glm::value_ptr(cave->ambientLight->intensity*cave->ambientLight->color));
    glUniform3fv(glGetUniformLocation(fishader->id, "lightPos"), 1, glm::value_ptr(cave->light->pos));
	glUniform3fv(glGetUniformLocation(fishader->id, "lightColor"), 1, glm::value_ptr(cave->light->color));
    glUniform3fv(glGetUniformLocation(fishader->id, "lightDir"), 1, glm::value_ptr(cave->light->dir));
    // model matrix
    model = glm::mat4(1.0f);

	modelLoc = glGetUniformLocation(fishader->id, "model");
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model)); 

    //projection matrix
    projection = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.1f, 1000.0f);
	projectionLoc = glGetUniformLocation(fishader->id, "projection");
	glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));


    // ----- BOID GRAPHICS
    glGenVertexArrays(1, &boidsys->VAO);
    glGenBuffers(1, &boidsys->VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(boidsys->VAO);
    glBindBuffer(GL_ARRAY_BUFFER, boidsys->VBO);
    glBufferData(
        GL_ARRAY_BUFFER,
        boidsys->verts.size() * sizeof(float), 
        static_cast<void*>(boidsys->verts.data()), 
        GL_STATIC_DRAW
    );

    /*
        vertex data attributes
        |   position    |   normal  |   texmap  |
        |   x, y, z     |   x, y, z |   x, y    |
    */
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture mapping attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6*sizeof(float)));
    glEnableVertexAttribArray(2);
    // safety unbinding
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

	// capture mouse
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    //glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    //glEnable(GL_CULL_FACE);
    //glCullFace(GL_BACK);
    glEnable(GL_MULTISAMPLE);
    return 0;
}

int getCaustic(double time){
    return (int)(30*time)%32+2;
}

void App::render(){
	// rendering commands
    glm::vec3 ambient = cave->ambientLight->color*(cave->ambientLight->intensity*0.2f);
	glClearColor(ambient.x, ambient.y, ambient.z, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
	/* 
	trans = glm::mat4(1.0f);
	trans = glm::translate(trans, glm::vec3(0.5f, -0.5f, 0.0f));
	trans = glm::rotate(trans, (float)glfwGetTime(), glm::vec3(0.0f, 0.0f, 1.0f));*/
	//glUniformMatrix4fv(transformLoc, 1, GL_FALSE, glm::value_ptr(trans));
    
	shader->use();
    //model matrix
    glm::mat4 model = glm::mat4(1.0f);

	int modelLoc = glGetUniformLocation(shader->id, "model");
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    // view matrix
	int viewLoc = glGetUniformLocation(shader->id, "view");
	glm::mat4 view = camera.view();
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    fishader->use();
    glUniformMatrix4fv(glGetUniformLocation(fishader->id, "view"), 1, GL_FALSE, glm::value_ptr(view));
    unsigned int currentCaustic = getCaustic(glfwGetTime());
    shader->use();
    // draw cave
    glUniform1i(glGetUniformLocation(shader->id, "causticTex"), currentCaustic);
	glBindVertexArray(VAO);
	glDrawArrays(GL_TRIANGLES, 0, cave->vertices.size()/5);
    // draw fishies
    fishader->use();
    glBindVertexArray(boidsys->VAO);
	glUniform3fv(glGetUniformLocation(fishader->id, "objectColor"), 1, glm::value_ptr(glm::vec3(1.0f)));
    glUniform1i(glGetUniformLocation(fishader->id, "causticTex"), currentCaustic);
    for (int i=0; i< boidsys->boids.size(); i++){
        Boid* boid = boidsys->boids[i];
        glUniform1f(glGetUniformLocation(fishader->id, "phase"), glfwGetTime()+boid->phase);
        glm::mat4 fishmodel = glm::mat4(1.0f);
        fishmodel = glm::translate(fishmodel, glm::vec3(-1.5f, -1.5f, -2.5f));
        glm::mat4 trans = glm::translate(glm::mat4(1.0f), boid->pos);
        glm::mat4 look = glm::inverse(glm::lookAt(glm::vec3(0.0f), glm::normalize(-boid->vel), glm::vec3(0.0f, 1.0f, 0.0f)));
        fishmodel = trans*look*boid->scale*fishmodel;
        int modelLoc = glGetUniformLocation(fishader->id, "model");
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(fishmodel));
        glDrawArrays(GL_TRIANGLES, 0, boidsys->verts.size()/8);
    }
	
	glUniform3fv(glGetUniformLocation(fishader->id, "objectColor"), 1, glm::value_ptr(glm::vec3(0.8f, 0.5f, 0.5f)));
    for (int i=0; i< boidsys->killers.size(); i++){
        Killer* boid = boidsys->killers[i];
        glUniform1f(glGetUniformLocation(fishader->id, "phase"), glfwGetTime()+boid->phase);
        glm::mat4 fishmodel = glm::mat4(1.0f);
        fishmodel = glm::translate(fishmodel, glm::vec3(-1.5f, -1.5f, -2.5f));
        glm::mat4 trans = glm::translate(glm::mat4(1.0f), boid->pos);
        glm::mat4 look = glm::inverse(glm::lookAt(glm::vec3(0.0f), glm::normalize(-boid->vel), glm::vec3(0.0f, 1.0f, 0.0f)));
        fishmodel = trans*look*boid->scale*fishmodel;
        int modelLoc = glGetUniformLocation(fishader->id, "model");
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(fishmodel));
        glDrawArrays(GL_TRIANGLES, 0, boidsys->verts.size()/8);
    }
    //cave->light->pos = boidsys->killers[0]->pos;
	// check and call events and swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;
}

void App::input(){
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
	
    const float cameraSpeed = 50.0f*deltaTime;
    
    glm::vec3 right = glm::normalize(glm::cross(camera.front, camera.up));
    
    if (camera.fps){
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            camera.pos += cameraSpeed * camera.front;
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            camera.pos -= cameraSpeed * camera.front;
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            camera.pos -= right * cameraSpeed;
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            camera.pos += right * cameraSpeed;
    } else {
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            camera.pos += cameraSpeed * glm::normalize(glm::cross(camera.front, -right));
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            camera.pos -= cameraSpeed * glm::normalize(glm::cross(camera.front, -right));
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            camera.pos -= right * cameraSpeed;
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            camera.pos += right * cameraSpeed;
        glm::vec3 direction;
        direction = -camera.pos+glm::vec3(cave->sizex/2, cave->sizey/3.0, cave->sizez/2);
        camera.front = glm::normalize(direction);
    }
}

void App::update(){
    //camera.update(deltaTime, cave);
    boidsys->update(deltaTime);
}
