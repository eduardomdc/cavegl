#ifndef BOIDS_HPP
#define BOIDS_HPP

#include <vector>
#include "glm/ext/matrix_transform.hpp"
#include "glm/ext/vector_float3.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "cave.hpp"

class BoidSystem;

class Boid {
public:
    Boid();
    glm::vec3 pos;
    glm::vec3 vel;
    glm::vec3 up;
    float delta;
    BoidSystem* system;
    float phase; // wiggle phase
    void update(float delta, BoidSystem* bs);
};

class BoidSystem{
public:
    BoidSystem(int amount, Cave* cave);
    Cave* cave;
    std::vector<Boid*> boids;
    std::vector<float> verts;
    unsigned int VBO, VAO;
    void update(float delta);
    float neighboorRadius;
    float repellingForce;
    float cohesionForce;
    float alignmentForce;
    glm::mat4 scale;
};

#endif
