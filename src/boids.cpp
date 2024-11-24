#include "boids.hpp"
#include "glm/ext/matrix_transform.hpp"
#include "glm/ext/quaternion_geometric.hpp"
#include "glm/geometric.hpp"
#include "glm/gtc/random.hpp"
#include <iostream>
#include <vector>

Boid::Boid(){
    up = glm::vec3(0.0f, 1.0f, 0.0f);
    phase = 1.5*randomf();
}

Killer::Killer(){
    up = glm::vec3(0.0f, 1.0f, 0.0f);
    phase = 1.5*randomf();
}

bool inVisionCone(Boid* boid, Boid* neigh){
    // returns true if neigh is in vision cone of boid
    // we can check the dot product between normalized boid->neigh and boid velocity
    // the closer to -1 the more neigh is directly behind boid.
    glm::vec3 dir = glm::normalize(neigh->pos-boid->pos);
    if (glm::dot(dir, glm::normalize(boid->vel)) < -0.2f) return false;
    else return true;
}

glm::vec3 collisionResponse(glm::vec3 pos, BoidSystem* bs, bool* hit){
    *hit = true;
    if (pos.x < 0.0f){
        return glm::vec3(1.0f, 0.0f, 0.0f);
    }
    if (pos.x > bs->cave->sizex -0.0f){
        return glm::vec3(-1.0f, 0.0f, 0.0f);
    }
    if (pos.y < 0.0f){
        return glm::vec3(0.0f, 1.0f, 0.0f);
    }
    if (pos.y > bs->cave->sizey-0.0f){
        return glm::vec3(0.0f, -1.0f, 0.0f);
    }
    if (pos.z < 0.0f){
        return glm::vec3(0.0f, 0.0f, 1.0f);
    }
    if (pos.z > bs->cave->sizez-0.0f){
        return glm::vec3(0.0f, 0.0f, -1.0f);
    }
    // terrain avoidance
    if (bs->cave->field(pos.x, pos.y, pos.z) < bs->cave->level){
        return glm::normalize(bs->cave->fieldGrad(pos));
    }
    *hit = false;
    return glm::vec3(0.0f, 0.0f, 0.0f);
}

void Boid::update(float delta){
    
    unsigned int count = 0; // number of neighbors
    glm::vec3 repelling = glm::vec3(0.0f, 0.0f, 0.0f); // repelling component of velocity
    glm::vec3 cohesion = glm::vec3(0.0f, 0.0f, 0.0f); // cohesion component of velocity
    glm::vec3 alignment = glm::vec3(0.0f, 0.0f, 0.0f); // alignment component
    glm::vec3 averagePos = glm::vec3(0.0f, 0.0f, 0.0f); // average position of neighbors
    glm::vec3 averageVel = glm::vec3(0.0f, 0.0f, 0.0f); // average velocity of neighbors
    bool predator = false;
    
    for (int i=0; i<bs->boids.size(); i++){
        if (bs->boids[i] == this) continue; // no self interaction!
        float distance = glm::distance(bs->boids[i]->pos, pos);

        // find close boids
        if (distance < bs->neighboorRadius && inVisionCone(this, bs->boids[i])){
            count++;
            //repelling component
            glm::vec3 toBoid = bs->boids[i]->pos - pos;
            repelling += (1/distance)*bs->repellingForce*glm::normalize(toBoid);
            // averaging position
            averagePos += bs->boids[i]->pos;
            // averaging velocity
            averageVel += bs->boids[i]->vel;
        } else continue;
    }

    // predator avoidance
    glm::vec3 runaway = glm::vec3(0.0f);
    for (int i=0; i<bs->killers.size(); i++){
        float distance = glm::distance(bs->killers[i]->pos, pos);
        // find close predators 
        if (distance < 2*bs->neighboorRadius){
            //repelling component
            predator = true;
            glm::vec3 toKiller = bs->killers[i]->pos - pos;
            runaway += glm::normalize(-toKiller);
        } else continue;
    }

    // collision avoidance
    // boundaries
    glm::vec3 avoidance = glm::vec3(0.0f, 0.0f , 0.0f); // avoidance component
    int maxdistance = 5;
    glm::vec3 dir = glm::normalize(vel);
    glm::vec3 raycast = pos; // raycast starts at boid pos
    int terrain = 0; // is there a predicted terrain collision
    bool hit = false;
    
    //raycast loop increment distance each iteration
    for (int distance = 0; distance < maxdistance and not terrain; distance++){
        raycast += dir;
        // side of box detection
        if (raycast.x < 5.0f){
            terrain = 1;
            avoidance += 0.5f*glm::vec3(1.0f, 0.0f, 0.0f);
            break;
        }
        if (raycast.x > bs->cave->sizex -5.0f){
            terrain = 1;
            avoidance -= 0.5f*glm::vec3(1.0f, 0.0f, 0.0f);
            break;
        }
        if (raycast.y < 5.0f){
            terrain = 1;
            avoidance += .5f*glm::vec3(0.0f, 1.0f, 0.0f);
            break;
        }
        if (raycast.y > bs->cave->sizey-5.0f){
            terrain = 1;
            avoidance -= .5f*glm::vec3(0.0f, 1.0f, 0.0f);
            break;
        }
        if (raycast.z < 5.0f){
            terrain = 1;
            avoidance += .5f*glm::vec3(0.0f, 0.0f, 1.0f);
            break;
        }
        if (raycast.z > bs->cave->sizez-5.0f){
            terrain = 1;
            avoidance -= .5f*glm::vec3(0.0f, 0.0f, 1.0f);
            break;
        }
        // terrain avoidance
        if (bs->cave->field(raycast.x, raycast.y, raycast.z) < bs->cave->level){
            terrain = 1; // im going to hit terrain!
            glm::vec3 grad = bs->cave->fieldGrad(raycast);
            avoidance += 1.0f*grad;
            break;
        }
    }

    if (count > 0){
        averagePos = averagePos/(float)count;
        cohesion = bs->cohesionForce*glm::normalize(averagePos-pos);
        averageVel = averageVel/(float)count;
        alignment = bs->alignmentForce*glm::normalize(averageVel-vel);
    }
    if (!terrain) vel += repelling + cohesion + alignment + avoidance + 2.0f*runaway;
    else vel += avoidance;
    vel = 6.0f*glm::normalize(vel);
    glm::vec3 newpos =  pos+vel*delta;
    glm::vec3 avoidvel = collisionResponse(pos, bs, &hit);
    if (hit){
        vel = 6.0f*avoidvel; 
        pos += vel*delta;
    } else {
        pos = newpos;
    }
}

void Killer::update(float delta){
    unsigned int count = 0; // number of neighbors
    float radius = 10.0f;
    glm::vec3 hunting = glm::vec3(0.0f); // hunting component of velocity
    glm::vec3 averagePos = glm::vec3(0.0f);
    
    for (int i=0; i<bs->boids.size(); i++){
        float distance = glm::distance(bs->boids[i]->pos, pos);

        // find close boids
        if (distance < radius){
            count++;
            //repelling component
            glm::vec3 toBoid = bs->boids[i]->pos - pos;
            // averaging position
            averagePos += bs->boids[i]->pos;
        } else continue;
    }

    // collision avoidance
    // boundaries
    glm::vec3 avoidance = glm::vec3(0.0f, 0.0f , 0.0f); // avoidance component
    int maxdistance = 5;
    glm::vec3 dir = glm::normalize(vel);
    glm::vec3 raycast = pos; // raycast starts at boid pos
    int terrain = 0; // is there a predicted terrain collision
    bool hit = false;
    
    //raycast loop increment distance each iteration
    for (int distance = 0; distance < maxdistance and not terrain; distance++){
        raycast += dir;
        // side of box detection
        if (raycast.x < 5.0f){
            terrain = 1;
            avoidance += 0.5f*glm::vec3(1.0f, 0.0f, 0.0f);
        }
        if (raycast.x > bs->cave->sizex -5.0f){
            terrain = 1;
            avoidance -= 0.5f*glm::vec3(1.0f, 0.0f, 0.0f);
        }
        if (raycast.y < 5.0f){
            terrain = 1;
            avoidance += .5f*glm::vec3(0.0f, 1.0f, 0.0f);
        }
        if (raycast.y > bs->cave->sizey-5.0f){
            terrain = 1;
            avoidance -= .5f*glm::vec3(0.0f, 1.0f, 0.0f);
        }
        if (raycast.z < 5.0f){
            terrain = 1;
            avoidance += .5f*glm::vec3(0.0f, 0.0f, 1.0f);
        }
        if (raycast.z > bs->cave->sizez-5.0f){
            terrain = 1;
            avoidance -= .5f*glm::vec3(0.0f, 0.0f, 1.0f);
        }
        // terrain avoidance
        if (bs->cave->field(raycast.x, raycast.y, raycast.z) < bs->cave->level){
            terrain = 1; // im going to hit terrain!
            glm::vec3 grad = bs->cave->fieldGrad(raycast);
            avoidance += 1.0f*grad;
        }
    }

    if (count > 0){
        averagePos = averagePos/(float)count;
        hunting = glm::normalize(averagePos - pos);
    }
    if (!terrain) vel += hunting + avoidance;
    else vel += avoidance;
    vel = 6.0f*glm::normalize(vel);
    glm::vec3 newpos =  pos+vel*delta;
    glm::vec3 avoidvel = collisionResponse(pos, bs, &hit);
    if (hit){
        vel = 6.0f*avoidvel; 
        pos += vel*delta;
    } else {
        pos = newpos;
    }
}

glm::vec3 randomCavePoint(Cave* cave){
    glm::vec3 point;
    do {
        point.x = cave->sizex*(randomf()*0.5+0.5);
        point.y = cave->sizey/2.0;
        point.z = cave->sizez*(randomf()*0.5+0.5);
    } while (cave->field(point.x, point.y, point.z) < cave->level);
    return point;
}

BoidSystem* Boid::bs;
BoidSystem* Killer::bs;
glm::mat4 Boid::scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.3f,0.5f,0.4f));
glm::mat4 Killer::scale = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f));

BoidSystem::BoidSystem(int amount, int amountkillers, Cave* cave){
    this->cave = cave;
    Boid::bs = this;
    Killer::bs = this;
    for (int i=0; i<amount; i++){
        Boid* boid = new Boid();
        boid->pos = randomCavePoint(cave);
        boid->vel = glm::sphericalRand(1.0f);
        boids.push_back(boid);
    }
    for (int i=0; i<amountkillers; i++){
        Killer* killer = new Killer();
        killer->pos = randomCavePoint(cave);
        killer->vel = glm::sphericalRand(1.0f);
        killers.push_back(killer);
    }
    verts =  {
2, 1.625, 1, 0.553452, 0.735419, 0.390961, 0, 0,
2, 2, 0.769231, 0.531098, 0.794179, 0.295321, 1, 0,
1.625, 2, 1, 0.465057, 0.806483, 0.365113, 0.5, 1,
2, 1.625, 1, 0.553452, 0.735419, 0.390961, 0, 0,
1.76923, 2, 2, 0.452244, 0.674538, 0.583501, 1, 0,
2, 1.76923, 2, 0.498196, 0.627943, 0.597904, 0.5, 1,
1.625, 2, 1, 0.465057, 0.806483, 0.365113, 0, 0,
1.76923, 2, 2, 0.452244, 0.674538, 0.583501, 0.5, 1,
2, 1.625, 1, 0.553452, 0.735419, 0.390961, 1, 0,
2, 1.76923, 2, 0.498196, 0.627943, 0.597904, 0, 0,
1.4, 2, 3, 0.346854, 0.559238, 0.752957, 1, 0,
2, 1.4, 3, 0.449462, 0.439744, 0.777566, 0.5, 1,
1.76923, 2, 2, 0.452244, 0.674538, 0.583501, 0, 0,
1.4, 2, 3, 0.346854, 0.559238, 0.752957, 0.5, 1,
2, 1.76923, 2, 0.498196, 0.627943, 0.597904, 1, 0,
2, 1.4, 3, 0.449462, 0.439744, 0.777566, 0, 0,
1.66667, 2, 4, 0.338714, 0.457164, 0.822359, 1, 0,
2, 1.66667, 4, 0.387752, 0.400327, 0.830293, 0.5, 1,
1.4, 2, 3, 0.346854, 0.559238, 0.752957, 0, 0,
1.66667, 2, 4, 0.338714, 0.457164, 0.822359, 0.5, 1,
2, 1.4, 3, 0.449462, 0.439744, 0.777566, 1, 0,
2, 1.66667, 4, 0.387752, 0.400327, 0.830293, 0, 0,
1.66667, 2, 4, 0.338714, 0.457164, 0.822359, 1, 0,
2, 2, 4.33333, 0.367425, 0.427693, 0.82588, 0.5, 1,
1.625, 2, 1, 0.465057, 0.806483, 0.365113, 0, 0,
2, 2, 0.769231, 0.531098, 0.794179, 0.295321, 1, 0,
2, 2.375, 1, 0.501467, 0.803654, 0.320423, 0.5, 1,
2, 2.23077, 2, 0.476418, 0.692015, 0.542348, 0, 0,
1.625, 2, 1, 0.465057, 0.806483, 0.365113, 1, 0,
2, 2.375, 1, 0.501467, 0.803654, 0.320423, 0.5, 1,
1.76923, 2, 2, 0.452244, 0.674538, 0.583501, 0, 0,
1.625, 2, 1, 0.465057, 0.806483, 0.365113, 0.5, 1,
2, 2.23077, 2, 0.476418, 0.692015, 0.542348, 1, 0,
2, 2.6, 3, 0.417642, 0.629508, 0.655205, 0, 0,
1.76923, 2, 2, 0.452244, 0.674538, 0.583501, 1, 0,
2, 2.23077, 2, 0.476418, 0.692015, 0.542348, 0.5, 1,
1.4, 2, 3, 0.346854, 0.559238, 0.752957, 0, 0,
1.76923, 2, 2, 0.452244, 0.674538, 0.583501, 0.5, 1,
2, 2.6, 3, 0.417642, 0.629508, 0.655205, 1, 0,
2, 2.33333, 4, 0.377769, 0.501562, 0.778284, 0, 0,
1.4, 2, 3, 0.346854, 0.559238, 0.752957, 1, 0,
2, 2.6, 3, 0.417642, 0.629508, 0.655205, 0.5, 1,
1.66667, 2, 4, 0.338714, 0.457164, 0.822359, 0, 0,
1.4, 2, 3, 0.346854, 0.559238, 0.752957, 0.5, 1,
2, 2.33333, 4, 0.377769, 0.501562, 0.778284, 1, 0,
2, 2, 4.33333, 0.367425, 0.427693, 0.82588, 0, 0,
1.66667, 2, 4, 0.338714, 0.457164, 0.822359, 1, 0,
2, 2.33333, 4, 0.377769, 0.501562, 0.778284, 0.5, 1,
2.375, 2, 1, 0.583181, 0.738605, 0.338176, 0, 0,
2, 2, 0.769231, 0.531098, 0.794179, 0.295321, 1, 0,
2, 1.625, 1, 0.553452, 0.735419, 0.390961, 0.5, 1,
2.375, 2, 1, 0.583181, 0.738605, 0.338176, 0, 0,
2, 1.76923, 2, 0.498196, 0.627943, 0.597904, 1, 0,
2.23077, 2, 2, 0.521021, 0.648616, 0.554829, 0.5, 1,
2, 1.625, 1, 0.553452, 0.735419, 0.390961, 0, 0,
2, 1.76923, 2, 0.498196, 0.627943, 0.597904, 0.5, 1,
2.375, 2, 1, 0.583181, 0.738605, 0.338176, 1, 0,
2.23077, 2, 2, 0.521021, 0.648616, 0.554829, 0, 0,
2, 1.4, 3, 0.449462, 0.439744, 0.777566, 1, 0,
2.6, 2, 3, 0.515137, 0.525437, 0.677162, 0.5, 1,
2, 1.76923, 2, 0.498196, 0.627943, 0.597904, 0, 0,
2, 1.4, 3, 0.449462, 0.439744, 0.777566, 0.5, 1,
2.23077, 2, 2, 0.521021, 0.648616, 0.554829, 1, 0,
2.6, 2, 3, 0.515137, 0.525437, 0.677162, 0, 0,
2, 1.66667, 4, 0.387752, 0.400327, 0.830293, 1, 0,
2.33333, 2, 4, 0.425662, 0.447509, 0.786478, 0.5, 1,
2, 1.4, 3, 0.449462, 0.439744, 0.777566, 0, 0,
2, 1.66667, 4, 0.387752, 0.400327, 0.830293, 0.5, 1,
2.6, 2, 3, 0.515137, 0.525437, 0.677162, 1, 0,
2, 2, 4.33333, 0.367425, 0.427693, 0.82588, 0, 0,
2.33333, 2, 4, 0.425662, 0.447509, 0.786478, 1, 0,
2, 1.66667, 4, 0.387752, 0.400327, 0.830293, 0.5, 1,
2.375, 2, 1, 0.583181, 0.738605, 0.338176, 0, 0,
2, 2.375, 1, 0.501467, 0.803654, 0.320423, 1, 0,
2, 2, 0.769231, 0.531098, 0.794179, 0.295321, 0.5, 1,
2.23077, 2, 2, 0.521021, 0.648616, 0.554829, 0, 0,
2, 2.375, 1, 0.501467, 0.803654, 0.320423, 1, 0,
2.375, 2, 1, 0.583181, 0.738605, 0.338176, 0.5, 1,
2, 2.23077, 2, 0.476418, 0.692015, 0.542348, 0, 0,
2, 2.375, 1, 0.501467, 0.803654, 0.320423, 0.5, 1,
2.23077, 2, 2, 0.521021, 0.648616, 0.554829, 1, 0,
2.6, 2, 3, 0.515137, 0.525437, 0.677162, 0, 0,
2, 2.23077, 2, 0.476418, 0.692015, 0.542348, 1, 0,
2.23077, 2, 2, 0.521021, 0.648616, 0.554829, 0.5, 1,
2, 2.6, 3, 0.417642, 0.629508, 0.655205, 0, 0,
2, 2.23077, 2, 0.476418, 0.692015, 0.542348, 0.5, 1,
2.6, 2, 3, 0.515137, 0.525437, 0.677162, 1, 0,
2.33333, 2, 4, 0.425662, 0.447509, 0.786478, 0, 0,
2, 2.6, 3, 0.417642, 0.629508, 0.655205, 1, 0,
2.6, 2, 3, 0.515137, 0.525437, 0.677162, 0.5, 1,
2, 2.33333, 4, 0.377769, 0.501562, 0.778284, 0, 0,
2, 2.6, 3, 0.417642, 0.629508, 0.655205, 0.5, 1,
2.33333, 2, 4, 0.425662, 0.447509, 0.786478, 1, 0,
2, 2, 4.33333, 0.367425, 0.427693, 0.82588, 0, 0,
2, 2.33333, 4, 0.377769, 0.501562, 0.778284, 1, 0,
2.33333, 2, 4, 0.425662, 0.447509, 0.786478, 0.5, 1 };   
    neighboorRadius = 5.0f;
    repellingForce = -0.1f;
    cohesionForce = 0.1f;
    alignmentForce = 0.05f;
}

void BoidSystem::update(float delta){
    for (int i=0; i<boids.size(); i++){
        boids[i]->update(delta);
    }
    for (int i=0; i<killers.size(); i++){
        killers[i]->update(delta);
    }
}
