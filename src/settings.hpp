#pragma once

#include <iostream>
#include "raylib.h"
#include "raymath.h"
#include <cmath>
#include <cstdint>
#include <vector>

const int WIDTH =  800;
const int HEIGHT = 800;

//math shit
const double PM = 1.00727647;
const double NM = 1.008665;

const double CONSTE = 2.71828;

const int EC_RES = 8; //electron cloud resolution, DONT MAKE IT LESS THAN 8 IF YOU WANT FAST!!!!
const double UPP = 0.01; //units per pixel
const double BOHR_RAD = 0.529;
const double BOND_THRESHOLD = 0.01;

void draw_grid(int size, Color color, int spacing, Vector2 position) {
    Vector2 curr_pos = position;
    for (int y = -(size + spacing); y < (WIDTH / (size + spacing)) + (size + spacing); y++) {
        for (int x = -(size + spacing); x < (WIDTH / (size + spacing)) + (size + spacing); x++) {
            DrawRectangle(curr_pos.x, curr_pos.y, size, size, color);
            curr_pos.x += size + spacing;
        }
        curr_pos.x = position.x;
        curr_pos.y += size + spacing;
    }
}