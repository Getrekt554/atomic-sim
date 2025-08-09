#include "settings.hpp"

// one unit one angstrom

float get_ideal_bond_angle(int Z, int neighbor_count) {
    if (Z == 8 && neighbor_count == 2) return 104.5f * 3.14159f / 180.0f;
    if (Z == 6) {
        if (neighbor_count == 2) return 180.0f * 3.14159f / 180.0f;
        if (neighbor_count == 3) return 120.0f * 3.14159f / 180.0f;
        if (neighbor_count == 4) return 109.5f * 3.14159f / 180.0f;
    }
    if (Z == 7 && neighbor_count == 3) return 107.0f * 3.14159f / 180.0f;
    return 180.0f * 3.14159f / 180.0f;
}

class Nucleus {
    public:
        double amu;
        int Z;
        int N;
        double radius;
        Vector2 pos;
        Vector2 vel;

        Nucleus(int N, int Z, Vector2 pos) {
            this->amu = (Z * PM) + (N * NM);
            this->Z = Z;
            this->N = N;
            this->radius = estimate_radius(amu);
            this->pos = pos;
            this->vel = {0, 0};
        }
        Nucleus() {

        }
    private:
        double estimate_radius(double atomic_mass) {
            double cbrt_atomic_mass = cbrt(atomic_mass);
            double fm_rad = cbrt_atomic_mass * 1.2; //1.2 fm
            double AU_rad = fm_rad / 100000; // oh yea, we convertin to angstroms in my part of town son
            return AU_rad;
        }
};

class Atom {
public:
    bool is_full_shell() const {

        return (shell_1s == 2 && shell_2s == 2 && shell_2p == 6) ||
               (shell_1s == 2 && shell_2s == 0 && shell_2p == 0 && nucleus.Z == 2);
    }
    Nucleus nucleus;
    int E;
    int shell_1s = 0;
    int shell_2s = 0;
    int shell_2p = 0;
    int valence = 0;
    int max_bonds = 0;
    int current_bonds = 0;
    int charge = 0;
    Atom(Vector2 pos, int Z, int N, int E) {
        this->nucleus = Nucleus(N, Z, pos);
        this->E = E;
        this->charge = Z - E; 
        switch (Z) {
            case 1: valence = 1; max_bonds = 1; break; // H
            case 2: valence = 0; max_bonds = 0; break; // He
            case 3: valence = 1; max_bonds = 1; break; // Li
            case 4: valence = 2; max_bonds = 2; break; // Be
            case 5: valence = 3; max_bonds = 3; break; // B
            case 6: valence = 4; max_bonds = 4; break; // C
            case 7: valence = 3; max_bonds = 3; break; // N
            case 8: valence = 2; max_bonds = 2; break; // O
            case 9: valence = 1; max_bonds = 1; break; // F
            case 10: valence = 0; max_bonds = 0; break; // Ne
            default: valence = 1; max_bonds = 1; break;
        }
    }
        void tick() {
            const float maxVel = 1.0f;
            float speed = Vector2Length(nucleus.vel);
            if (speed > maxVel) {
                nucleus.vel = Vector2Scale(Vector2Normalize(nucleus.vel), maxVel);
            }
            nucleus.pos = Vector2Add(nucleus.vel, nucleus.pos);
            nucleus.vel = Vector2Scale(nucleus.vel, 0.6f);
            form_clouds();
        }
       double p(Vector2 coord) const {
            Vector2 relative = coord - nucleus.pos;
            double r = sqrt(relative.x * relative.x + relative.y * relative.y);
            double theta = atan2(relative.y, relative.x);
            double Z = (double)nucleus.Z;

            double s_1s = 0;
            double s_2s = 0;
            double s_2p = 0;
            if (E >= 1) {
                s_1s = pow(Z, 3) * exp(-2 * Z * r);
            }
            if (E >= 3) {
                double term = (1 - Z * r / 2);
                s_2s = pow(Z, 3) * term * term * exp(-Z * r);
            }
            if (E >= 5) {
                double angular = cos(theta) * cos(theta);
                s_2p = pow(Z, 3) * r * r * angular * exp(-Z * r);
            }

            return (s_1s + s_2s + s_2p) / static_cast<double>(E);
        }
 
        double get_effective_radius() const {
            int Z = nucleus.Z;
            
            if (Z == 1) return 0.37; // H, pm
            if (Z == 2) return 0.32; // He
            if (Z == 3) return 1.28; // Li
            if (Z == 4) return 0.96; // Be
            if (Z == 5) return 0.84; // B
            if (Z == 6) return 0.76; // C
            if (Z == 7) return 0.71; // N
            if (Z == 8) return 0.66; // O
            if (Z == 9) return 0.57; // F
            if (Z == 10) return 0.58; // Ne

            return 1.0;
        }
        double get_electronegativity() const {
            int Z = nucleus.Z;
            
            switch (Z) {
                case 1: return 2.20;
                case 2: return 0.00; // noble gas, doesn't bond
                case 3: return 0.98;
                case 4: return 1.57;
                case 5: return 2.04;
                case 6: return 2.55;
                case 7: return 3.04;
                case 8: return 3.44;
                case 9: return 3.98;
                case 10: return 0.00;
                default: return 1.0;
            }
}



    private:
        void form_clouds() {
            int temp_e_amount = E;
            shell_1s = 0;
            shell_2s = 0;
            shell_2p = 0;
            while (temp_e_amount != 0) {
                temp_e_amount--;
                if (shell_1s < 2) {
                    shell_1s++;
                    continue;
                }
                else if (shell_2s < 2) {
                    shell_2s++;
                    continue;
                }
                else if (shell_2p < 6) {
                    shell_2p++;
                    continue;
                }
            }
        }
};

void tick_atoms(std::vector<Atom>& Atoms) {
    for (auto& atom : Atoms) {
        atom.tick();
    }
}
void draw_atoms(const std::vector<Atom>& Atoms) {
    for (float y = 0; y < HEIGHT; y += EC_RES) {
        for (float x = 0; x < WIDTH; x += EC_RES) {
            Vector2 sample = {x * UPP, y * UPP};
            double total_density = 0.0;

            for (const Atom& atom : Atoms) {
                total_density += atom.p(sample);
            }

            // Scale and clamp
            double scaled = total_density;
            if (scaled > 1.0) scaled = 1.0;

            uint8_t brightness = static_cast<uint8_t>(scaled * 255);

            DrawRectangle(static_cast<int>(x), static_cast<int>(y), EC_RES, EC_RES,
                          Color{brightness, brightness, brightness, brightness});
        }
    }
}



double get_equilibrium_distance(const Atom& a1, const Atom& a2) {
    double r1 = a1.get_effective_radius();
    double r2 = a2.get_effective_radius();
    return r1 + r2;
}

double get_epsilon(const Atom& a1, const Atom& a2) {
    double en1 = a1.get_electronegativity();
    double en2 = a2.get_electronegativity();
    double average_en = (en1 + en2) / 2.0;
    return average_en * 0.004;
}



double compute_overlap(const Atom& a1, const Atom& a2) {
    double total = 0.0;
    for (float y = 0; y < HEIGHT; y += EC_RES) {
        for (float x = 0; x < WIDTH; x += EC_RES) {
            Vector2 world = {x * UPP, y * UPP};

            double p1 = a1.p(world);
            double p2 = a2.p(world);
            total += p1 * p2;
        }
    }
    return total;
}


double get_covalent_radius(int atomic_number) {
    switch (atomic_number) {
        case 1: return 0.31f; // H
        case 2: return 0.37f; // He
        case 3: return 1.28f; // Li
        case 4: return 0.96f; // Be
        case 5: return 0.84f; // B
        case 6: return 0.76f; // C
        case 7: return 0.71f; // N
        case 8: return 0.66f; // O
        case 9: return 0.57f; // F
        case 10: return 0.58f; // Ne
        default: return 1.0f; // Default for unknown elements
    }
}

double calculate_ideal_bond_length(Atom* a1, Atom* a2) {
    double r1 = get_covalent_radius(a1->nucleus.Z);
    double r2 = get_covalent_radius(a2->nucleus.Z);
    return r1 + r2;
}
class Bond {
public:
    int atomAIndex;
    int atomBIndex;
    double bondLength;
    double bondStrength;
    double bondDamping;

    Bond(int a, int b, double length) : atomAIndex(a), atomBIndex(b), bondLength(length), bondStrength(100), bondDamping(3) {}

    void apply_bond_force(std::vector<Atom>& atoms) {
        Atom& a = atoms[atomAIndex];
        Atom& b = atoms[atomBIndex];

        Vector2 delta = Vector2Subtract(b.nucleus.pos, a.nucleus.pos);
        double dist = Vector2Length(delta);
        if (dist == 0) return;

        double diff = dist - bondLength;
        Vector2 dir = Vector2Scale(delta, 1.0 / dist);

        double total_mass = a.nucleus.amu + b.nucleus.amu;
        double ratio_a = b.nucleus.amu / total_mass;
        double ratio_b = a.nucleus.amu / total_mass;

        a.nucleus.pos = Vector2Add(a.nucleus.pos, Vector2Scale(dir, diff * ratio_a));
        b.nucleus.pos = Vector2Subtract(b.nucleus.pos, Vector2Scale(dir, diff * ratio_b));

        Vector2 rel_vel = Vector2Subtract(b.nucleus.vel, a.nucleus.vel);
        double rel_vel_along_dir = rel_vel.x * dir.x + rel_vel.y * dir.y;
        Vector2 correction_vel = Vector2Scale(dir, rel_vel_along_dir * 0.5);
        a.nucleus.vel = Vector2Add(a.nucleus.vel, correction_vel);
        b.nucleus.vel = Vector2Subtract(b.nucleus.vel, correction_vel);
    }

};

void apply_forces_from_clouds(std::vector<Atom>& atoms, std::vector<Bond>& bonds){
    for (int i = 0; i < atoms.size(); i++) {
        for (int j = i + 1; j < atoms.size(); j++) {
            double overlap = compute_overlap(atoms[i], atoms[j]);

            if (overlap < 0.00001) continue;
            bool bonded = false;
            for (const Bond& bond : bonds) {
                if ((bond.atomAIndex == i && bond.atomBIndex == j) ||
                    (bond.atomAIndex == j && bond.atomBIndex == i)) {
                    bonded = true;
                    break;
                }
            }
            if (bonded) continue;

            Atom& a1 = atoms[i];
            Atom& a2 = atoms[j];


            Vector2 delta = Vector2Subtract(a2.nucleus.pos, a1.nucleus.pos);
            double dist = Vector2Length(delta);
            dist = fmax(dist, 0.1);
            if (dist == 0) continue;

            Vector2 dir = Vector2Normalize(delta);

            double r0 = get_equilibrium_distance(a1, a2);
            double epsilon = get_epsilon(a1, a2);

            double lj_force = 24 * epsilon * (2 * pow(r0 / dist, 13) - pow(r0 / dist, 7)) / dist;
            double max_lj_force = 0.2;
            if (lj_force > max_lj_force) lj_force = max_lj_force;
            if (lj_force < -max_lj_force) lj_force = -max_lj_force;
            Vector2 force = Vector2Scale(dir, lj_force);
            a1.nucleus.vel = Vector2Add(a1.nucleus.vel, force);
            a2.nucleus.vel = Vector2Subtract(a2.nucleus.vel, force);
        }
    }
}


int main() {
    auto can_form_bond = [&](int i, int j, const std::vector<Atom>& atoms, const std::vector<Bond>& bonds) -> bool {
        if (i == j) return false;
        const Atom& a = atoms[i];
        const Atom& b = atoms[j];
        if (a.is_full_shell() || b.is_full_shell()) return false;
        if (a.valence <= 0 || b.valence <= 0) return false;
        if (a.current_bonds >= a.max_bonds || b.current_bonds >= b.max_bonds) return false;
        int bond_count = 0;
        for (const Bond& bond : bonds) {
            if ((bond.atomAIndex == i && bond.atomBIndex == j) || (bond.atomAIndex == j && bond.atomBIndex == i))
                bond_count++;
        }
        if ((a.nucleus.Z == 6 && b.nucleus.Z == 7) || (a.nucleus.Z == 7 && b.nucleus.Z == 6)) {
            if (bond_count >= 3) return false;
        } else if (a.nucleus.Z == 1 && b.nucleus.Z == 1) {
            if (bond_count >= 1) return false;
        } else {
            if (bond_count >= 2) return false;
        }
        if (a.charge != 0 && b.charge != 0 && (a.charge * b.charge > 0)) return false; // both same sign
        return true;
    };
    InitWindow(WIDTH, HEIGHT, "atoms");

    SetTargetFPS(60);

    float centerX = WIDTH * 0.5f * UPP;
    float centerY = HEIGHT * 0.5f * UPP;
    float bondLen = 1.16f; 
    std::vector<Atom> atoms = {
        Atom(Vector2{centerX, centerY}, 6, 6, 6),
    };

    std::vector<Bond> bonds;



    while (WindowShouldClose() == false) {
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            Vector2 mouse = GetMousePosition();
            Vector2 world = { mouse.x * UPP, mouse.y * UPP };
            atoms.emplace_back(world, 7, 7, 7);
        }
        bool all_angles_stable = true;
        for (int i = 0; i < atoms.size(); ++i) {
            if (!(atoms[i].nucleus.Z == 6 || atoms[i].nucleus.Z == 8)) {
                std::vector<int> neighbors;
                for (const Bond& bond : bonds) {
                    if (bond.atomAIndex == i) neighbors.push_back(bond.atomBIndex);
                    else if (bond.atomBIndex == i) neighbors.push_back(bond.atomAIndex);
                }
                if (neighbors.size() == 1) {
                    int j = neighbors[0];
                    // Only apply if both atoms have 1 bond (diatomic)
                    if (atoms[i].current_bonds == 1 && atoms[j].current_bonds == 1) {
                        double m1 = atoms[i].nucleus.amu;
                        double m2 = atoms[j].nucleus.amu;
                        double total_mass = m1 + m2;
                        Vector2 com = {
                            (atoms[i].nucleus.pos.x * m1 + atoms[j].nucleus.pos.x * m2) / total_mass,
                            (atoms[i].nucleus.pos.y * m1 + atoms[j].nucleus.pos.y * m2) / total_mass
                        };
                        Vector2 r1 = {atoms[i].nucleus.pos.x - com.x, atoms[i].nucleus.pos.y - com.y};
                        Vector2 r2 = {atoms[j].nucleus.pos.x - com.x, atoms[j].nucleus.pos.y - com.y};
                        Vector2 v1 = atoms[i].nucleus.vel;
                        Vector2 v2 = atoms[j].nucleus.vel;
                        double Lz = m1 * (r1.x * v1.y - r1.y * v1.x) + m2 * (r2.x * v2.y - r2.y * v2.x);
                        double I = m1 * (r1.x * r1.x + r1.y * r1.y) + m2 * (r2.x * r2.x + r2.y * r2.y);
                        if (I > 0.0) {
                            double omega = Lz / I;
                            double damping = 0.2; // Local damping for diatomics
                            double damp_omega = omega * damping;
                            Vector2 tangential1 = {-r1.y, r1.x};
                            Vector2 tangential2 = {-r2.y, r2.x};
                            atoms[i].nucleus.vel.x -= damp_omega * tangential1.x;
                            atoms[i].nucleus.vel.y -= damp_omega * tangential1.y;
                            atoms[j].nucleus.vel.x -= damp_omega * tangential2.x;
                            atoms[j].nucleus.vel.y -= damp_omega * tangential2.y;
                        }
                    }
                }
                continue;
            }
            std::vector<int> neighbors;
            for (const Bond& bond : bonds) {
                if (bond.atomAIndex == i) neighbors.push_back(bond.atomBIndex);
                else if (bond.atomBIndex == i) neighbors.push_back(bond.atomAIndex);
            }
            if (neighbors.size() < 2) continue;
            for (size_t a = 0; a < neighbors.size(); ++a) {
                for (size_t b = a + 1; b < neighbors.size(); ++b) {
                    int j = neighbors[a];
                    int k = neighbors[b];
                    Vector2 v1 = Vector2Normalize(Vector2Subtract(atoms[j].nucleus.pos, atoms[i].nucleus.pos));
                    Vector2 v2 = Vector2Normalize(Vector2Subtract(atoms[k].nucleus.pos, atoms[i].nucleus.pos));
                    float dot = v1.x * v2.x + v1.y * v2.y;
                    dot = fmaxf(-1.0f, fminf(1.0f, dot));
                    float angle = acosf(dot);
                    float target_angle = get_ideal_bond_angle(atoms[i].nucleus.Z, (int)neighbors.size());
                    float angle_error = angle - target_angle;
                    if (fabs(angle_error) >= 0.03f) all_angles_stable = false;
                }
            }
        }
        tick_atoms(atoms);
        // Only skip force/bond updates if stable AND at least 3 individual bonds exist
        if (!all_angles_stable || bonds.size() < 3) {
            apply_forces_from_clouds(atoms, bonds);


            for (int i = 0; i < atoms.size(); i++) {
                for (int j = i + 1; j < atoms.size(); j++) {
                    double dx = atoms[j].nucleus.pos.x - atoms[i].nucleus.pos.x;
                    double dy = atoms[j].nucleus.pos.y - atoms[i].nucleus.pos.y;
                    double dist = sqrt(dx * dx + dy * dy);
                    double idealDist = calculate_ideal_bond_length(&atoms[i], &atoms[j]);

                    if (dist < idealDist * 1.3 && can_form_bond(i, j, atoms, bonds)) {
                        bonds.emplace_back(i, j, idealDist);
                        atoms[i].current_bonds++;
                        atoms[j].current_bonds++;
                        atoms[i].valence--;
                        atoms[j].valence--;
                        printf("Bond formed between atom %d (Z=%d, charge=%d) and atom %d (Z=%d, charge=%d) at dist=%.2f\n",
                            i, atoms[i].nucleus.Z, atoms[i].charge,
                            j, atoms[j].nucleus.Z, atoms[j].charge,
                            dist);
                    }
                }
            }

            for (Bond& bond : bonds) {
                bond.apply_bond_force(atoms);
            }
        }



        // --- Angle correction for atoms with 2+ bonds ---
        // (Already handled above, remove duplicate block)

        // --- True angular momentum damping (whole molecule, always active) ---
        // Compute center of mass
        Vector2 com = {0, 0};
        double total_mass = 0.0;
        for (const Atom& atom : atoms) {
            com.x += atom.nucleus.pos.x * atom.nucleus.amu;
            com.y += atom.nucleus.pos.y * atom.nucleus.amu;
            total_mass += atom.nucleus.amu;
        }
        if (total_mass > 0.0) {
            com.x /= total_mass;
            com.y /= total_mass;
        }
        // Compute total angular momentum (about center of mass)
        double Lz = 0.0;
        double I = 0.0;
        for (const Atom& atom : atoms) {
            Vector2 r = {atom.nucleus.pos.x - com.x, atom.nucleus.pos.y - com.y};
            Vector2 v = atom.nucleus.vel;
            // 2D cross product (z-component): r.x * v.y - r.y * v.x
            Lz += atom.nucleus.amu * (r.x * v.y - r.y * v.x);
            I += atom.nucleus.amu * (r.x * r.x + r.y * r.y);
        }
        if (I > 0.0) {
            double omega = Lz / I; // angular velocity (scalar, z-axis)
            double damping = 0.1; // Lower damping to prevent overcorrection/flapping
            double damp_omega = omega * damping;
            // Subtract angular velocity from each atom's velocity
            for (auto& atom : atoms) {
                Vector2 r = {atom.nucleus.pos.x - com.x, atom.nucleus.pos.y - com.y};
                // Perpendicular vector (for 2D rotation): (-r.y, r.x)
                Vector2 tangential = {-r.y, r.x};
                atom.nucleus.vel.x -= damp_omega * tangential.x;
                atom.nucleus.vel.y -= damp_omega * tangential.y;
            }
        }

        // --- Lock internal geometry of any fully bonded diatomic fragment (e.g., CN-) ---
        if (all_angles_stable) {
            for (const Bond& bond : bonds) {
                int i = bond.atomAIndex;
                int j = bond.atomBIndex;
                // Both atoms must have only one bond, and both must be at max bonds (fully satisfied diatomic)
                if (atoms[i].current_bonds == 1 && atoms[j].current_bonds == 1 &&
                    atoms[i].current_bonds == atoms[i].max_bonds && atoms[j].current_bonds == atoms[j].max_bonds) {
                    // Lock their bond length and orientation, but preserve center-of-mass
                    double m1 = atoms[i].nucleus.amu;
                    double m2 = atoms[j].nucleus.amu;
                    double total_mass = m1 + m2;
                    Vector2 com = {
                        (atoms[i].nucleus.pos.x * m1 + atoms[j].nucleus.pos.x * m2) / total_mass,
                        (atoms[i].nucleus.pos.y * m1 + atoms[j].nucleus.pos.y * m2) / total_mass
                    };
                    // Desired bond vector (current direction)
                    Vector2 bond_vec = Vector2Subtract(atoms[j].nucleus.pos, atoms[i].nucleus.pos);
                    double bond_len = sqrt(bond_vec.x * bond_vec.x + bond_vec.y * bond_vec.y);
                    // Target bond length (from bond)
                    double target_len = bond.bondLength;
                    // Normalize bond_vec
                    if (bond_len > 0.0) {
                        bond_vec.x /= bond_len;
                        bond_vec.y /= bond_len;
                    } else {
                        bond_vec = {1, 0}; // arbitrary direction if degenerate
                    }
                    // Place i and j at correct positions along bond_vec, centered at com
                    Vector2 pi = {com.x - bond_vec.x * target_len * m2 / total_mass,
                                  com.y - bond_vec.y * target_len * m2 / total_mass};
                    Vector2 pj = {com.x + bond_vec.x * target_len * m1 / total_mass,
                                  com.y + bond_vec.y * target_len * m1 / total_mass};
                    atoms[i].nucleus.pos = pi;
                    atoms[j].nucleus.pos = pj;
                }
            }
            // Only anchor and zero velocities if the system is exactly two atoms (standalone diatomic)
            if (atoms.size() == 2) {
                // Zero center-of-mass velocity
                Vector2 v_cm = {0, 0};
                double total_mass = 0.0;
                for (const Atom& atom : atoms) {
                    v_cm.x += atom.nucleus.vel.x * atom.nucleus.amu;
                    v_cm.y += atom.nucleus.vel.y * atom.nucleus.amu;
                    total_mass += atom.nucleus.amu;
                }
                if (total_mass > 0.0) {
                    v_cm.x /= total_mass;
                    v_cm.y /= total_mass;
                    for (auto& atom : atoms) {
                        atom.nucleus.vel.x -= v_cm.x;
                        atom.nucleus.vel.y -= v_cm.y;
                    }
                }
                // Anchor center-of-mass position to screen center
                Vector2 com = {0, 0};
                total_mass = 0.0;
                for (const Atom& atom : atoms) {
                    com.x += atom.nucleus.pos.x * atom.nucleus.amu;
                    com.y += atom.nucleus.pos.y * atom.nucleus.amu;
                    total_mass += atom.nucleus.amu;
                }
                if (total_mass > 0.0) {
                    com.x /= total_mass;
                    com.y /= total_mass;
                    float centerX = WIDTH * 0.5f * UPP;
                    float centerY = HEIGHT * 0.5f * UPP;
                    Vector2 shift = {centerX - com.x, centerY - com.y};
                    for (auto& atom : atoms) {
                        atom.nucleus.pos.x += shift.x;
                        atom.nucleus.pos.y += shift.y;
                    }
                }
                // Forcibly zero all atom velocities to guarantee no drift
                for (auto& atom : atoms) {
                    atom.nucleus.vel = {0, 0};
                }
            }
        }

        BeginDrawing();
        DrawFPS(0, 0);
        ClearBackground({0, 0, 0});
        draw_atoms(atoms);
        for (const Bond& bond : bonds) {
            Vector2 posA = atoms[bond.atomAIndex].nucleus.pos;
            Vector2 posB = atoms[bond.atomBIndex].nucleus.pos;
            // Count how many bonds exist between this pair
            int bond_count = 0;
            for (const Bond& b2 : bonds) {
                if ((b2.atomAIndex == bond.atomAIndex && b2.atomBIndex == bond.atomBIndex) ||
                    (b2.atomAIndex == bond.atomBIndex && b2.atomBIndex == bond.atomAIndex)) {
                    bond_count++;
                }
            }
            // Only draw each bond once (avoid duplicate lines)
            if (bond.atomAIndex < bond.atomBIndex) {
                Vector2 screenA = { posA.x / UPP, posA.y / UPP };
                Vector2 screenB = { posB.x / UPP, posB.y / UPP };
                float dx = screenB.x - screenA.x;
                float dy = screenB.y - screenA.y;
                float len = sqrtf(dx*dx + dy*dy);
                float offx = -dy / len * 5.0f; // 5 pixels offset
                float offy = dx / len * 5.0f;
                if (bond_count == 1) {
                    DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x, (int)screenB.y, RED);
                } else if (bond_count == 2) {
                    // Draw two parallel lines for double bond
                    DrawLine((int)(screenA.x + offx), (int)(screenA.y + offy), (int)(screenB.x + offx), (int)(screenB.y + offy), RED);
                    DrawLine((int)(screenA.x - offx), (int)(screenA.y - offy), (int)(screenB.x - offx), (int)(screenB.y - offy), RED);
                } else if (bond_count == 3) {
                    // Draw three parallel lines for triple bond
                    DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x, (int)screenB.y, RED);
                    DrawLine((int)(screenA.x + offx), (int)(screenA.y + offy), (int)(screenB.x + offx), (int)(screenB.y + offy), RED);
                    DrawLine((int)(screenA.x - offx), (int)(screenA.y - offy), (int)(screenB.x - offx), (int)(screenB.y - offy), RED);
                }
            }
        }
        EndDrawing();
    }

    CloseWindow();

    return 0;
}

