#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <random>
#include <numeric> // for std::iota
#include <memory>  // for std::unique_ptr

// Constants for customization
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;
const int TERRAIN_WIDTH = 100;       // Number of points along X axis
const int TERRAIN_HEIGHT = 100;      // Number of points along Z axis
const float NOISE_SCALE = 4.0f;      // Higher = more rugged terrain
const float NOISE_HEIGHT = 100.0f;   // Height multiplier
const int NOISE_OCTAVES = 3;         // More octaves = more detail
const float VIEW_ANGLE = 70.0f;      // Elevation angle in degrees
const float VIEW_DISTANCE = 1000.0f;  // Perspective distance
const int GRID_SPACING = 1;          // Higher = sparser grid lines
const sf::Color BACKGROUND_COLOR = sf::Color::Black;
const sf::Color WIREFRAME_COLOR = sf::Color::White;

class PerlinNoise {
private:
    std::vector<int> p;
    
public:
    PerlinNoise(unsigned int seed) {
        // Initialize the permutation vector with random values based on seed
        p.resize(256);
        std::iota(p.begin(), p.end(), 0);
        std::default_random_engine engine(seed);
        std::shuffle(p.begin(), p.end(), engine);
        
        // Duplicate the permutation vector
        p.insert(p.end(), p.begin(), p.end());
    }
    
    double noise(double x, double y, double z = 0) {
        // Find the unit cube that contains the point
        int X = (int)floor(x) & 255;
        int Y = (int)floor(y) & 255;
        int Z = (int)floor(z) & 255;
        
        // Find relative x, y, z of point in cube
        x -= floor(x);
        y -= floor(y);
        z -= floor(z);
        
        // Compute fade curves for each of x, y, z
        double u = fade(x);
        double v = fade(y);
        double w = fade(z);
        
        // Hash coordinates of the 8 cube corners
        int A = p[X] + Y;
        int AA = p[A] + Z;
        int AB = p[A + 1] + Z;
        int B = p[X + 1] + Y;
        int BA = p[B] + Z;
        int BB = p[B + 1] + Z;
        
        // Add blended results from 8 corners of cube
        double res = lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),
                                          grad(p[BA], x-1, y, z)),
                                    lerp(u, grad(p[AB], x, y-1, z),
                                          grad(p[BB], x-1, y-1, z))),
                              lerp(v, lerp(u, grad(p[AA+1], x, y, z-1),
                                          grad(p[BA+1], x-1, y, z-1)),
                                    lerp(u, grad(p[AB+1], x, y-1, z-1),
                                          grad(p[BB+1], x-1, y-1, z-1))));
        return (res + 1.0) / 2.0;
    }
    
private:
    double fade(double t) {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }
    
    double lerp(double t, double a, double b) {
        return a + t * (b - a);
    }
    
    double grad(int hash, double x, double y, double z) {
        int h = hash & 15;
        double u = h < 8 ? x : y;
        double v = h < 4 ? y : h == 12 || h == 14 ? x : z;
        return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }
};

class WireframeTerrain {
private:
    sf::VertexArray vertices;
    std::vector<std::vector<float>> heightMap;
    std::unique_ptr<PerlinNoise> pn;
    
public:
    WireframeTerrain() {
        // Generate random seed
        std::random_device rd;
        unsigned int seed = rd();
        
        // Initialize Perlin noise with random seed
        pn = std::make_unique<PerlinNoise>(seed);
        
        heightMap.resize(TERRAIN_WIDTH, std::vector<float>(TERRAIN_HEIGHT));
        generateTerrain();
        updateVertices();
    }
    
    void generateTerrain() {
        for (int i = 0; i < TERRAIN_WIDTH; ++i) {
            for (int j = 0; j < TERRAIN_HEIGHT; ++j) {
                double x = (double)i / TERRAIN_WIDTH * NOISE_SCALE;
                double y = (double)j / TERRAIN_HEIGHT * NOISE_SCALE;
                
                // Generate fractal noise by adding multiple octaves
                double value = 0.0;
                double amplitude = 1.0;
                double frequency = 1.0;
                
                for (int octave = 0; octave < NOISE_OCTAVES; octave++) {
                    value += pn->noise(x * frequency, y * frequency) * amplitude;
                    amplitude *= 0.5;
                    frequency *= 2.0;
                }
                
                heightMap[i][j] = value * NOISE_HEIGHT;
            }
        }
    }
    
    void updateVertices() {
        vertices.clear();
        vertices.setPrimitiveType(sf::Lines);
        
        // Convert view angle to radians
        const float viewAngleRad = VIEW_ANGLE * 3.14159265f / 180.0f;
        const float sinAngle = sin(viewAngleRad);
        const float cosAngle = cos(viewAngleRad);
        
        const float centerX = TERRAIN_WIDTH / 2.0f;
        const float centerY = TERRAIN_HEIGHT / 2.0f;
        const float screenCenterX = WINDOW_WIDTH / 2.0f;
        const float screenCenterY = WINDOW_HEIGHT / 2.0f;
        const float scaleFactor = std::min(WINDOW_WIDTH, WINDOW_HEIGHT) / std::max(TERRAIN_WIDTH, TERRAIN_HEIGHT);
        
        // Draw horizontal lines
        for (int i = 0; i < TERRAIN_WIDTH; i += GRID_SPACING) {
            for (int j = 0; j < TERRAIN_HEIGHT - 1; ++j) {
                // Project 3D point to 2D with perspective
                float x1 = (i - centerX) * scaleFactor;
                float y1 = heightMap[i][j];
                float z1 = (j - centerY) * scaleFactor;
                
                float x2 = (i - centerX) * scaleFactor;
                float y2 = heightMap[i][j+1];
                float z2 = ((j+1) - centerY) * scaleFactor;
                
                // Apply perspective projection
                float scale1 = VIEW_DISTANCE / (VIEW_DISTANCE - z1 * sinAngle + y1 * cosAngle);
                float scale2 = VIEW_DISTANCE / (VIEW_DISTANCE - z2 * sinAngle + y2 * cosAngle);
                
                float projX1 = x1 * scale1 + screenCenterX;
                float projY1 = (y1 * sinAngle + z1 * cosAngle) * scale1 + screenCenterY;
                
                float projX2 = x2 * scale2 + screenCenterX;
                float projY2 = (y2 * sinAngle + z2 * cosAngle) * scale2 + screenCenterY;
                
                vertices.append(sf::Vertex(sf::Vector2f(projX1, projY1), WIREFRAME_COLOR));
                vertices.append(sf::Vertex(sf::Vector2f(projX2, projY2), WIREFRAME_COLOR));
            }
        }
        
        // Draw vertical lines
        for (int j = 0; j < TERRAIN_HEIGHT; j += GRID_SPACING) {
            for (int i = 0; i < TERRAIN_WIDTH - 1; ++i) {
                // Project 3D point to 2D with perspective
                float x1 = (i - centerX) * scaleFactor;
                float y1 = heightMap[i][j];
                float z1 = (j - centerY) * scaleFactor;
                
                float x2 = ((i+1) - centerX) * scaleFactor;
                float y2 = heightMap[i+1][j];
                float z2 = (j - centerY) * scaleFactor;
                
                // Apply perspective projection
                float scale1 = VIEW_DISTANCE / (VIEW_DISTANCE - z1 * sinAngle + y1 * cosAngle);
                float scale2 = VIEW_DISTANCE / (VIEW_DISTANCE - z2 * sinAngle + y2 * cosAngle);
                
                float projX1 = x1 * scale1 + screenCenterX;
                float projY1 = (y1 * sinAngle + z1 * cosAngle) * scale1 + screenCenterY;
                
                float projX2 = x2 * scale2 + screenCenterX;
                float projY2 = (y2 * sinAngle + z2 * cosAngle) * scale2 + screenCenterY;
                
                vertices.append(sf::Vertex(sf::Vector2f(projX1, projY1), WIREFRAME_COLOR));
                vertices.append(sf::Vertex(sf::Vector2f(projX2, projY2), WIREFRAME_COLOR));
            }
        }
    }
    
    void draw(sf::RenderWindow& window) {
        window.draw(vertices);
    }
};

int main() {
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Perlin Noise Wireframe Mountains");
    window.setFramerateLimit(60);
    
    WireframeTerrain terrain;
    
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        
        window.clear(BACKGROUND_COLOR);
        terrain.draw(window);
        window.display();
    }
    
    return 0;
}