class Vector3ff
{
public:
    float x, y, z;
    Vector3ff(float x = 0, float y = 0, float z = 0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    operator Vector3f() const { return Vector3f(this->x, this->y, this->z); }
    Vector3ff operator+(Vector3ff v2)
    {
        return Vector3ff(this->x + v2.x, this->y + v2.y, this->z + v2.z);
    }
    Vector3ff operator-(Vector3ff v2)
    {
        return Vector3ff(this->x - v2.x, this->y - v2.y, this->z - v2.z);
    }
    Vector3ff operator*(float a)
    {
        return Vector3ff(this->x * a, this->y * a, this->z * a);
    }
    Vector3ff operator/(float a)
    {
        return Vector3ff(this->x / a, this->y / a, this->z / a);
    }
    Vector3ff &operator+=(Vector3ff a)
    {
        this->x += a.x;
        this->y += a.y;
        this->z += a.z;
        return *this;
    }
    Vector3ff &operator=(Vector3ff a)
    {
        this->x = a.x;
        this->y = a.y;
        this->z = a.z;
        return *this;
    }
    Vector3ff operator-()
    {
        return Vector3ff(-this->x, -this->y, -this->z);
    }
    friend ostream &operator<<(ostream &out, const Vector3ff &vec)
    {
        out << vec.x << ", " << vec.y << ", " << vec.z;
        return out;
    }
    float magnitude()
    {
        return sqrtf(this->x * this->x + this->y * this->y + this->z * this->z);
    }
    float dot(Vector3ff b)
    {
        return this->x * b.x + this->y * b.y + this->z * b.z;
    }
    Vector3ff reflect(Vector3ff normal)
    {
        return (*this) - normal * 2.0f * (*this).dot(normal);
    }
    Vector3ff cross(Vector3ff b)
    {
        return Vector3ff(
            -this->z * b.y + this->y * b.z,
            this->z * b.x - this->x * b.z,
            -this->y * b.x + this->x * b.y);
    }
    Vector3ff normalize()
    {
        return Vector3ff(this->x, this->y, this->z) / this->magnitude();
    }
};
Vector3ff cross(Vector3ff a, Vector3ff b)
{
    return Vector3ff(
        -a.z * b.y + a.y * b.z,
        a.z * b.x - a.x * b.z,
        -a.y * b.x + a.x * b.y);
}
Vector3ff to3ff(Vector3f a)
{
    return Vector3ff(a.x, a.y, a.z);
}
float magnitude(Vector2f a)
{
    return sqrtf(a.x * a.x + a.y * a.y);
}
float dot(Vector2f a, Vector2f b)
{
    return a.x * b.x + a.y * b.y;
}
Vector2f normalize(Vector2f a)
{
    float l = magnitude(a);
    if (l == 0)
        return Vector2f(0, 0);
    return Vector2f(a.x, a.y) / l;
}
Vector2f Lerp(Vector2f a, Vector2f b, float t)
{
    return a + (b - a) * t;
}
class Complex
{
public:
    float a, b;
    Complex(float a = 0, float b = 0)
    {
        this->a = a;
        this->b = b;
    }
};
class NullDivide
{
private:
    int getFloatsCount(float a)
    {
        float b = a;
        while (fmodf(b, 10.0f) != 0)
        {
            b *= 10;
        }
        return b;
    }

public:
    float a;
    NullDivide(float a = 0)
    {
        this->a = a;
    }
    // friend const NullDivide &operator+(const float &left, const NullDivide &right);
    NullDivide operator+(float a)
    {
        float n = this->getFloatsCount(a);
        float c = powf(10.0f, n);
        float b = a * c;

        return NullDivide(c * a);
    }
};

float random()
{
    return rand() / (float)RAND_MAX;
}
float srandom(unsigned int seed)
{
    float v = 0;
    srand(seed);
    v = random();
    srand(time(NULL));
    return v;
}
Vector3ff MatrixMult3x3(Vector3ff v, float m[3][3])
{
    return Vector3ff(
        v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
        v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
        v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2]);
}

class Camera
{
public:
    Vector2f position;
    float scale, scrollSensivity = 0.1f, moveSensivity = 1;
    Camera(Vector2f StartPosition = Vector2f(0, 0), float Scale = 1)
    {
        position = StartPosition;
        scale = Scale;
    }
    void onMove(Vector2f deltaPosition)
    {
        position += (Vector2f)deltaPosition * this->moveSensivity / this->scale;
    }
    Vector2f toViewport(Vector2f globalPosition)
    {
        return (globalPosition + position) * scale + Vector2f(WIDTH, HEIGHT) / 2.0f;
    }
    Vector2f toWorld(Vector2f viewportPosition)
    {
        return (viewportPosition - Vector2f(WIDTH, HEIGHT) / 2.0f) / scale - position;
    }
};
class Camera3d
{
private:
    // Shader depthBufferDrawShader;

public:
    Vector3ff position, rotation;
    Vector3ff forward, right, up;
    float fov;
    // RenderTexture depthBuffer;
    Camera3d(Vector3ff StartPosition = Vector3ff(0, 0, 0), float Fov = 1)
    {
        this->position = StartPosition;
        this->fov = Fov;

        this->rotation = Vector3ff(0, 0, 0);
        this->forward = Vector3ff(0, 0, 1);
        this->right = Vector3ff(1, 0, 0);
        this->up = Vector3ff(0, 1, 0);

        // this->depthBuffer.create(WIDTH, HEIGHT);
        // this->depthBufferDrawShader.loadFromFile("DepthBufferDrawerFrag.glsl", Shader::Type::Fragment);
        // this->depthBufferDrawShader.setUniform("u_resolution", Glsl::Vec2(WIDTH, HEIGHT));
    }
    void Rotate(Vector3ff rotation)
    {
        this->rotation = rotation;
        this->forward = Vector3ff(0, 0, 1);
        this->right = Vector3ff(1, 0, 0);
        this->up = Vector3ff(0, 1, 0);

        // rotation.x += PI;
        // this->right = Vector3ff(cosf(rotation.x + PI / 2), 0, sinf(rotation.x + PI / 2));
        // this->forward = Vector3ff(cosf(rotation.y) * cosf(rotation.x), sinf(rotation.y), sinf(rotation.x) * cosf(rotation.y));
        // this->up = Vector3ff(cosf(rotation.y + PI / 2) * cosf(rotation.x), sinf(rotation.y + PI / 2), sinf(rotation.x) * cosf(rotation.y + PI / 2));

        this->forward = Vector3ff(cosf(rotation.y) * sinf(rotation.x), sinf(rotation.y), cosf(rotation.x) * cosf(rotation.y));
        this->right = Vector3ff(sinf(rotation.x + PI / 2.0f), 0.0f, cosf(rotation.x + PI / 2.0f));
        this->up = Vector3ff(cosf(rotation.y + PI / 2.0f) * sinf(rotation.x), sinf(rotation.y + PI / 2.0f), cosf(rotation.x) * cosf(rotation.y + PI / 2.0f));
        // cout << this->forward << endl;
        // cout << this->forward.dot(this->right) << ", " << this->right.dot(this->up) << ", " << this->up.dot(this->forward) << endl;
        return;

        float xRotMat[3][3] = {
            {1, 0, 0},
            {0, cosf(rotation.x), -sinf(rotation.x)},
            {0, sinf(rotation.x), cosf(rotation.x)}};
        float yRotMat[3][3] = {
            {cosf(rotation.y), 0, sinf(rotation.y)},
            {0, 1, 0},
            {-sinf(rotation.y), 0, cosf(rotation.y)}};
        float zRotMat[3][3] = {
            {cosf(rotation.z), -sinf(rotation.z), 0},
            {sinf(rotation.z), cosf(rotation.z), 0},
            {0, 0, 1}};

        this->forward = MatrixMult3x3(this->forward, xRotMat);
        this->forward = MatrixMult3x3(this->forward, yRotMat);
        this->forward = MatrixMult3x3(this->forward, zRotMat);

        this->right = MatrixMult3x3(this->right, xRotMat);
        this->right = MatrixMult3x3(this->right, yRotMat);
        this->right = MatrixMult3x3(this->right, zRotMat);

        this->up = MatrixMult3x3(this->up, xRotMat);
        this->up = MatrixMult3x3(this->up, yRotMat);
        this->up = MatrixMult3x3(this->up, zRotMat);
    }
    Vector2f ProjectToCanvas(Vector3ff point)
    {
        float size = min(WIDTH, HEIGHT);

        Vector3ff delta = point - this->position;
        Vector3ff nDelta = delta.normalize();
        // float t = this->forward.dot(this->forward) / nDelta.dot(this->forward);
        if (nDelta.dot(this->forward) < 0)
            return Vector2f(-1, -1);
        Vector3ff v = (nDelta / nDelta.dot(this->forward) - this->forward) * size * this->fov;
        return Vector2f(v.dot(this->right) + WIDTH / 2, -v.dot(this->up) + HEIGHT / 2);
    }
    bool isOnScreen(Vector2f p)
    {
        if (p.x > WIDTH || p.y > HEIGHT)
        {
            return false;
        }
        if (p.x < 0 || p.y < 0)
        {
            return false;
        }
        return true;
    }
    Vector2f ProjectToCanvas(Vector3f point)
    {
        float size = min(WIDTH, HEIGHT);

        Vector3ff delta = Vector3ff(point.x - this->position.x, point.y - this->position.y, point.z - this->position.z);
        Vector3ff nDelta = delta.normalize();
        // float t = this->forward.dot(this->forward) / nDelta.dot(this->forward);
        Vector3ff v = (nDelta / nDelta.dot(this->forward) - this->forward) * size * this->fov;
        return Vector2f(v.dot(this->right) + WIDTH / 2, -v.dot(this->up) + HEIGHT / 2);
    }
    bool DrawToDepthBuffer(Vector3ff p1, Vector3ff p2, Vector3ff p3)
    {
        // return false;
        // VertexArray v(Triangles, 3);
        // v[0] = Vertex(this->ProjectToCanvas(p1), Vector2f(0, 0));
        // v[1] = Vertex(this->ProjectToCanvas(p2), Vector2f(1, 0));
        // v[2] = Vertex(this->ProjectToCanvas(p3), Vector2f(0, 1));

        // this->depthBufferDrawShader.setUniform("cameraPosition", Glsl::Vec3(this->position.x, this->position.y, this->position.z));

        // Vector3ff transformativeP1 = this->right * p1.x + this->up * p1.y + this->forward * p1.z;
        // Vector3ff transformativeP2 = this->right * p2.x + this->up * p2.y + this->forward * p2.z;
        // Vector3ff transformativeP3 = this->right * p3.x + this->up * p3.y + this->forward * p3.z;
        // this->depthBufferDrawShader.setUniform("p1", Glsl::Vec3(transformativeP1.x, transformativeP1.y, transformativeP1.z));
        // this->depthBufferDrawShader.setUniform("p2", Glsl::Vec3(transformativeP2.x, transformativeP2.y, transformativeP2.z));
        // this->depthBufferDrawShader.setUniform("p3", Glsl::Vec3(transformativeP3.x, transformativeP3.y, transformativeP3.z));

        // Vector2f lp1 = this->ProjectToCanvas(p1);
        // Vector2f lp2 = this->ProjectToCanvas(p2);
        // Vector2f lp3 = this->ProjectToCanvas(p3);

        // this->depthBufferDrawShader.setUniform("lp1", Glsl::Vec2(lp1.x, lp1.y));
        // this->depthBufferDrawShader.setUniform("lp2", Glsl::Vec2(lp2.x, lp2.y));
        // this->depthBufferDrawShader.setUniform("lp3", Glsl::Vec2(lp3.x, lp3.y));

        // RenderStates state;
        // state.shader = &(this->depthBufferDrawShader);
        // this->depthBuffer.draw(v, state);

        return true;
    }
};
Vector2f VecWithoutZ(Vector3ff v)
{
    return Vector2f(v.x, v.y);
}
class Particle
{
private:
public:
    float mass, radius, restitution;
    Vector3ff position, velocity;
    CircleShape drawer;
    Color color;
    bool isDisabled = false;
    Particle(Vector3ff startPosition = Vector3ff(0, 0, 0), Vector3ff startVelocity = Vector3ff(0, 0, 0), float Mass = 1, float Radius = 5, Color Clr = Color::Cyan, float Restitution = 1.6)
    {
        position = startPosition;
        velocity = startVelocity;
        mass = Mass;
        radius = Radius;
        color = Clr;
        restitution = Restitution;

        drawer = CircleShape(radius);

        drawer.setOrigin(radius, radius);
        drawer.setFillColor(color);
        drawer.setPosition(VecWithoutZ(position));
    }
    void UpdatePosition()
    {
        if (!this->isDisabled)
            position += velocity;
        // drawer.setRadius(radius * defaultCam.scale);
        // drawer.setOrigin(radius * defaultCam.scale, radius * defaultCam.scale);
        // drawer.setPosition(defaultCam.toViewport(position));
    }
};

class Mesh
{
private:
    vector<Vector3ff> trianglesMidPoints;
    void calucalteTrianglesMidPoints()
    {
        this->trianglesMidPoints.clear();
        for (int i = 0; i < this->triangles.size(); i += 3)
        {
            Vector3ff pos = (this->vertices[this->triangles[i]] + this->vertices[this->triangles[i + 1]] + this->vertices[this->triangles[i + 2]]) / 3.0f;
            this->trianglesMidPoints.push_back(pos);
        }
    }
    void sortTriangles(Camera3d *CamToSort)
    {
        this->calucalteTrianglesMidPoints();
        bool direction = false;
        while (true)
        {
            int swapsCount = 0;
            if (direction)
            {
                for (int i = 1; i < trianglesMidPoints.size(); i++)
                {
                    float l1 = (trianglesMidPoints[i] - CamToSort->position).magnitude();
                    float l2 = ((trianglesMidPoints[i - 1] - CamToSort->position).magnitude());

                    if (l2 < l1)
                    {
                        swap(trianglesMidPoints[i], trianglesMidPoints[i - 1]);

                        swap(this->triangles[i * 3], this->triangles[(i - 1) * 3]);
                        swap(this->triangles[i * 3 + 1], this->triangles[(i - 1) * 3 + 1]);
                        swap(this->triangles[i * 3 + 2], this->triangles[(i - 1) * 3 + 2]);

                        swap(this->normals[i], this->normals[i - 1]);

                        swapsCount++;
                    }
                }
            }
            else
            {
                for (int i = trianglesMidPoints.size() - 1 - 1; i >= 0; i--)
                {
                    float l1 = (trianglesMidPoints[i] - CamToSort->position).magnitude();
                    float l2 = ((trianglesMidPoints[i + 1] - CamToSort->position).magnitude());

                    if (l2 > l1)
                    {
                        swap(trianglesMidPoints[i], trianglesMidPoints[i + 1]);

                        swap(this->triangles[i * 3], this->triangles[(i + 1) * 3]);
                        swap(this->triangles[i * 3 + 1], this->triangles[(i + 1) * 3 + 1]);
                        swap(this->triangles[i * 3 + 2], this->triangles[(i + 1) * 3 + 2]);

                        swap(this->normals[i], this->normals[i + 1]);

                        swapsCount++;
                    }
                }
            }
            if (swapsCount == 0)
                break;
            direction = !direction;
        }
    }
    VertexArray generateVertexArray(Camera3d *cameraToRender, Vector3ff sunDirection)
    {
        int trianglesCount = 0;
        if (this->triangles.size() % 3 == 0)
            int trianglesCount = this->triangles.size() / 3;
        else
        {
            cout << "Triangles count not dividable 3: " << trianglesCount << endl;
            return VertexArray();
        }
        VertexArray V(Triangles, trianglesCount);

        for (int i = 0, p = 0; i < this->triangles.size(); i += 3, p++)
        {
            Color color = Color::White;
            if (this->normals.size() > 0)
            {
                // float reflectedV =
                // float d = -sunDirection.dot(this->normals[this->triangles[i]]);
                // d = (d < 0.4) ? 0.4f : d;

                float E = this->normals[p].dot(sunDirection);
                if (E < 0.1)
                    E = 0.1;
                // E = 1.0f;
                color = Color(E * 255, E * 255, E * 255, 255);
                // float V = (float)i / this->triangles.size();
                // color = Color(V * 255, 0, 255 - V * 255, 255);
            }
            Vertex v1 = Vertex(cameraToRender->ProjectToCanvas(vertices[this->triangles[i]]), color);
            Vertex v2 = Vertex(cameraToRender->ProjectToCanvas(vertices[this->triangles[i + 1]]), color);
            Vertex v3 = Vertex(cameraToRender->ProjectToCanvas(vertices[this->triangles[i + 2]]), color);

            V.append(v1);
            V.append(v2);
            V.append(v3);

            cameraToRender->DrawToDepthBuffer(vertices[this->triangles[i]], vertices[this->triangles[i + 1]], vertices[this->triangles[i + 2]]);
        }

        return V;
    }
    VertexArray generateGridVertexArray(Camera3d *cameraToRender, Color clr)
    {
        if (this->normals.size() == 0)
            this->recalculateNormals();

        int trianglesCount = 0;
        if (this->triangles.size() % 3 == 0)
            int trianglesCount = this->triangles.size() / 3;
        else
        {
            cout << "Triangles count not dividable 3: " << trianglesCount << endl;
            return VertexArray();
        }
        VertexArray V(Lines, trianglesCount);
        for (int i = 0; i < this->triangles.size(); i += 3)
        {
            Color color = clr;
            Vertex v1 = Vertex(cameraToRender->ProjectToCanvas(vertices[this->triangles[i]]), color);
            Vertex v2 = Vertex(cameraToRender->ProjectToCanvas(vertices[this->triangles[i + 1]]), color);
            Vertex v3 = Vertex(cameraToRender->ProjectToCanvas(vertices[this->triangles[i + 2]]), color);

            V.append(v1);
            V.append(v2);

            V.append(v2);
            V.append(v3);

            V.append(v3);
            V.append(v1);
        }
        return V;
    }
    VertexArray generatePointsVertexArray(Camera3d *cameraToRender)
    {

        VertexArray V(Points, this->vertices.size());
        for (int i = 0; i < this->vertices.size(); i++)
        {
            Color color = Color::White;
            Vertex v1 = Vertex(cameraToRender->ProjectToCanvas(vertices[i]), color);
            V.append(v1);
        }
        return V;
    }

public:
    vector<Vector3ff> vertices = vector<Vector3ff>();
    vector<int> triangles = vector<int>();
    vector<Vector2f> uv = vector<Vector2f>();
    vector<Vector3ff> normals = vector<Vector3ff>();

    VertexArray vertexArray;

    Mesh()
    {
        // vertices
    }

    VertexArray &getVertexArray(Camera3d *camToRender, Vector3ff sunDirection = Vector3ff(0, -1, 0))
    {
        this->sortTriangles(camToRender);
        this->vertexArray = generateVertexArray(camToRender, sunDirection);
        return (this->vertexArray);
    }
    VertexArray &getGridVertexArray(Camera3d *camToRender, Color clr = Color::White)
    {
        // this->sortTriangles(camToRender);
        this->vertexArray = generateGridVertexArray(camToRender, clr);
        return (this->vertexArray);
    }
    VertexArray &getPointsVertexArray(Camera3d *camToRender, Vector3ff sunDirection = Vector3ff(0, -1, 0))
    {
        // this->sortTriangles(camToRender);
        this->vertexArray = generatePointsVertexArray(camToRender);
        return (this->vertexArray);
    }
    bool recalculateNormals()
    {
        normals.clear();
        for (int i = 0; i < this->triangles.size(); i += 3)
        {
            Vector3ff delta1 = this->vertices[this->triangles[i + 1]] - this->vertices[this->triangles[i]];
            Vector3ff delta2 = this->vertices[this->triangles[i + 2]] - this->vertices[this->triangles[i]];
            normals.push_back(cross(delta1, delta2));
        }
        return true;
    }
};
void zeroAxisLines(RenderWindow *window, Camera3d *camera)
{
    Vector3f p1 = Vector3f(-1, 0, 0);
    Vector3f p2 = Vector3f(0, -1, 0);
    Vector3f p3 = Vector3f(0, 0, -1);
    Vector2f vp1 = camera->ProjectToCanvas(p1);
    Vector2f vp2 = camera->ProjectToCanvas(p2);
    Vector2f vp3 = camera->ProjectToCanvas(p3);
    Vector2f mvp1 = camera->ProjectToCanvas(-p1);
    Vector2f mvp2 = camera->ProjectToCanvas(-p2);
    Vector2f mvp3 = camera->ProjectToCanvas(-p3);

    VertexArray Lines = VertexArray(PrimitiveType::Lines, 6);
    Lines.append(Vertex(vp1, Color::Red));
    Lines.append(Vertex(mvp1, Color::Red));
    Lines.append(Vertex(vp2, Color::Green));
    Lines.append(Vertex(mvp2, Color::Green));
    Lines.append(Vertex(vp3, Color::Blue));
    Lines.append(Vertex(mvp3, Color::Blue));

    window->draw(Lines);
}

class Plane
{
private:
    VertexArray triangles;

public:
    vector<vector<Vector3ff>> gridPoints;
    Vector2i dotsCount;
    Vector2f size;

    Vector3ff position, sunDirection;
    Color color;
    float lightness;
    Plane(Vector3ff startPosition = Vector3ff(0, 0, 0), Vector2i dotsCount = Vector2i(0, 0), Vector2f size = Vector2f(0, 0))
    {
        this->position = startPosition;
        this->dotsCount = dotsCount;
        this->size = size;

        this->sunDirection = Vector3ff(-0.1, -1, -0.1);
        this->lightness = 1.0f;

        this->color = Color(200, 20, 20);

        // VertexArray triangles(Triangles, this->dotsCount.x * this->dotsCount.y * 3 * 2);
    }
    void GenerateGrid()
    {
        gridPoints.clear();
        for (int y = 0; y < this->dotsCount.y; y++)
        {
            vector<Vector3ff> gridRow;
            for (int x = 0; x < this->dotsCount.x; x++)
            {
                // float lx = ((float)x / (float)this->dotsCount.x - 0.5f) * 10;
                // float ly = ((float)y / (float)this->dotsCount.y - 0.5f) * 10;
                // Vector3ff l = Vector3ff(lx, ly, 0);
                // float n = sqrtf(powf((l - Vector3ff(3, 2, 0)).magnitude(), 2) - powf((Vector3ff(6, 9, 0) - Vector3ff(3, 2, 0)).normalize().dot(l - Vector3ff(3, 2, 0)), 2));
                // float v = logf(cosf(lx) + sinf(ly));

                // float v = sin(lx + ly);
                // if (x > 0 && y > 0)
                //     v = (gridRow[x - 1].y + gridPoints[y - 1][x].y) / 2 + (-1 + srandom(x * (y + this->dotsCount.x)) * 2) * 0.02;
                gridRow.push_back(Vector3ff((float)x / (float)this->dotsCount.x * this->size.x, 0, (float)y / (float)this->dotsCount.y * this->size.y));
            }
            gridPoints.push_back(gridRow);
        }
    }
    void SetGridByArray(vector<vector<float>> arr)
    {
        for (int i = 0; i < this->dotsCount.x; i++)
        {
            for (int j = 0; j < this->dotsCount.y; j++)
            {
                this->gridPoints[i][j].y = arr[i][j];
            }
        }
    }
    void Draw(RenderWindow *window, Camera3d *RenderCamera)
    {
        int verticesCount = this->dotsCount.x * this->dotsCount.y * 3 * 2;
        VertexArray triangles(Triangles, this->dotsCount.x * this->dotsCount.y * 3 * 2);
        for (int y = 1; y < this->gridPoints.size(); y++)
        {
            for (int x = 1; x < this->gridPoints[y].size(); x++)
            {
                Vector2f screenPosition1 = RenderCamera->ProjectToCanvas(this->position + this->gridPoints[y][x]);
                Vector2f screenPosition2 = RenderCamera->ProjectToCanvas(this->position + this->gridPoints[y][x - 1]);
                Vector2f screenPosition3 = RenderCamera->ProjectToCanvas(this->position + this->gridPoints[y - 1][x]);
                Vector2f screenPosition4 = RenderCamera->ProjectToCanvas(this->position + this->gridPoints[y - 1][x - 1]);

                Vector3ff viewDirection = ((this->position + this->gridPoints[y][x - 1]) - RenderCamera->position).normalize();

                Vector3ff normal1 = (this->gridPoints[y - 1][x] - this->gridPoints[y][x]).cross(this->gridPoints[y][x - 1] - this->gridPoints[y][x]).normalize();
                Vector3ff normal2 = (this->gridPoints[y][x - 1] - this->gridPoints[y - 1][x - 1]).cross(this->gridPoints[y - 1][x] - this->gridPoints[y - 1][x - 1]).normalize();

                float v1 = normal1.dot(-sunDirection.normalize());
                float v2 = normal2.dot(-sunDirection.normalize());

                v1 = powf(v1, lightness);
                v2 = powf(v2, lightness);

                if (v1 < 0)
                    v1 = 0;
                if (v2 < 0)
                    v2 = 0;

                Color clr1 = Color(this->color.r * v1, this->color.g * v1, this->color.b * v1, 255);
                Color clr2 = Color(this->color.r * v2, this->color.g * v2, this->color.b * v2, 255);

                // clr1 = Color(normal1.x * 255, normal1.y * 255, normal1.z * 255, 255);

                triangles.append(Vertex(screenPosition1, clr1));
                triangles.append(Vertex(screenPosition2, clr1));
                triangles.append(Vertex(screenPosition3, clr1));

                triangles.append(Vertex(screenPosition4, clr2));
                triangles.append(Vertex(screenPosition2, clr2));
                triangles.append(Vertex(screenPosition3, clr2));

                // cout << screenPosition1.x << ", " << screenPosition1.y << endl;
            }
        }
        window->draw(triangles);
    }
};