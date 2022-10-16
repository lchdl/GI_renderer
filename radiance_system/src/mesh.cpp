#include "mesh.h"
#include <stdlib.h>
#include <string.h>

Mesh::OBJFaceFormat Mesh::getOBJFaceFormat(const char * file)
{
    FILE* fp = fopen(file, "r");
    if (fp == NULL) return OBJFaceFormat::InvalidFace;

    char buf[256];

    OBJFaceFormat fmt = OBJFaceFormat::InvalidFace;

    /* start parsing */
    while (txtGetWord(fp, buf, 256, " #\n\r\t/", '#'))
    {
        if (strcmp(buf, "f") == 0)
        {
            txtGetWord(fp, buf, 256, " #\n\r\t", '#');
            int n = 0;
            for (int i = 0; i < strlen(buf); i++) if (buf[i] == '/') n++;
            if (n == 0) { fmt = OBJFaceFormat::V; break; }
            else if (n == 1) { fmt = OBJFaceFormat::V_Vt; break; }
            else if (n == 2) {
                if (strstr(buf, "//") != NULL) { fmt = OBJFaceFormat::V_Vn; break; }
                else { fmt = OBJFaceFormat::V_Vt_Vn; break; }
            }
            else { fmt = OBJFaceFormat::InvalidFace; break; }
        }
    }

    fclose(fp);
    return fmt;
}

Mesh::Mesh()
{
    usage = GL_STATIC_DRAW;
}

Mesh::Mesh(const Mesh & that)
{
    this->p = that.p;
    this->n = that.n;
    this->t = that.t;
    this->f = that.f;
    this->usage = that.usage;

    if (!isNull()) {
        /* model is loaded, now we need to upload data to GPU */
        /* vertex format is (X,Y,Z, U,V) */
        vbuf.gen(); /* initialize vertex fBuffer */
        packVertices();
        vbuf.loadData((void*)(packedVertices.data()), sizeof(GL_VERTEX_XYZ_RGB_UV) * packedVertices.size(), this->usage);
    }
}

Mesh::~Mesh()
{
    unload();
    vbuf.unload();
}

Mesh & Mesh::operator=(const Mesh & that)
{
    if (this == &that)
        return (*this);
    this->unload();

    this->p = that.p;
    this->n = that.n;
    this->t = that.t;
    this->f = that.f;
    this->usage = that.usage;

    if (!isNull()) {
        /* model is loaded, now we need to upload data to GPU */
        /* vertex format is (X,Y,Z, U,V) */
        vbuf.gen(); /* initialize vertex fBuffer */
        packVertices();
        vbuf.loadData((void*)(packedVertices.data()), sizeof(GL_VERTEX_XYZ_RGB_UV) * packedVertices.size(), this->usage);
    }

    return (*this);
}

void Mesh::packVertices()
{
    packedVertices.clear();

    int iv, it, in;
    GL_VERTEX_XYZ_RGB_UV v;
    for (int i = 0; i < this->f.size(); i++) {
        iv = f[i].x;
        it = f[i].y;
        in = f[i].z;
        v.position[0] = float(this->p[iv].x);
        v.position[1] = float(this->p[iv].y);
        v.position[2] = float(this->p[iv].z);
        v.normal[0] = float(this->n[in].x);
        v.normal[1] = float(this->n[in].y);
        v.normal[2] = float(this->n[in].z);
        v.texcoord[0] = float(this->t[it].x);
        v.texcoord[1] = float(this->t[it].y);
        packedVertices.append(v);
    }
}

bool Mesh::loadOBJ(const char * file, GLenum usage)
{
    unload();

    FILE* fp = fopen(file, "r");
    if (fp == NULL) {
        printf("error, cannot open file '%s' for loading.", file);
        return false;
    }

    OBJFaceFormat fmt = getOBJFaceFormat(file);
    if (!(fmt == OBJFaceFormat::V_Vt_Vn || fmt == OBJFaceFormat::V_Vn)) {
        printf("error, invalid OBJ face format. "
            "Only support face format \"V_Vt_Vn\" and \"V_Vn\".");
        fclose(fp);
        return false;
    }

    char buf[256];

    VEC3 v;
    INT2 z2;
    INT3 z3;


    /* start parsing */
    while (txtGetWord(fp, buf, 256, " #\n\r\t/", '#'))
    {
        if (strcmp(buf, "v") == 0) { /* load vector */
            if (!txtGetVec3(fp, &v)) return false;
            this->p.append(v);
        }
        else if (strcmp(buf, "vn") == 0) {
            if (!txtGetVec3(fp, &v)) return false;
            this->n.append(v);
        }
        else if (strcmp(buf, "vt") == 0) {
            if (!txtGetVec3(fp, &v)) return false;
            this->t.append(v);
        }
        else if (strcmp(buf, "f") == 0) {
            /* index starts with 0 but obj index starts with 1,     */
            /* so we need to decrease each vertex index by 1 to get */
            /* the correct offset. */
            if (fmt == V_Vt_Vn) {
                if (!txtGetInt3(fp, &z3)) return false;
                z3.x--; z3.y--; z3.z--; /* f v/t/n */
                this->f.append(z3);
                if (!txtGetInt3(fp, &z3)) return false;
                z3.x--; z3.y--; z3.z--;
                this->f.append(z3);
                if (!txtGetInt3(fp, &z3)) return false;
                z3.x--; z3.y--; z3.z--;
                this->f.append(z3);
            }
            else if (fmt == V_Vn) {
                if (!txtGetInt2(fp, &z2)) return false;
                z2.x--; z2.y--; /* f v//n */
                this->f.append(INT3(z2.x, 0, z2.y));
                if (!txtGetInt2(fp, &z2)) return false;
                z2.x--; z2.y--;
                this->f.append(INT3(z2.x, 0, z2.y));
                if (!txtGetInt2(fp, &z2)) return false;
                z2.x--; z2.y--;
                this->f.append(INT3(z2.x, 0, z2.y));
                if (this->t.size() == 0)
                    this->t.append(VEC3()); /* add a dummy coordinate */
            }
        }
    }

    fclose(fp);

    /* we dont need to allocate vertex buffer for this mesh if it doesnt */
    /* contain any valid data to draw on screen */
    if (!isNull()) {
        /* model is loaded, now we need to upload data to GPU */
        /* vertex format is (X,Y,Z, U,V) */
        vbuf.gen(); /* initialize vertex buffer */
        packVertices();
        this->usage = usage;
        vbuf.loadData((void*)(packedVertices.data()),
            sizeof(GL_VERTEX_XYZ_RGB_UV) * packedVertices.size(), usage);
    }

    return true;
}

void Mesh::draw()
{
    if (!isNull()) {
        /* DO NOT FORGET TO SET MODEL, VIEW AND PROJECTION MATRICES */
        /* BEFORE CALLING THIS FUNCTION !!! */
        vbuf.draw(GL_TRIANGLES, 0, packedVertices.size());
    }
}

void Mesh::unload()
{
    p.clear();
    n.clear();
    t.clear();
    f.clear();
    vbuf.unload(); /* unload all buffers allocated in GPU in case something bad happens */
    packedVertices.clear();
}

void Mesh::upload()
{
    vbuf.gen();
    packVertices();
    vbuf.loadData((void*)(packedVertices.data()), sizeof(GL_VERTEX_XYZ_RGB_UV) * packedVertices.size(), usage);
    //vbuf.updateData((void*)(packedVertices.data()), 0, sizeof(VERTEX_XYZ_RGB_UV) * packedVertices.size());
}

const Array<VEC3>& Mesh::getPositionArray() const { return p; }
VEC3 Mesh::getVertex(const int & index) const
{
    if (index < 0 || index >= p.size())
        return VEC3::zeros();
    return p[index];
}
bool Mesh::setVertex(const int & index, const VEC3 & v)
{
    if (index < 0 || index >= p.size())
        return false;
    p[index] = v;
    return true;
}
const Array<VEC3>& Mesh::getNormalArray() const { return n; }
const Array<VEC3>& Mesh::getTexCoordArray() const { return t; }
const Array<INT3>& Mesh::getFaceIndexArray() const { return f; }

int Mesh::bytes()
{
    int raw = sizeof(VEC3) * (p.size() + n.size() + t.size()) + sizeof(INT3) * f.size();
    int packed = sizeof(GL_VERTEX_XYZ_RGB_UV) * packedVertices.size();
    return raw + packed;
}

bool Mesh::saveOBJ(const char * file)
{
    FILE* fp = fopen(file, "w");
    if (fp == NULL) {
        return false;
    }
    char buffer[256];
    const char* fixedString[] = {
        "# vertex positions\n",
        "# vertex normals\n",
        "# vertex texture coordinates\n",
        "# face indices (v/vt/vn)\n"
    };

    fwrite(fixedString[0], 1, strlen(fixedString[0]), fp);
    for (int i = 0; i < p.size(); i++) {
        VEC3 v = p[i];
        sprintf(buffer, "v %.6f %.6f %.6f\n", v.x, v.y, v.z);
        fwrite(buffer, sizeof(char), strlen(buffer), fp);
    }

    fwrite(fixedString[1], 1, strlen(fixedString[1]), fp);
    for (int i = 0; i < n.size(); i++) {
        VEC3 v = n[i];
        sprintf(buffer, "vn %.6f %.6f %.6f\n", v.x, v.y, v.z);
        fwrite(buffer, sizeof(char), strlen(buffer), fp);
    }

    fwrite(fixedString[2], 1, strlen(fixedString[2]), fp);
    for (int i = 0; i < t.size(); i++) {
        VEC3 v = t[i];
        sprintf(buffer, "vt %.6f %.6f %.6f\n", v.x, v.y, v.z);
        fwrite(buffer, sizeof(char), strlen(buffer), fp);
    }

    fwrite(fixedString[3], 1, strlen(fixedString[3]), fp);
    int nFaces = f.size() / 3;
    for (int i = 0; i < nFaces; i++) {
        INT3 v1 = f[i * 3 + 0];
        INT3 v2 = f[i * 3 + 1];
        INT3 v3 = f[i * 3 + 2];
        sprintf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
            v1.x + 1, v1.y + 1, v1.z + 1,
            v2.x + 1, v2.y + 1, v2.z + 1,
            v3.x + 1, v3.y + 1, v3.z + 1); /* face index starts with 1 */
        fwrite(buffer, sizeof(char), strlen(buffer), fp);
    }

    fclose(fp);
    return true;
}

void Mesh::rotate(QUAT & q)
{
    MAT3x3 mRot = quatToMatrix(normalize(q));
    /* rotate both vertex positions and normals */
    for (int i = 0; i < this->p.size(); i++) {
        this->p[i] = mRot * this->p[i];
    }
    for (int i = 0; i < this->n.size(); i++) {
        this->n[i] = mRot * this->n[i];
    }
}

void Mesh::translate(VEC3 & offset)
{
    /* translate vertex positions */
    for (int i = 0; i < this->p.size(); i++) {
        this->p[i] = this->p[i] + offset;
    }
}

void Mesh::assemble(
    const Array<VEC3>& p, const Array<VEC3>& n, 
    const Array<VEC3>& t, const Array<INT3>& f, GLenum usage)
{
    unload();
    this->p = p;
    this->n = n;
    this->t = t;
    this->f = f;
    this->usage = usage;

    if (!isNull()) {
        /* model is loaded, now we need to upload data to GPU */
        /* vertex format is (X,Y,Z, U,V) */
        vbuf.gen(); /* initialize vertex fBuffer */
        packVertices();
        vbuf.loadData((void*)(packedVertices.data()), 
            sizeof(GL_VERTEX_XYZ_RGB_UV) * packedVertices.size(), 
            this->usage);
    }
}

bool Mesh::isNull() const
{
    if (this->f.size() < 3) /* cannot form a triangle */
        return true; 
    else return false;
}

triangle Mesh::getTriangle(const int & ind) const
{
    triangle t;
    t.p[0] = p[f[ind].x];
    t.p[1] = p[f[ind + 1].x];
    t.p[2] = p[f[ind + 2].x];
    return t;
}

triangleEx Mesh::getTriangleEx(const int & ind) const
{
    triangleEx t;
    t.p[0] = this->p[f[ind].x];
    t.p[1] = this->p[f[ind + 1].x];
    t.p[2] = this->p[f[ind + 2].x];
    t.t[0] = this->t[f[ind].y];
    t.t[1] = this->t[f[ind + 1].y];
    t.t[2] = this->t[f[ind + 2].y];
    t.n[0] = this->n[f[ind].z];
    t.n[1] = this->n[f[ind + 1].z];
    t.n[2] = this->n[f[ind + 2].z];
    return t;
}

aabb Mesh::getBoundingBox() const
{
    VEC3 bmin = VEC3(REAL(REAL_MAX), REAL(REAL_MAX), REAL(REAL_MAX));
    VEC3 bmax = VEC3(REAL(-REAL_MAX), REAL(-REAL_MAX), REAL(-REAL_MAX));

    for (int i = 0; i < f.size(); i+=3) {
        triangle t = getTriangle(i);
        for (int j = 0; j < 3; j++) {
            if (bmin.x > t.p[j].x) bmin.x = t.p[j].x;
            if (bmin.y > t.p[j].y) bmin.y = t.p[j].y;
            if (bmin.z > t.p[j].z) bmin.z = t.p[j].z;
            if (bmax.x < t.p[j].x) bmax.x = t.p[j].x;
            if (bmax.y < t.p[j].y) bmax.y = t.p[j].y;
            if (bmax.z < t.p[j].z) bmax.z = t.p[j].z;
        }
    }

    aabb box;
    box.bmin = bmin;
    box.bmax = bmax;

    return box;
}

Array<VEC3> Mesh::sampleMesh(const int& count) const
{
    Array<VEC3> sampled;
    if (isNull())
        return sampled;
    Array<REAL> pSums; /* partial sums */
    pSums.append(ZERO);
    REAL pSum = ZERO;
    for (int i = 0; i < f.size(); i+=3) {
        triangle t = getTriangle(i);
        VEC3& a = t.p[0], &b = t.p[1], &c = t.p[2];
        REAL tArea = REAL(0.5) * len(cross(a - c, b - c));
        pSum += tArea;
        pSums.append(pSum);
    }

    for (int i = 0; i < count; i++) {
        REAL r = uniform(ZERO, pSum);
        int mid = 0, left = 0, right = f.size() / 3-1;
        while (true) {
            mid = (left + right) / 2;
            REAL tx = pSums[mid + 1];
            REAL tn = pSums[mid];
            if (r > tx) { left = mid + 1; continue; }
            if (r < tn) { right = mid - 1; continue; }
            /* found! */
            left = right = mid;
            break;
        }

        REAL r1 = uniform(), r2 = uniform();
        triangle t = getTriangle(mid * 3);
        VEC3& a = t.p[0], &b = t.p[1], &c = t.p[2];
        /* https://www.cs.princeton.edu/~funk/tog02.pdf  */
        VEC3 samp = (1 - sqrt(r1)) * a + sqrt(r1)*(1 - r2)*b + sqrt(r1)*r2*c;
        sampled.append(samp);
    }

    return sampled;
}

bool Mesh::sampleMesh(const int& count, const char* savefile) const
{
    Array<VEC3> samples = sampleMesh(count);
    FILE* fp = NULL;
    if ((fp = fopen(savefile, "w")) == NULL) {
        return false;
    }
    char lineBuf[512];
    for (int i = 0; i < samples.size(); i++)
    {
        VEC3 p = samples[i];
        sprintf(lineBuf, "%.4f %.4f %.4f\n", float(p.x), float(p.y), float(p.z));
        fwrite(lineBuf, 1, strlen(lineBuf), fp);
    }
    fclose(fp);
    return true;
}

MorphMesh::MorphMesh()
{
}

MorphMesh::MorphMesh(const MorphMesh & that)
{
    this->meshArray = that.meshArray;
    this->meshInterp = that.meshInterp;
}

MorphMesh::~MorphMesh()
{
    unload();
}

MorphMesh & MorphMesh::operator=(const MorphMesh & that)
{
    unload();

    this->meshArray = that.meshArray;
    this->meshInterp = that.meshInterp;

    return (*this);
}

unsigned int MorphMesh::loadMorphTarget(const char * file)
{
    Mesh mesh;
    if (!mesh.loadOBJ(file, GL_STREAM_DRAW)) return 0;
    meshArray.append(mesh);
    unsigned int handle = meshArray.size(); /* note that handle = actual_index + 1 */

    /* create interpolated mesh when the first morph target is successfully loaded */
    if (handle == 1) { 
        /* just copy the first mesh data */
        /* in this process the GL context will also be prepared */
        meshInterp = meshArray[0];
    }
    return handle; 
}

bool MorphMesh::interpMeshes(int n, unsigned int * meshHandles, REAL * mixWeights)
{    
    /* interpolate vertex positions */
    int Nv = meshArray[meshHandles[0] - 1].getPositionArray().size(); /* note that handle = actual_index + 1 */
    for (int Vi = 0; Vi < Nv; Vi++) {
        VEC3 v = VEC3::zeros();
        for (int Mi = 0; Mi < n; Mi++) {
            VEC3 p = meshArray[meshHandles[Mi] - 1].getVertex(Vi);
            v = v + p * mixWeights[Mi];
        }
        meshInterp.setVertex(Vi, v);
    }
    /* upload interpolated mesh fBuffer data */
    meshInterp.upload();

    return false;
}

void MorphMesh::unload()
{
    meshArray.clear();
    meshInterp.unload();
}

void MorphMesh::draw()
{
    /* DO NOT FORGET TO SET MODEL, VIEW AND PROJECTION MATRICES */
    /* BEFORE CALLING THIS FUNCTION !!! */
    meshInterp.draw();
}

Mesh buildConvexMesh(Mesh & m)
{
    /* calculate convex hull based on point cloud */
    Array<VEC3> p = m.getPositionArray();

    convex cvex = buildConvex3d(p);
    if (cvex.vertices.size() == 0) {
        /* cannot build convex hull shape for this mesh,
           return empty mesh instead. */
        return Mesh();
    }

    Array<VEC3>& vertices = cvex.vertices;
    Array<INT3> faces; /* v/vt/vn */
    Array<VEC3> normals;
    Array<VEC3> texcoords;

    texcoords.append(VEC3(0, 0, 0));

    /* calculate surface normals for the convex shape */
    /* triangle winding is counter-clockwise */
    for (int i = 0; i < cvex.faces.size(); i++) {
        INT3 vert_idxs = cvex.faces[i];
        VEC3 v1 = cvex.vertices[vert_idxs.x];
        VEC3 v2 = cvex.vertices[vert_idxs.y];
        VEC3 v3 = cvex.vertices[vert_idxs.z];
        VEC3 v1v2 = v2 - v1;
        VEC3 v1v3 = v3 - v1;
        VEC3 n = normalize(cross(v1v2, v1v3));
        normals.append(n);

        /* build mesh */
        INT3 f;
        f.x = vert_idxs.x;        /* v */
        f.y = 0;                  /* vt */
        f.z = normals.size() - 1; /* vn: always the last (newly added) element */
        faces.append(f);
        f.x = vert_idxs.y;        /* v */
        f.y = 0;                  /* vt */
        f.z = normals.size() - 1; /* vn: always the last (newly added) element */
        faces.append(f);
        f.x = vert_idxs.z;        /* v */
        f.y = 0;                  /* vt */
        f.z = normals.size() - 1; /* vn: always the last (newly added) element */
        faces.append(f);
    }

    Mesh mesh;
    mesh.assemble(vertices, normals, texcoords, faces, GL_STATIC_DRAW);

    return mesh;
}

Mesh MeshMaker::box(const VEC3 & size)
{
    VEC3 hsize = REAL(0.5) * size;
    Array<VEC3> p;
    p.append(VEC3(-hsize.x, -hsize.y, -hsize.z)); /* v0 */
    p.append(VEC3(+hsize.x, -hsize.y, -hsize.z)); /* v1 */
    p.append(VEC3(+hsize.x, +hsize.y, -hsize.z)); /* v2 */
    p.append(VEC3(-hsize.x, +hsize.y, -hsize.z)); /* v3 */
    p.append(VEC3(-hsize.x, -hsize.y, +hsize.z)); /* v4 */
    p.append(VEC3(+hsize.x, -hsize.y, +hsize.z)); /* v5 */
    p.append(VEC3(+hsize.x, +hsize.y, +hsize.z)); /* v6 */
    p.append(VEC3(-hsize.x, +hsize.y, +hsize.z)); /* v7 */
    Array<VEC3> n;
    n.append(VEC3(REAL(0), REAL(-1), REAL(0))); /* n0 */
    n.append(VEC3(REAL(+1), REAL(0), REAL(0))); /* n1 */
    n.append(VEC3(REAL(0), REAL(+1), REAL(0))); /* n2 */
    n.append(VEC3(REAL(-1), REAL(0), REAL(0))); /* n3 */
    n.append(VEC3(REAL(0), REAL(0), REAL(-1))); /* n4 */
    n.append(VEC3(REAL(0), REAL(0), REAL(+1))); /* n5 */
    Array<VEC3> t;
    t.append(VEC3(REAL(0), REAL(0), REAL(0))); /* t0 */
    t.append(VEC3(REAL(1), REAL(0), REAL(0))); /* t1 */
    t.append(VEC3(REAL(1), REAL(1), REAL(0))); /* t2 */
    t.append(VEC3(REAL(0), REAL(1), REAL(0))); /* t3 */
    Array<INT3> f; /* p/t/n */
    f.append(INT3(0, 0, 0)); f.append(INT3(1, 1, 0)); f.append(INT3(5, 2, 0));
    f.append(INT3(0, 0, 0)); f.append(INT3(5, 2, 0)); f.append(INT3(4, 3, 0));
    f.append(INT3(1, 0, 1)); f.append(INT3(2, 1, 1)); f.append(INT3(6, 2, 1));
    f.append(INT3(1, 0, 1)); f.append(INT3(6, 2, 1)); f.append(INT3(5, 3, 1));
    f.append(INT3(2, 0, 2)); f.append(INT3(3, 1, 2)); f.append(INT3(7, 2, 2));
    f.append(INT3(2, 0, 2)); f.append(INT3(7, 2, 2)); f.append(INT3(6, 3, 2));
    f.append(INT3(3, 0, 3)); f.append(INT3(0, 1, 3)); f.append(INT3(4, 2, 3));
    f.append(INT3(3, 0, 3)); f.append(INT3(4, 2, 3)); f.append(INT3(7, 3, 3));
    f.append(INT3(0, 0, 4)); f.append(INT3(3, 1, 4)); f.append(INT3(2, 2, 4));
    f.append(INT3(0, 0, 4)); f.append(INT3(2, 2, 4)); f.append(INT3(1, 3, 4));
    f.append(INT3(5, 0, 5)); f.append(INT3(6, 1, 5)); f.append(INT3(7, 2, 5));
    f.append(INT3(5, 0, 5)); f.append(INT3(7, 2, 5)); f.append(INT3(4, 3, 5));
    Mesh m;
    m.assemble(p, n, t, f, GL_STATIC_DRAW);
    return m;
}

Mesh MeshMaker::stitchMeshes(const Array<Mesh*>& meshes)
{
   
    int length = int(Sqrt(REAL(meshes.size())) + REAL(0.5));
    if (length * length < meshes.size()) length++;
    int xMesh = 0, yMesh = 0;
    REAL unitTexOff = ONE / REAL(length);

    Array<VEC3> pS;
    Array<VEC3> nS;
    Array<VEC3> tS;
    Array<INT3> fS;
    int pOff = 0, nOff = 0, tOff = 0;
    for (int i = 0; i < meshes.size(); i++) {
        Mesh* mesh = meshes[i];
        Array<VEC3> p = mesh->getPositionArray();
        pS.append(p);
        nS.append(mesh->getNormalArray());
        Array<VEC3> t = mesh->getTexCoordArray();
        for (int ti = 0; ti < t.size(); ti++) {
            VEC3 ot = t[ti];
            /* wrap coord to (0,1) */
            ot.x = Fmod(ot.x, ONE); ot.y = Fmod(ot.y, ONE); ot.z = Fmod(ot.z, ONE);
            VEC3 dt = (ot + VEC3(REAL(xMesh), REAL(length - 1 - yMesh), ZERO)) * unitTexOff;
            tS.append(dt);
        }
        Array<INT3> f = mesh->getFaceIndexArray();
        /* v/t/n */
        for (int fi = 0; fi < f.size(); fi++)
            fS.append(f[fi] + INT3(pOff, tOff, nOff));

        xMesh++;
        if (xMesh >= length) { yMesh++; xMesh = 0; }
        pOff = pS.size(), nOff = nS.size(), tOff = tS.size();
    }

    Mesh mesh;
    mesh.assemble(pS, nS, tS, fS, GL_STATIC_DRAW);
    
    return mesh;
}

BSTMesh::BSTMeshNode::BSTMeshNode()
{
    left = right = NULL;
    depth = n = 0;
    idxs = NULL;
}

BSTMesh::BSTMeshNode::~BSTMeshNode()
{
    if (left) { delete left; left = NULL; }
    if (right) { delete right; right = NULL; }
    if (idxs) { free(idxs); idxs = NULL; }
}

bool BSTMesh::BSTMeshNode::isLeaf() const
{
    return (left == NULL || right == NULL);
}


BSTMesh::BSTMeshBuildSettings BSTMesh::defaultBuildSettings()
{
    BSTMeshBuildSettings settings;
    settings.maxDepth = 24;
    settings.sahResolution = 16;
    settings.splitThreshold = 4;
    return settings;
}

void BSTMesh::buildBST_splitBox(const aabb * box, const char axis, const REAL & value, aabb * left, aabb * right)
{
    switch (axis) {
    case 'x':
        left->bmin = box->bmin;
        left->bmax = VEC3(value, box->bmax.y, box->bmax.z);
        right->bmin = VEC3(value, box->bmin.y, box->bmin.z);
        right->bmax = box->bmax;
        break;
    case 'y':
        left->bmin = box->bmin;
        left->bmax = VEC3(box->bmax.x, value, box->bmax.z);
        right->bmin = VEC3(box->bmin.x, value, box->bmin.z);
        right->bmax = box->bmax;
        break;
    case 'z':
        left->bmin = box->bmin;
        left->bmax = VEC3(box->bmax.x, box->bmax.y, value);
        right->bmin = VEC3(box->bmin.x, box->bmin.y, value);
        right->bmax = box->bmax;
        break;
    default:
        break;
    }
}

void BSTMesh::buildBST_sahCost(
    const Array<int>& idxs, 
    const aabb& left, const aabb& right,
    const char& axis, const REAL& value, 
    REAL * cost)
{
    int leftCount = 0, rightCount = 0;
    for (int i = 0; i < idxs.size(); i++) {
        triangle t = getTriangle(idxs[i]);
        if (isIntersect(t, left)) leftCount++;
        if (isIntersect(t, right)) rightCount++;
    }
    *cost = surfaceArea(left) * REAL(leftCount) + surfaceArea(right) * REAL(rightCount);
}

void BSTMesh::buildBST_recur(BSTMeshNode * node, const Array<int>& idxs)
{
    bool shouldSplit = true;
    aabb lBox, rBox;
    Array<int> lIdxs, rIdxs;
    BSTMeshNodeType type = BSTMeshNodeType::Invalid;

    /* determine if this node should split */
    if (node->depth < this->bstSettings.maxDepth && 
        idxs.size() > this->bstSettings.splitThreshold)
    {
        /* find best split */
        VEC3 dxyz = node->box.bmax - node->box.bmin;

        /* init axes and splits (three directions) */
        int& sahResolution = this->bstSettings.sahResolution;
        REAL* costs = new REAL[sahResolution * 3];
        REAL* splits = new REAL[sahResolution * 3];
        char* axes = new char[sahResolution * 3];
        for (int i = 0; i < sahResolution * 3; i++) {
            int D = i / sahResolution, S = i % sahResolution;
            REAL step = dxyz.e[D] / REAL(sahResolution + 1);
            axes[i] = 'x' + D;
            splits[i] = node->box.bmin.e[D] + REAL(S + 1) * step;
        }

        /* find minimum cost for all dimensions */
        for (int i = 0; i < sahResolution * 3; i++) {
            buildBST_splitBox(&(node->box), axes[i], splits[i], &lBox, &rBox);
            buildBST_sahCost(idxs, lBox, rBox, axes[i], splits[i], &(costs[i]));
        }
        int ax = argmin(costs, 3 * sahResolution);
        char axis = axes[ax];
        REAL split = splits[ax];

        buildBST_splitBox(&(node->box), axis, split, &lBox, &rBox);

        for (int i = 0; i < idxs.size(); i++) {
            triangle t = getTriangle(idxs[i]);

            if (isIntersect(t, lBox)) 
                lIdxs.append(idxs[i]);
            if (isIntersect(t, rBox)) 
                rIdxs.append(idxs[i]);

        }
        /* calculate volume overlap */
        const REAL maxIncreaseRatio = REAL(0.40);
        REAL increaseRatio = ONE - REAL(lIdxs.size() + rIdxs.size()) / REAL(idxs.size());
        if (node->depth > this->bstSettings.maxDepth * REAL(0.75) &&
            increaseRatio > maxIncreaseRatio) {
            shouldSplit = false;
            type = BSTMeshNodeType::AtomLeaf;
        }
        else
        {
            shouldSplit = true;
            type = BSTMeshNodeType::NonLeaf;
        }

        delete[] costs;
        delete[] splits;
        delete[] axes;
    }
    else {
        shouldSplit = false;
        if (idxs.size() <= this->bstSettings.splitThreshold)
            type = BSTMeshNodeType::OrdinaryLeaf;
        else if (node->depth == this->bstSettings.maxDepth)
            type = BSTMeshNodeType::DeepLeaf;
    }

    node->type = type;

    /* * * split topNode into two child nodes or turn it to leaf topNode * * */

    if (shouldSplit){
        /* build nodes */
        node->left = new BSTMeshNode();
        node->right = new BSTMeshNode();
        node->left->depth = node->depth + 1;
        node->right->depth = node->depth + 1;
        node->left->box = lBox;
        node->right->box = rBox;

        /* resursive build */
        buildBST_recur(node->left, lIdxs);
        buildBST_recur(node->right, rIdxs);
    }
    else {
        node->n = idxs.size();
        node->idxs = (int*)malloc(sizeof(int)*(node->n));
        for (int i = 0; i < node->n; i++) {
            node->idxs[i] = idxs[i];
        }
    }
}

void BSTMesh::buildBST_boxFit(const Array<int>& idxs, VEC3* bmin, VEC3* bmax)
{
    VEC3 _bmin = VEC3(REAL(REAL_MAX), REAL(REAL_MAX), REAL(REAL_MAX));
    VEC3 _bmax = VEC3(REAL(-REAL_MAX), REAL(-REAL_MAX), REAL(-REAL_MAX));

    for (int i = 0; i < idxs.size(); i++) {
        triangle t = getTriangle(idxs[i]);
        for (int j = 0; j < 3; j++) {
            if (_bmin.x > t.p[j].x) _bmin.x = t.p[j].x;
            if (_bmin.y > t.p[j].y) _bmin.y = t.p[j].y;
            if (_bmin.z > t.p[j].z) _bmin.z = t.p[j].z;
            if (_bmax.x < t.p[j].x) _bmax.x = t.p[j].x;
            if (_bmax.y < t.p[j].y) _bmax.y = t.p[j].y;
            if (_bmax.z < t.p[j].z) _bmax.z = t.p[j].z;
        }
    }

    *bmin = _bmin;
    *bmax = _bmax;
}

BSTMesh::BSTMesh()
{
    bstRoot = NULL;
}

BSTMesh::~BSTMesh()
{
    unload();
}

void BSTMesh::buildBST(const BSTMeshBuildSettings& settings)
{
    clearBST();

    /* collect all faces */
    Array<int> allFaces;
    for (int i = 0; i < f.size(); i += 3)
        allFaces.append(i);

    bstRoot = new BSTMeshNode();
    bstRoot->depth = 1;
    buildBST_boxFit(allFaces, &(bstRoot->box.bmin), &(bstRoot->box.bmax));
    buildBST_recur(bstRoot, allFaces);
}

void BSTMesh::clearBST_recur(BSTMeshNode * node)
{
    if (node == NULL) return;
    if (node->left) {
        clearBST_recur(node->left);
        node->left = NULL;
    }
    if (node->right) {
        clearBST_recur(node->right);
        node->right = NULL;
    }
    if (node->idxs) {
        free(node->idxs);
        node->idxs = NULL;
    }
    delete node;
}

void BSTMesh::clearBST()
{
    clearBST_recur(bstRoot);
    bstRoot = NULL;
}

void BSTMesh::bstAdvInit()
{
    bstShader.gen(
        /* vertex shader */
        "#version 330 core                                                   \n"
        "layout(location = 0) in vec3 iPosition;                             \n"
        "uniform mat4 model;                                                 \n"
        "uniform mat4 view;                                                  \n"
        "uniform mat4 projection;                                            \n"
        "void main()                                                         \n"
        "{                                                                   \n"
        "	gl_Position = projection * view * model * vec4(iPosition, 1.0);  \n"
        "}                                                                   \n",
        /* fragment shader */
        "#version 330 core                                                   \n"
        "out vec4 FragColor;                                                 \n"
        "uniform vec3 Color;                                                 \n"
        "void main()                                                         \n"
        "{                                                                   \n"
        "	FragColor = vec4(Color, 1.0);                                    \n"
        "}                                                                   \n"
    );
    bstVB.gen();
    Array<VEC3> p = MeshMaker::box(VEC3(ONE, ONE, ONE)).getPositionArray();
    Array<GL_VERTEX_XYZ> V;
    for (int i = 0; i < 8; i++) {
        GL_VERTEX_XYZ v;
        v.position[0] = float(p[i].x);
        v.position[1] = float(p[i].y);
        v.position[2] = float(p[i].z);
        V.append(v);
    }
    int indices[] = { 0,1,1,2,2,3,3,0,4,5,5,6,6,7,7,4,0,4,1,5,2,6,3,7 };
    bstVB.loadData(V.data(), sizeof(float) * 3 * V.size(), indices, sizeof(int) * 24, GL_STATIC_DRAW);

}

void BSTMesh::bstAdvFree()
{
    bstShader.unload();
    bstVB.unload();
}

bool BSTMesh::loadOBJ(const char * file, const BSTMeshBuildSettings& settings, GLenum usage)
{
    this->bstSettings = settings;
    bool success = Mesh::loadOBJ(file, usage);
    if (success)
        buildBST(settings);
    bstAdvInit();
    return success;
}

void BSTMesh::assemble(const Array<VEC3>& p, const Array<VEC3>& n, 
    const Array<VEC3>& t, const Array<INT3>& f, 
    const BSTMeshBuildSettings& settings, GLenum usage)
{
    this->bstSettings = settings;
    Mesh::assemble(p, n, t, f, usage);
    buildBST(settings);
    bstAdvInit();
}

void BSTMesh::unload()
{
    Mesh::unload();
    clearBST();
    bstAdvFree();
}

void BSTMesh::rotate(QUAT & q)
{
    Mesh::rotate(q);
    /* rebuild BST */
    if (isNull() == false)
        buildBST(this->bstSettings);
}

void BSTMesh::translate(VEC3 & offset)
{
    Mesh::translate(offset);
    /* rebuild BST */
    if (isNull() == false)
        buildBST(this->bstSettings);
}

void BSTMesh::saveNode_recur(BSTMeshNode * node, FILE * fp)
{
    /* box */
    int iBuffer;
    float fBuffer;
    fBuffer = float(node->box.bmin.x); fwrite(&fBuffer, 4, 1, fp);
    fBuffer = float(node->box.bmin.y); fwrite(&fBuffer, 4, 1, fp);
    fBuffer = float(node->box.bmin.z); fwrite(&fBuffer, 4, 1, fp);
    fBuffer = float(node->box.bmax.x); fwrite(&fBuffer, 4, 1, fp);
    fBuffer = float(node->box.bmax.y); fwrite(&fBuffer, 4, 1, fp);
    fBuffer = float(node->box.bmax.z); fwrite(&fBuffer, 4, 1, fp);
    /* depth */
    iBuffer = node->depth; fwrite(&iBuffer, 4, 1, fp);
    /* n */
    iBuffer = node->n; fwrite(&iBuffer, 4, 1, fp);
    /* type */
    iBuffer = int(node->type); fwrite(&iBuffer, 4, 1, fp);
    /* idxs */
    if (node->isLeaf()) {
        fwrite("#", 1, 1, fp);        
        for (int i = 0; i < node->n; i++) {
            iBuffer = node->idxs[i]; fwrite(&iBuffer, 4, 1, fp);
        }
    }
    else {
        fwrite("@", 1, 1, fp);
        saveNode_recur(node->left, fp);
        saveNode_recur(node->right, fp);
    }
}

void BSTMesh::loadNode_recur(BSTMeshNode ** node, FILE * fp)
{
    int iBuffer;
    float fBuffer;
    char cBuffer;

    *node = new BSTMeshNode();
    fread(&fBuffer, 4, 1, fp); (*node)->box.bmin.x = REAL(fBuffer);
    fread(&fBuffer, 4, 1, fp); (*node)->box.bmin.y = REAL(fBuffer);
    fread(&fBuffer, 4, 1, fp); (*node)->box.bmin.z = REAL(fBuffer);
    fread(&fBuffer, 4, 1, fp); (*node)->box.bmax.x = REAL(fBuffer);
    fread(&fBuffer, 4, 1, fp); (*node)->box.bmax.y = REAL(fBuffer);
    fread(&fBuffer, 4, 1, fp); (*node)->box.bmax.z = REAL(fBuffer);

    fread(&iBuffer, 4, 1, fp); (*node)->depth = iBuffer;
    fread(&iBuffer, 4, 1, fp); (*node)->n = iBuffer;
    fread(&iBuffer, 4, 1, fp); (*node)->type = BSTMeshNodeType(iBuffer);
    fread(&cBuffer, 1, 1, fp);
    bool isLeaf = cBuffer == '#' ? true : false;
    if (isLeaf) {
        if ((*node)->n > 0)
            (*node)->idxs = (int*)malloc(sizeof(int)*(*node)->n);

        for (int i = 0; i < (*node)->n; i++) {
            fread(&iBuffer, 4, 1, fp); (*node)->idxs[i] = iBuffer; 
        }
    }
    else {
        loadNode_recur(&((*node)->left), fp);
        loadNode_recur(&((*node)->right), fp);
    }

}

bool BSTMesh::saveAsset(const char * file)
{
    if (isNull() || bstRoot == NULL)
        return false;
    FILE* fp = NULL;
    if ((fp = fopen(file, "wb")) == NULL)
        return false;

    const char magicString[64] = "rkJYme3skGz11EQtpwTuCh";
    const char className[64] = "BSTMesh";
    fwrite(magicString, 1, 64, fp);
    fwrite(className, 1, 64, fp);
    int iBuffer;
    float fBuffer;

    /* write BST build settings */
    fwrite(&this->bstSettings, sizeof(BSTMeshBuildSettings), 1, fp);

    /* write all vertices (save as single precision float) */
    iBuffer = p.size(); fwrite(&iBuffer, 4, 1, fp);
    for (int i = 0; i < p.size(); i++) {
        VEC3 data = p[i];
        fBuffer = float(data.x); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.y); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.z); fwrite(&fBuffer, 4, 1, fp);
    }
    /* save normals */
    iBuffer = n.size(); fwrite(&iBuffer, 4, 1, fp);
    for (int i = 0; i < n.size(); i++) {
        VEC3 data = n[i];
        fBuffer = float(data.x); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.y); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.z); fwrite(&fBuffer, 4, 1, fp);
    }
    /* save texture coordinates */
    iBuffer = t.size(); fwrite(&iBuffer, 4, 1, fp);
    for (int i = 0; i < t.size(); i++) {
        VEC3 data = t[i];
        fBuffer = float(data.x); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.y); fwrite(&fBuffer, 4, 1, fp);
        fBuffer = float(data.z); fwrite(&fBuffer, 4, 1, fp);
    }
    /* save faces */
    iBuffer = f.size(); fwrite(&iBuffer, 4, 1, fp);
    for (int i = 0; i < f.size(); i++) {
        INT3 data = f[i]; /* v/t/n */
        iBuffer = data.x; fwrite(&iBuffer, 4, 1, fp);
        iBuffer = data.y; fwrite(&iBuffer, 4, 1, fp);
        iBuffer = data.z; fwrite(&iBuffer, 4, 1, fp);
    }
    /* save BST */
    saveNode_recur(bstRoot, fp);

    fclose(fp);
    return true;
}

bool BSTMesh::loadAsset(const char * file)
{
    unload();

    FILE* fp = NULL;
    if ((fp = fopen(file, "rb")) == NULL)
        return false;

    char cBuffer[64] = { 0 };
    const char magicString[64] = "rkJYme3skGz11EQtpwTuCh";
    const char className[64] = "BSTMesh";
    fread(cBuffer, 1, 64, fp);
    if (strcmp(cBuffer, magicString) != 0) {
        fclose(fp);
        return false;
    }
    fread(cBuffer, 1, 64, fp);
    if (strcmp(cBuffer, className) != 0) {
        fclose(fp);
        return false;
    }
    int iBuffer;
    float fBuffer;


    /* read BST build settings */
    fread(&this->bstSettings, sizeof(BSTMeshBuildSettings), 1, fp);


    /* p */
    fread(&iBuffer, 4, 1, fp);
    for (int i = 0; i < iBuffer; i++) {
        VEC3 v;
        fread(&fBuffer, 4, 1, fp); v.x = REAL(fBuffer);
        fread(&fBuffer, 4, 1, fp); v.y = REAL(fBuffer);
        fread(&fBuffer, 4, 1, fp); v.z = REAL(fBuffer);
        p.append(v);
    }
    /* n */
    fread(&iBuffer, 4, 1, fp);
    for (int i = 0; i < iBuffer; i++) {
        VEC3 v;
        fread(&fBuffer, 4, 1, fp); v.x = REAL(fBuffer);
        fread(&fBuffer, 4, 1, fp); v.y = REAL(fBuffer);
        fread(&fBuffer, 4, 1, fp); v.z = REAL(fBuffer);
        n.append(v);
    }
    /* t */
    fread(&iBuffer, 4, 1, fp);
    for (int i = 0; i < iBuffer; i++) {
        VEC3 v;
        fread(&fBuffer, 4, 1, fp); v.x = REAL(fBuffer);
        fread(&fBuffer, 4, 1, fp); v.y = REAL(fBuffer);
        fread(&fBuffer, 4, 1, fp); v.z = REAL(fBuffer);
        t.append(v);
    }
    /* f */
    fread(&iBuffer, 4, 1, fp);
    for (int i = 0; i < iBuffer; i++) {
        INT3 v;
        fread(&v.x, 4, 1, fp);
        fread(&v.y, 4, 1, fp);
        fread(&v.z, 4, 1, fp);
        f.append(v);
    }

    loadNode_recur(&bstRoot, fp);

    upload(); /* upload data to GPU */

    bstAdvInit();

    return true;
}

void BSTMesh::drawNodeBoundary(GLCamera* camera, const VEC3& color, int minDepth, int maxDepth)
{
    if (bstRoot == NULL) return;

    GLint lastPolyModes[2]; /* front and back */

    glGetIntegerv(GL_POLYGON_MODE, lastPolyModes);
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    bstShader.use();
    bstShader.setUniformMatrix4f("view", camera->getViewMatrix());
    bstShader.setUniformMatrix4f("projection", camera->getProjectionMatrix());

    drawNodeBoundary_recur(camera, bstRoot, minDepth, maxDepth);

    /* restore polygon mode */
    glPolygonMode(GL_FRONT_AND_BACK, lastPolyModes[0]);
}

void BSTMesh::drawNodeBoundary_recur(GLCamera* camera,
    BSTMeshNode * node, int minDepth, int maxDepth)
{
    VEC3 p = (node->box.bmin + node->box.bmax) / REAL(2.0);
    VEC3 l = node->box.bmax - node->box.bmin;

    glm::mat4 model = glm::translate(glm::vec3(p.x, p.y, p.z)) * glm::scale(glm::vec3(l.x, l.y, l.z));

    if (node->depth >= minDepth && node->depth <= maxDepth) {
        VEC3 color;
        switch (node->type)
        {
        case BSTMeshNodeType::NonLeaf:
            color = VEC3(ONE, ONE, ZERO); /* yellow */
            break;
        case BSTMeshNodeType::OrdinaryLeaf:
            color = VEC3(ZERO, ONE, ZERO); /* green */
            break;
        case BSTMeshNodeType::AtomLeaf:
            color = VEC3(ONE, ZERO, ZERO); /* red */
            break;
        case BSTMeshNodeType::DeepLeaf:
            color = VEC3(ZERO, ZERO, ONE); /* blue */
            break;
        default:
            color = VEC3(ONE, ONE, ONE); /* white */
        }
        bstShader.setUniform3f("Color", color);
        bstShader.setUniformMatrix4f("model", model);
        bstVB.draw(GL_LINES, 24, GL_UNSIGNED_INT, NULL);
    }

    if (node->left) drawNodeBoundary_recur(camera, node->left, minDepth, maxDepth);
    if (node->right) drawNodeBoundary_recur(camera, node->right, minDepth, maxDepth);
}


void BSTMesh::doStatistics(char** log)
{
    /* collect all faces */
    Array<int> allFaces;
    for (int i = 0; i < f.size(); i += 3)
        allFaces.append(i);

    if (bstRoot == NULL)
        return;

    bstStat.clear();
    doStatistics_recur(bstRoot);
    bstStat.overlapRatio = REAL(bstStat.nLeafTrisSum) / REAL(allFaces.size());

    if (*log != NULL) {
        /* User passed a non NULL pointer, which is unsafe. It */
        /* can be a wild pointer or a pointer pointing to a    */
        /* valid memory block. So simply overwrite its value   */
        /* is not permitted. */
        return;
    }
    else {
        *log = (char*)malloc(BSTMESH_DBG_LOGBUFFER_SIZE);
        memset(*log, 0, BSTMESH_DBG_LOGBUFFER_SIZE);
    }

    char linebuf[BSTMESH_DBG_LOGBUFFER_SIZE];

    int nTotalLeaf = bstStat.nAtomLeaf + bstStat.nDeepLeaf + bstStat.nOrdinaryLeaf;
    strcat(*log, "* * * * * * * BST Mesh Statistics * * * * * * *\n");
    sprintf(linebuf, ">>> Size:  %d bytes\n", bstStat.nBytes);
    strcat(*log, linebuf);
    sprintf(linebuf, ">>> Nodes:     total = %d\n", bstStat.nTotalNodes);
    strcat(*log, linebuf);
    sprintf(linebuf, "            non-leaf = %d\n", bstStat.nNonLeaf);
    strcat(*log, linebuf);
    sprintf(linebuf, "                leaf = %d\n", nTotalLeaf);
    strcat(*log, linebuf);
    if (bstStat.nInvalidNodes > 0) {
        sprintf(linebuf, "             invalid = %d\n", bstStat.nInvalidNodes);
        strcat(*log, linebuf);
    }
    sprintf(linebuf, ">>> Leaves: ordinary = %d (%3.2f%%)\n", bstStat.nOrdinaryLeaf,
        REAL(bstStat.nOrdinaryLeaf) / REAL(nTotalLeaf) * REAL(100.0));
    strcat(*log, linebuf);
    sprintf(linebuf, "                deep = %d (%3.2f%%)\n", bstStat.nDeepLeaf,
        REAL(bstStat.nDeepLeaf) / REAL(nTotalLeaf) * REAL(100.0));
    strcat(*log, linebuf);
    sprintf(linebuf, "                atom = %d (%3.2f%%)\n", bstStat.nAtomLeaf,
        REAL(bstStat.nAtomLeaf) / REAL(nTotalLeaf) * REAL(100.0));
    strcat(*log, linebuf);
    sprintf(linebuf, "            max_tris = %d\n", bstStat.nMaxLeafTris);
    strcat(*log, linebuf);

    strcat(*log, "* * * * * * * * * * * * * * * * * * * * * * * *\n");

}

void BSTMesh::doRaytracing(const ray & r, ray_hit& rtResult, const REAL& searchDist)
{
    /* here i don't use recursion, instead i maintain a node pointer stack */
    Stack<BSTMeshNode*> nodeStack;
    nodeStack.reserve(bstSettings.maxDepth + 1);
    nodeStack.push(bstRoot); /* initialize stack and trace */

    REAL t0 = REAL_MAX, u0, v0;
    int I = -1;

    while (nodeStack.isEmpty() == false) {
        BSTMeshNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                triangle A = getTriangle(topNode->idxs[i]);
                REAL u, v, t;
                if (isIntersect(r, A, u, v, t)) {
                    if (t0 > t && t < searchDist) {
                        t0 = t;
                        u0 = u;
                        v0 = v;
                        I = topNode->idxs[i];
                    }
                }
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(r, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(r, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }

    if (I != -1) {
        triangleEx t = getTriangleEx(I);
        rtResult.hit = true;
        rtResult.point = t.p[0] * u0 + t.p[1] * v0 + t.p[2] * (1 - u0 - v0);
        rtResult.normal = normalize(t.n[0] * u0 + t.n[1] * v0 + t.n[2] * (1 - u0 - v0));
        rtResult.uvw = t.t[0] * u0 + t.t[1] * v0 + t.t[2] * (1 - u0 - v0);
        rtResult.dist = t0;
    }
    else {
        rtResult.hit = false;
    }
}

void BSTMesh::doRaytracing(const ray & r, Array<ray_hit>& rtResults, const REAL& searchDist)
{
    /* here i don't use recursion, instead i maintain a node pointer stack */
    Stack<BSTMeshNode*> nodeStack;
    nodeStack.reserve(bstSettings.maxDepth + 1);
    nodeStack.push(bstRoot); /* initialize stack and trace */

    while (nodeStack.isEmpty() == false) {
        BSTMeshNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                int triangleId = topNode->idxs[i];
                triangleEx A = getTriangleEx(triangleId);
                REAL u, v, dist;
                if (isIntersect(r, A, u, v, dist)) {
                    if (dist < searchDist) {
                        /* a triangle may be found in multiple leaf nodes */
                        /* so a ray may hit a triangle "multiple" times, so */
                        /* here i will do an extra check to ensure each */
                        /* triangle will only have one record */
                        bool found_duplicate = false;
                        for (int i = 0; i < rtResults.size(); i++) {
                            if (rtResults[i].triangleId == triangleId) {
                                found_duplicate = true;
                                break;
                            }
                        }
                        if (!found_duplicate) {
                            ray_hit rtResult;
                            rtResult.hit = true;
                            rtResult.dist = dist;
                            rtResult.point = A.p[0] * u + A.p[1] * v + A.p[2] * (1 - u - v);
                            rtResult.normal = normalize(A.n[0] * u + A.n[1] * v + A.n[2] * (1 - u - v));
                            rtResult.uvw = A.t[0] * u + A.t[1] * v + A.t[2] * (1 - u - v);
                            rtResult.triangleId = triangleId;
                            rtResults.append(rtResult);
                        }
                    }
                }
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(r, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(r, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
}

void BSTMesh::doRaytracing(const VEC3 & pos, const VEC3 & look, const VEC3 & up, 
    const REAL & fov, const int & w, const int & h, const char * image)
{
    /* calculate raytrace parameters */
    VEC3 F = normalize(look - pos);   /* front */
    VEC3 R = normalize(cross(F, up)); /* right */
    VEC3 U = normalize(cross(R, F));  /* up */
    VEC3 X = 2 * Tan(fov / TWO) * R;
    VEC3 Y = len(X) * REAL(h) / REAL(w) * (-U);
    VEC3 O = pos + F - X / TWO - Y / TWO;
    VEC3 dX = X / REAL(w);
    VEC3 dY = Y / REAL(h);

    /* allocate render fBuffer */
    BSTMeshRenderBuffer im;
    im.pixels = (VEC3*)malloc(sizeof(VEC3) * w * h);
    im.w = w;
    im.h = h;

    /* start rendering for each pixel (single thread) */
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            ray_hit rt;
            VEC3 c; /* color */
            VEC3 T = O + dX * (REAL(x) + REAL(0.5)) + dY * (REAL(y) + REAL(0.5));
            ray r(pos, T - pos);
            doRaytracing(r, rt, REAL_MAX);
            if (rt.hit) {
                c.x = c.y = c.z = (ONE + dot(VEC3(0, 1, 0), rt.normal)) / TWO;
            }
            im.pixels[y*w + x] = c;
        }
    }

    /* save image to disk (default format PNG) */
    BYTE* rawbuf = (BYTE*)malloc(sizeof(BYTE) * 3 * w * h);
    for (int i = 0; i < w * h; i++) {
        rawbuf[i * 3] = doRT_byteClamp(im.pixels[i].x * REAL(255.0));
        rawbuf[i * 3 + 1] = doRT_byteClamp(im.pixels[i].y * REAL(255.0));
        rawbuf[i * 3 + 2] = doRT_byteClamp(im.pixels[i].z * REAL(255.0));
    }
    stbi_write_png(image, w, h, 3, rawbuf, w * 3);

    free(rawbuf);
    free(im.pixels);

}

void BSTMesh::doSphereCrossing(const sphere & s, Array<int>& idxs)
{
    /* here i don't use recursion, instead i maintain a node pointer stack */
    Stack<BSTMeshNode*> nodeStack;
    nodeStack.reserve(bstSettings.maxDepth + 1);
    nodeStack.push(bstRoot); /* initialize stack and trace */

    while (nodeStack.isEmpty() == false) {
        BSTMeshNode* topNode = NULL;
        nodeStack.pop(topNode);
        if (topNode->n > 0) {
            /* leaf node */
            for (int i = 0; i < topNode->n; i++) {
                triangle A = getTriangle(topNode->idxs[i]);
                if (isIntersect(s, A))
                    idxs.append(topNode->idxs[i]);
            }
        }
        else {
            /* non leaf node */
            if (topNode->left && isIntersect(s, topNode->left->box))
                nodeStack.push(topNode->left);
            if (topNode->right && isIntersect(s, topNode->right->box))
                nodeStack.push(topNode->right);
        }
    }
}

void BSTMesh::doStatistics_recur(BSTMeshNode * node)
{
    bstStat.nTotalNodes++;
    bstStat.nBytes += sizeof(BSTMeshNode) + node->n * sizeof(int);
    if (node->n > bstStat.nMaxLeafTris)
        bstStat.nMaxLeafTris = node->n;
    switch (node->type) {
    case BSTMeshNodeType::OrdinaryLeaf:
        bstStat.nOrdinaryLeaf++;
        bstStat.nLeafTrisSum += node->n;
        bstStat.leafTris.append(node->n);
        break;
    case BSTMeshNodeType::AtomLeaf:
        bstStat.nAtomLeaf++;
        bstStat.nLeafTrisSum += node->n;
        bstStat.leafTris.append(node->n);
        break;
    case BSTMeshNodeType::DeepLeaf:
        bstStat.nDeepLeaf++;
        bstStat.nLeafTrisSum += node->n;
        bstStat.leafTris.append(node->n);
        break;
    case BSTMeshNodeType::Invalid:
        bstStat.nInvalidNodes++;
        break;
    case BSTMeshNodeType::NonLeaf:
        bstStat.nNonLeaf++;
        break;
    }
    if (node->left) doStatistics_recur(node->left);
    if (node->right) doStatistics_recur(node->right);
}

BYTE BSTMesh::doRT_byteClamp(const REAL & value)
{
    if (value < 0) return BYTE(0);
    else if (value > 255) return BYTE(255);
    else return BYTE(value);
}

BSTMesh::BSTMeshStatistics::BSTMeshStatistics()
{
    clear();
}

void BSTMesh::BSTMeshStatistics::clear()
{
    nTotalNodes = 0;
    nInvalidNodes = 0;
    nNonLeaf = 0;
    nOrdinaryLeaf = 0;
    nAtomLeaf = 0;
    nDeepLeaf = 0;
    nLeafTrisSum = 0;
    overlapRatio = REAL(0.0);
    nMaxLeafTris = 0;
    nBytes = 0;
}

BSTMesh::BSTMeshBuildSettings::BSTMeshBuildSettings()
{
    maxDepth = 0;
    splitThreshold = 0;
    sahResolution = 0;
}

bool BSTMesh::BSTMeshBuildSettings::isValid()
{
    if (maxDepth < 1) return false;
    if (splitThreshold < 1) return false;
    if (sahResolution < 1) return false;
    return true;
}
