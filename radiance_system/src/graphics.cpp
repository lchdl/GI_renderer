#include "graphics.h"

/* some utility functions and base classes */

int _getCurrentFrameBuffer() {
    int fboID;
    glDebugCheck(glGetIntegerv(GL_FRAMEBUFFER_BINDING, &fboID));
    return fboID;
}

_GLWindowMargin::_GLWindowMargin() { windowWidth = 0; windowHeight = 0; }

void _GLWindowMargin::setWindowMargin(int w, int h) { windowWidth = w; windowHeight = h; }

glm::vec2 _GLWindowMargin::getWindowMargin() { return glm::vec2(windowWidth, windowHeight); }

bool _GLWindowMargin::isWindowMarginSet() {
    if (windowWidth <= 0 || windowHeight <= 0) return false;
    else return true;
}

GLVertexBuffer::GLVertexBuffer() { vboID = 0; vaoID = 0; eboID = 0; }

void GLVertexBuffer::gen() {
    /* in this function, you need to implement */
    /* the behavior when the buffer is being   */
    /* generated */
    _createHandles();
}

void GLVertexBuffer::loadData(const void * data, int size, GLenum usage) {
    /* defines the behavior when loading data into */
    /* vertex buffer (which is located in GPU) */
    /* NOTE: glBufferData() completely deletes the original buffer */
    /* data uploaded to the GPU and creates a new one */
    if (!isInit()) gen(); /* generate buffer if it is not created */
    glDebugCheck(glBindVertexArray(vaoID));
    glDebugCheck(glBindBuffer(GL_ARRAY_BUFFER, vboID));
    glDebugCheck(glBufferData(GL_ARRAY_BUFFER, size, data, usage));
    defineAttribute();
    /* unbind VBO and VAO for safety */
    glDebugCheck(glBindBuffer(GL_ARRAY_BUFFER, 0));
    glDebugCheck(glBindVertexArray(0));
}

void GLVertexBuffer::loadData(const void * data, int dataSize, const void * elem, int elemSize, GLenum usage) {
    /* defines the behavior when loading data into */
    /* vertex buffer (which is located in GPU) */
    /* NOTE: glBufferData() completely deletes the original buffer */
    /* data uploaded to the GPU and creates a new one */

    /* elem: pointer to index buffer (usually an array of integers) */
    /* elemSize: size of index buffer measured in bytes */

    if (!isInit()) gen(); /* generate buffer if it is not created */
    glBindVertexArray(vaoID);
    glBindBuffer(GL_ARRAY_BUFFER, vboID);
    glBufferData(GL_ARRAY_BUFFER, dataSize, data, usage);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, elemSize, elem, usage);
    defineAttribute();
    /* unbind VBO and VAO for safety */
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void GLVertexBuffer::updateData(const void * data, int offset, int size)
{
    glBindVertexArray(vaoID);
    glBindBuffer(GL_ARRAY_BUFFER, vboID);
    glBufferSubData(GL_ARRAY_BUFFER, offset, size, data);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void GLVertexBuffer::draw(GLenum mode, GLint first, GLsizei count) {
    /* draw mesh using glDrawArrays() */
    /* set this buffer object to be the current active buffer */
    /* using glBindBuffer() and draw onto screen using glDrawArrays() */
    glBindVertexArray(vaoID);
    glDrawArrays(mode, first, count);
}

void GLVertexBuffer::draw(GLenum mode, GLsizei count, GLenum type, const void * indices) {
    /* draw mesh using glDrawElements() */
    glBindVertexArray(vaoID);
    glDrawElements(mode, count, type, indices);
}

void GLVertexBuffer::unload() { _deleteHandles(); }

bool GLVertexBuffer::isInit() {
    if (vboID == 0 || vaoID == 0 || eboID == 0) return false;
    else return true;
}

void GLVertexBuffer::_createHandles() {
    if (vboID == 0) glGenBuffers(1, &vboID);
    if (vaoID == 0) glGenVertexArrays(1, &vaoID);
    if (eboID == 0) glGenBuffers(1, &eboID);
}

void GLVertexBuffer::_deleteHandles() {
    if (vboID != 0) { 
        glDeleteBuffers(1, &vboID); 
        vboID = 0; 
    }
    if (vaoID != 0) { 
        glDeleteVertexArrays(1, &vaoID); 
        vaoID = 0; 
    }
    if (eboID != 0) { 
        glDeleteBuffers(1, &eboID); 
        eboID = 0; 
    }
}

void GLVertexBuffer_XYZ::defineAttribute() {
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

}

void GLVertexBuffer_XYZ_RGB_UV::defineAttribute() {
    // position attribute (x,y,z)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute (r,g,b)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coord attribute (u,v)
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
}

void GLVertexBuffer_XYZ_UV::defineAttribute() {
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // texture coord attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

void GLVertexBuffer_XYZ_RGB::defineAttribute()
{
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // texture coord attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

GLTexture2D::GLTexture2D() {
    texID = 0; width = 0; height = 0; channels = 0;
    internalFormat = GL_INVALID_VALUE;
}

GLTexture2D::~GLTexture2D() { unload(); }

bool GLTexture2D::gen(int w, int h, GLint internalFormat, GLint sampling) {
    if (isInit()) return false;
    /* check if internalFormat is valid */
    const GLint acceptedInternalFormats[25] = {
        /* fixed: fixed point number */
        /* format (bits), type of each component     */
        /* R8G8B8,        fixed8 (0.0~1.0)           */
        GL_RGB,    GL_RGB8,
        /* R8G8B8A8,      fixed8 (0.0~1.0)           */
        GL_RGBA,   GL_RGBA8,
        /* R16G16B16,     uint16 or float16          */
        GL_RGB16UI,GL_RGB16F,
        /* R16G16B16A16,  fixed16 or float16         */
        GL_RGBA16, GL_RGBA16F,
        /* R32G32B32,     int32 or uint32 or float32 */
        GL_RGB32I, GL_RGB32UI, GL_RGB32F,
        /* R32G32B32A32,  int32 or uint32 or float32 */
        GL_RGBA32I, GL_RGBA32UI, GL_RGBA32F,
        /* R8,            fixed8                     */
        GL_R8, GL_RED,
        /* R16,           fixed16, float16           */
        GL_R16, GL_R16F,
        /* R32,           uint16, float32            */
        GL_R32UI, GL_R32F,
        /* others */
        GL_DEPTH_COMPONENT, GL_DEPTH_COMPONENT16, GL_DEPTH_COMPONENT24, GL_DEPTH_COMPONENT32F,
        GL_DEPTH_STENCIL
    };
    bool internalFormatValid = false;
    int internalFormatIdx = -1;
    for (int i = 0; i < 25; i++) {
        if (internalFormat == acceptedInternalFormats[i]) {
            internalFormatValid = true;
            internalFormatIdx = i;
        }
    }
    if (internalFormatValid == false) {
        printf("error, unknown internal texture format.\n");
        return false;
    }
    this->internalFormat = internalFormat;
    /* generate an empty texture */
    glDebugCheck(glGenTextures(1, &texID));
    glDebugCheck(glBindTexture(GL_TEXTURE_2D, texID));
    if (internalFormatIdx < 20) {
        glDebugCheck(glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL));
    }
    else { /* depth component must set format to "GL_DEPTH_COMPONENT" */
        glDebugCheck(glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, w, h, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL));
    }
    glDebugCheck(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT));
    glDebugCheck(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT));
    glDebugCheck(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, sampling));
    glDebugCheck(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, sampling));
    glDebugCheck(glBindTexture(GL_TEXTURE_2D, 0));
    width = w;
    height = h;
    if (internalFormatIdx == 0 || internalFormatIdx == 1 || internalFormatIdx == 4 ||
        internalFormatIdx == 5 || internalFormatIdx == 8 || internalFormatIdx == 9 ||
        internalFormatIdx == 10) channels = 3;
    else if (internalFormatIdx == 2 || internalFormatIdx == 3 || internalFormatIdx == 6 ||
        internalFormatIdx == 7 || internalFormatIdx == 11 || internalFormatIdx == 12 ||
        internalFormatIdx == 13) channels = 4;
    else if (internalFormatIdx == 14 || internalFormatIdx == 15) channels = 1;
    else if (internalFormatIdx == 16 || internalFormatIdx == 17) channels = 2;
    else if (internalFormatIdx == 18 || internalFormatIdx == 19) channels = 4;
    else if (internalFormatIdx == 20 || internalFormatIdx == 21) channels = 1;
    else channels = -1;
    /* update mipmap if texture type is not depth or stencil */
    if (internalFormatIdx < 20) {
        updateMipmap();
    }
    return true;
}

bool GLTexture2D::gen(const char * file, GLint sampling, bool vertical_flip) {
    if (isInit()) return false;
    /* check if images should be vertically flipped when loading */
    if (vertical_flip) stbi_set_flip_vertically_on_load(true);
    else stbi_set_flip_vertically_on_load(false);
    /* load image file using stb_image */
    int w, h, c; /* image width, height, channels */
    unsigned char* image_data = stbi_load(file, &w, &h, &c, 0);
    if (image_data == NULL) {
        printf("STB_IMAGE error: %s\n", stbi_failure_reason());
        return false;
    }
    /* set internal format as RGBA8 uint */
    internalFormat = GL_RGBA;
    /* generate texture handle */
    glGenTextures(1, &texID);
    glBindTexture(GL_TEXTURE_2D, texID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, sampling);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, sampling);
    if (c == 3)
        glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, image_data);
    else if (c == 4)
        glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);
    else {
        printf("Error, unknown # of channels.\n");
        unload();
        return false;
    }
    stbi_image_free(image_data);
    width = w;
    height = h;
    channels = c;
    updateMipmap();
    return true;
}

void GLTexture2D::update(const void* newData)
{
    if (!isInit()) return;
    glBindTexture(GL_TEXTURE_2D, texID);

    glDebugCheck(
        glTexSubImage2D(GL_TEXTURE_2D, 0,
            0, 0, 
            this->width, this->height, 
            this->internalFormat, GL_UNSIGNED_BYTE, newData)
    );

    glBindTexture(GL_TEXTURE_2D, 0);
    updateMipmap();
}

void GLTexture2D::updateMipmap() {
    /* generate or update mipmap */
    if (!isInit()) return;
    glDebugCheck(glBindTexture(GL_TEXTURE_2D, texID));
    glDebugCheck(glGenerateMipmap(GL_TEXTURE_2D));
    glDebugCheck(glBindTexture(GL_TEXTURE_2D, 0));
}

void GLTexture2D::bind(unsigned int slot) {
    /* NOTE: OpenGL guarantees there are at least 16 slots at a time. */
    /* So here slot need to less or equal than 15. */
    glActiveTexture(GL_TEXTURE0 + slot);
    glBindTexture(GL_TEXTURE_2D, texID);
}

unsigned int GLTexture2D::getTextureID() { return texID; }

int GLTexture2D::getWidth() { return width; }

int GLTexture2D::getHeight() { return height; }

int GLTexture2D::getChannels() { return channels; }

void GLTexture2D::unload() {
    if (texID != 0) { glDeleteTextures(1, &texID); texID = 0; }
    width = 0; height = 0; channels = 0;
    internalFormat = GL_INVALID_VALUE;
}

bool GLTexture2D::isInit() { return texID != 0; }

GLint GLTexture2D::getInternalFormat() { 
    return internalFormat; 
}

GLShader::GLShader() { shaderProgramID = 0; }

GLShader::~GLShader() { _deleteHandlers(); }

bool GLShader::gen(const char * vertexShaderSource, const char * fragmentShaderSource) {
    /* in this function you need to create a shader */
    /* including generate handlers, compiling sources */
    /* and setting vertex attributes */

    /* If shader is already generated before, it means gen() is called twice, */
    /* maybe there are some bugs. Just return */
    if (isInit()) return false;

    /* create shader handlers */
    _createHandlers();

    /* compile shader source codes */
    int vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    int fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
    int success;
    char info[1024];
    /* compile vertex shader */
    glShaderSource(vertexShaderID, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShaderID);
    glGetShaderiv(vertexShaderID, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShaderID, 1024, NULL, info);
        printf("vertex shader compile error:\n");
        printf("%s", info);
        _deleteHandlers();
        return false;
    }
    /* compile fragment shader */
    glShaderSource(fragmentShaderID, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShaderID);
    glGetShaderiv(fragmentShaderID, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShaderID, 1024, NULL, info);
        printf("fragment shader compile error:\n");
        printf("%s", info);
        _deleteHandlers();
        return false;
    }
    /* link shaders */
    glAttachShader(shaderProgramID, vertexShaderID);
    glAttachShader(shaderProgramID, fragmentShaderID);
    glLinkProgram(shaderProgramID);
    glGetProgramiv(shaderProgramID, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgramID, 1024, NULL, info);
        printf("shader linking error:\n");
        printf("%s", info);
        _deleteHandlers();
        return false;
    }
    /* delete unused shader objects */
    glDeleteShader(vertexShaderID);
    glDeleteShader(fragmentShaderID);

    return true;
}

void GLShader::use() {
    if (shaderProgramID != 0)
        glUseProgram(shaderProgramID);
}

bool GLShader::isInit() { return shaderProgramID != 0; }

void GLShader::unload() { _deleteHandlers(); }

/* set uniform variables */

void GLShader::setUniformMatrix2f(const char * varName, const glm::mat2 & value, bool needTranspose) {
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix2fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, glm::value_ptr(value));
}

void GLShader::setUniformMatrix3f(const char * varName, const glm::mat3 & value, bool needTranspose) {
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix3fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, glm::value_ptr(value));
}

void GLShader::setUniformMatrix4f(const char * varName, const glm::mat4 & value, bool needTranspose) {
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix4fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, glm::value_ptr(value));
}

void GLShader::setUniform1f(const char * varName, float value) {
    glUniform1f(glGetUniformLocation(shaderProgramID, varName), value);
}

void GLShader::setUniform2f(const char * varName, glm::vec2 value) {
    glUniform2f(glGetUniformLocation(shaderProgramID, varName), value.x, value.y);
}

void GLShader::setUniform3f(const char * varName, glm::vec3 value) {
    glUniform3f(glGetUniformLocation(shaderProgramID, varName), value.x, value.y, value.z);
}

void GLShader::setUniform4f(const char * varName, glm::vec4 value) {
    glUniform4f(glGetUniformLocation(shaderProgramID, varName), value.x, value.y, value.z, value.w);
}

void GLShader::setUniform1i(const char * varName, int value) {
    glUniform1i(glGetUniformLocation(shaderProgramID, varName), value);
}

void GLShader::setUniform2i(const char * varName, glm::ivec2 value) {
    glUniform2i(glGetUniformLocation(shaderProgramID, varName), value.x, value.y);
}

void GLShader::setUniform3i(const char * varName, glm::ivec3 value) {
    glUniform3i(glGetUniformLocation(shaderProgramID, varName), value.x, value.y, value.z);
}

void GLShader::setUniform4i(const char * varName, glm::ivec4 value) {
    glUniform4i(glGetUniformLocation(shaderProgramID, varName), value.x, value.y, value.z, value.w);
}

void GLShader::setUniformSampler2D(const char * samplerName, GLTexture2D & texture, int slot) {
    debugCheck(slot >= 0 && slot <= 15, "invalid slot id.");
    setUniform1i(samplerName, slot);
    texture.bind(slot);
}

void GLShader::setUniform3f(const char * varName, const VEC3 & value)
{
    glUniform3f(glGetUniformLocation(shaderProgramID, varName), float(value.x), float(value.y), float(value.z));
}

void GLShader::setUniformMatrix2f(const char * varName, const float * value, bool needTranspose)
{
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix2fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, value);
}

void GLShader::setUniformMatrix3f(const char * varName, const float * value, bool needTranspose)
{
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix3fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, value);
}

void GLShader::setUniformMatrix4f(const char * varName, const float * value, bool needTranspose)
{
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix4fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, value);
}

void GLShader::setUniformMatrix3f(const char * varName, const MAT3x3 & value, bool needTranspose)
{
    float e[9];
    for (int i = 0; i < 9; i++) e[i] = float(value.e[i]);
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix4fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, e);
}

void GLShader::setUniformMatrix4f(const char * varName, const MAT4x4 & value, bool needTranspose)
{
    float e[16];
    for (int i = 0; i < 16; i++) e[i] = float(value.e[i]);
    unsigned int varLocation = glGetUniformLocation(shaderProgramID, varName);
    glUniformMatrix4fv(varLocation, 1, needTranspose ? GL_TRUE : GL_FALSE, e);
}

void GLShader::_createHandlers() {
    /* internal function to create opengl shader handlers */
    if (shaderProgramID == 0) shaderProgramID = glCreateProgram();
}

void GLShader::_deleteHandlers() {
    if (shaderProgramID != 0) {
        glDeleteProgram(shaderProgramID);
        shaderProgramID = 0;
    }
}

GLTrueTypeFont::GLTrueTypeFont() : batchSize(8), numChars(128) {
    vaoID = 0; vboID = 0;
}

GLTrueTypeFont::~GLTrueTypeFont() {
    _free();
}

bool GLTrueTypeFont::loadFont(const char * path, int pixelSize) {

    if (isLoaded()) return true;

    /* hot initialization */
    FT_Library ft;
    if (FT_Init_FreeType(&ft)) {
        printf("FreeType error: could not initialize FreeType Library.\n");
        return false;
    }

    /* load font */
    FT_Face face;
    if (FT_New_Face(ft, path, 0, &face)) {
        printf("FreeType error: failed to load font \"%s\".\n", path);
        FT_Done_FreeType(ft);
        return false;
    }

    /* set font size */
    FT_Set_Pixel_Sizes(face, 0, pixelSize);

    /* load character glyph */
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // disable byte-alignment restriction
    for (unsigned char c = 0; c < 128; c++) {
        if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
            printf("FreeType error: failed to load glyph (char value = %d).\n", int(c));
            continue;
        }
        // generate texture
        unsigned int texID;
        glGenTextures(1, &texID);
        glBindTexture(GL_TEXTURE_2D, texID);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer
        );
        // set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // now store character for later use
        SingleChar _char = {
            texID,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            (unsigned int)(face->glyph->advance.x)
        };
        allChars[c] = _char; // store this char
    }
    glBindTexture(GL_TEXTURE_2D, 0);

    /* clear resources */
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    /* initialize shaders */
    shader.gen(
        "#version 330 core                                               \n"
        "layout(location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>    \n"
        "out vec2 TexCoords;                                             \n"
        "uniform mat4 projection;                                        \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);        \n"
        "	TexCoords = vertex.zw;                                       \n"
        "}                                                               \0"
        ,
        "#version 330 core                                               \n"
        "in vec2 TexCoords;                                              \n"
        "out vec4 color;                                                 \n"
        "uniform sampler2D text;                                         \n"
        "uniform vec3 textColor;                                         \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);\n"
        "	color = vec4(textColor, 1.0) * sampled;                      \n"
        "}                                                               \0"
    );

    /* initialize buffers */
    glGenVertexArrays(1, &vaoID);
    glGenBuffers(1, &vboID);
    glBindVertexArray(vaoID);
    glBindBuffer(GL_ARRAY_BUFFER, vboID);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4 * batchSize, NULL, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

void GLTrueTypeFont::drawText(const char * text, float x, float y, float scale, glm::vec3 color) {

    debugCheck(isWindowMarginSet() == true,
        "Error, window width and height are not manually set.");

    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // set uniform values
    shader.use();
    glm::mat4 projection = glm::ortho(0.0f, float(windowWidth), 0.0f, float(windowHeight));
    shader.setUniform3f("textColor", glm::vec3(color.x, color.y, color.z));
    shader.setUniformMatrix4f("projection", projection);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(vaoID);

    // iterate through all characters

    /* speed up text rendering by reducing the call of glBufferSubData() */
    int batch_counter = 0;
    float* batch_buffer = new float[batchSize * 24];
    char* ch_buffer = new char[batchSize];

    size_t textlen = strlen(text);
    for (size_t i = 0; i < textlen; i++) {
        SingleChar ch = allChars[text[i]];
        float xpos = x + ch.bearing.x * scale;
        float ypos = y - (ch.size.y - ch.bearing.y) * scale;

        float w = ch.size.x * scale;
        float h = ch.size.y * scale;
        // update VBO for each character
        float vertices[24] = {
            xpos,     ypos + h,   0.0f, 0.0f ,
            xpos,     ypos,       0.0f, 1.0f ,
            xpos + w, ypos,       1.0f, 1.0f ,

            xpos,     ypos + h,   0.0f, 0.0f ,
            xpos + w, ypos,       1.0f, 1.0f ,
            xpos + w, ypos + h,   1.0f, 0.0f
        };
        // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
        x += (ch.advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64)

                                        // fill_ in batch buffer
        memcpy(batch_buffer + 24 * batch_counter, vertices, sizeof(float) * 24);
        ch_buffer[batch_counter] = text[i];
        batch_counter++;
        if (batch_counter == batchSize || i == textlen - 1) {
            // render glyph texture over quad if batch is full or the last 
            // remaining length of string is smaller than batch size.
            // update content of VBO memory
            glBindBuffer(GL_ARRAY_BUFFER, vboID);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * 24 * batch_counter, batch_buffer);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            for (int j = 0; j < batch_counter; j++) {
                // render quad
                glBindTexture(GL_TEXTURE_2D, allChars[ch_buffer[j]].texID);
                glDrawArrays(GL_TRIANGLES, j * 6, 6);
            }
            batch_counter = 0;
        }
    }

    delete[] batch_buffer;
    delete[] ch_buffer;

    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glDisable(GL_BLEND);
}

void GLTrueTypeFont::drawText(const char * text, REAL x, REAL y, REAL scale, const VEC3 & color)
{
    return drawText(text, float(x), float(y), float(scale), 
        glm::vec3(float(color.x), float(color.y), float(color.z)));
}

bool GLTrueTypeFont::isLoaded() {
    if (vaoID != 0 && vboID != 0) return true;
    else return false;
}

void GLTrueTypeFont::unload() { _free(); }

void GLTrueTypeFont::_free() {
    if (vboID != 0) { glDeleteBuffers(1, &vboID); vboID = 0; }
    if (vaoID != 0) { glDeleteVertexArrays(1, &vaoID); vaoID = 0; }
    for (int i = 0; i < numChars; i++) {
        unsigned int texID = allChars[i].texID;
        glDeleteTextures(1, &texID);
    }
    shader.unload();
}

void GLVertexBuffer_XY_UV::defineAttribute() {
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
}

GLTextureBlitter2D_RGBA::GLTextureBlitter2D_RGBA() { isVerticalFlip = false; isHorizontalFlip = false; }

GLTextureBlitter2D_RGBA::~GLTextureBlitter2D_RGBA() {}

void GLTextureBlitter2D_RGBA::gen() {
    /* using a fixed shader to blit texture */
    texShader.gen(
        "#version 330 core                                               \n"
        "layout(location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>    \n"
        "uniform mat4 projection;                                        \n"
        "out vec2 texCoords;                                             \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);        \n"
        "	texCoords = vertex.zw;                                       \n"
        "}                                                               \0"
        ,
        "#version 330 core                                               \n"
        "in vec2 texCoords;                                              \n"
        "uniform sampler2D tex0;                                         \n"
        "out vec4 color;                                                 \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	color = texture(tex0, texCoords);                            \n"
        "}                                                               \0"
    );
}

void GLTextureBlitter2D_RGBA::draw(GLTexture2D & texture, float x, float y, float xScale, float yScale) {
    debugCheck(isWindowMarginSet() == true,
        "Error, window width and height are not manually set.");
    if (_checkInternalFormat(texture) == false) {
        printf("error, texture internal format is not supported.\n");
        return;
    }
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    float w = float(texture.getWidth());
    float h = float(texture.getHeight());
    float vertices[24] = {
        x,              y + h * yScale,   0.0f, 1.0f ,
        x,              y,                0.0f, 0.0f ,
        x + w * xScale, y,                1.0f, 0.0f ,

        x,              y + h * yScale,   0.0f, 1.0f ,
        x + w * xScale, y,                1.0f, 0.0f ,
        x + w * xScale, y + h * yScale,   1.0f, 1.0f
    };
    if (isVerticalFlip) {
        vertices[3] = 0.0f; vertices[7] = 1.0f; vertices[11] = 1.0f;
        vertices[15] = 0.0f; vertices[19] = 1.0f; vertices[23] = 0.0f;
    }
    if (isHorizontalFlip) {
        vertices[2] = 1.0f; vertices[6] = 1.0f; vertices[10] = 0.0f;
        vertices[14] = 1.0f; vertices[18] = 0.0f; vertices[22] = 0.0f;
    }
    vbuffer.loadData(vertices, 24 * sizeof(float), GL_DYNAMIC_DRAW);
    glm::mat4 projection = glm::ortho(0.0f, float(windowWidth), 0.0f, float(windowHeight));
    texShader.use();
    texShader.setUniformMatrix4f("projection", projection);
    texture.bind(0);
    vbuffer.draw(GL_TRIANGLES, 0, 6);

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}

void GLTextureBlitter2D_RGBA::draw(GLTexture2D & texture, float x, float y) { return draw(texture, x, y, 1.0f, 1.0f); }

void GLTextureBlitter2D_RGBA::draw(GLTexture2D & texture, float x, float y, float scale) { return draw(texture, x, y, scale, scale); }

void GLTextureBlitter2D_RGBA::setVerticalFlip(bool state) {
    isVerticalFlip = state;
}

void GLTextureBlitter2D_RGBA::setHorizontalFlip(bool state) {
    isHorizontalFlip = state;
}

bool GLTextureBlitter2D_RGBA::_checkInternalFormat(GLTexture2D & texture) {
    if (texture.getInternalFormat() != GL_RGBA) return false;
    else return true;
}

void GLTextureBlitter2D_Depth::gen() {
    /* using a fixed shader to blit texture */
    texShader.gen(
        "#version 330 core                                               \n"
        "layout(location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>    \n"
        "uniform mat4 projection;                                        \n"
        "out vec2 texCoords;                                             \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);        \n"
        "	texCoords = vertex.zw;                                       \n"
        "}                                                               \0"
        ,
        "#version 330 core                                               \n"
        "in vec2 texCoords;                                              \n"
        "uniform sampler2D tex0;                                         \n"
        "out vec4 color;                                                 \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	float depthValue = texture(tex0, texCoords).r;               \n"
        "	color = vec4(vec3(depthValue), 1.0);                         \n"
        "}                                                               \0"
    );
}

bool GLTextureBlitter2D_Depth::_checkInternalFormat(GLTexture2D & texture) {
    GLint texInternalFormat = texture.getInternalFormat();
    if (texInternalFormat != GL_DEPTH_COMPONENT && texInternalFormat != GL_DEPTH_COMPONENT16 &&
        texInternalFormat != GL_DEPTH_COMPONENT24 && texInternalFormat != GL_DEPTH_COMPONENT32F &&
        texInternalFormat != GL_DEPTH_STENCIL) {
        return false;
    }
    else return true;
}

void GLTextureFilterBlur2D::gen() {
    texShader.gen(
        "#version 330 core                                               \n"
        "layout(location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>    \n"
        "uniform mat4 projection;                                        \n"
        "out vec2 texCoords;                                             \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);        \n"
        "	texCoords = vertex.zw;                                       \n"
        "}                                                               \0"
        ,
        "#version 330 core                                               \n"
        "in vec2 texCoords;                                              \n"
        "uniform sampler2D tex0;                                         \n"
        "out vec4 color;                                                 \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "   const float offset = 1.0 / 300.0;                            \n"
        "	vec2 offsets[9] = vec2[](                                    \n"
        "		vec2(-offset, offset),  // top-left                      \n"
        "		vec2(0.0f, offset),     // top-center                    \n"
        "		vec2(offset, offset),   // top-right                     \n"
        "		vec2(-offset, 0.0f),    // center-left                   \n"
        "		vec2(0.0f, 0.0f),       // center-center                 \n"
        "		vec2(offset, 0.0f),     // center-right                  \n"
        "		vec2(-offset, -offset), // bottom-left                   \n"
        "		vec2(0.0f, -offset),    // bottom-center                 \n"
        "		vec2(offset, -offset)   // bottom-right                  \n"
        "	);                                                           \n"
        "	float kernel[9] = float[](                                   \n"
        "		1, 2, 1,                                                 \n"
        "		2, 4, 2,                                                 \n"
        "		1, 2, 1                                                  \n"
        "	);                                                           \n"
        "	vec3 sampleTex[9];                                           \n"
        "	for (int i = 0; i < 9; i++){                                 \n"
        "		sampleTex[i] = vec3(texture(tex0, texCoords + offsets[i])); \n"
        "	}                                                            \n"
        "	vec3 sum = vec3(0.0);                                        \n"
        "	for (int i = 0; i < 9; i++)                                  \n"
        "	    sum += sampleTex[i] * kernel[i];                         \n"
        "	color = vec4(sum/16.0f, 1.0f);                               \n"
        "}                                                               \0"
    );
}

void GLTextureFilterSharpen2D::gen() {
    texShader.gen(
        "#version 330 core                                               \n"
        "layout(location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>    \n"
        "uniform mat4 projection;                                        \n"
        "out vec2 texCoords;                                             \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);        \n"
        "	texCoords = vertex.zw;                                       \n"
        "}                                                               \0"
        ,
        "#version 330 core                                               \n"
        "in vec2 texCoords;                                              \n"
        "uniform sampler2D tex0;                                         \n"
        "out vec4 color;                                                 \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "   const float offset = 1.0 / 300.0;                            \n"
        "	vec2 offsets[9] = vec2[](                                    \n"
        "		vec2(-offset, offset),  // top-left                      \n"
        "		vec2(0.0f, offset),     // top-center                    \n"
        "		vec2(offset, offset),   // top-right                     \n"
        "		vec2(-offset, 0.0f),    // center-left                   \n"
        "		vec2(0.0f, 0.0f),       // center-center                 \n"
        "		vec2(offset, 0.0f),     // center-right                  \n"
        "		vec2(-offset, -offset), // bottom-left                   \n"
        "		vec2(0.0f, -offset),    // bottom-center                 \n"
        "		vec2(offset, -offset)   // bottom-right                  \n"
        "	);                                                           \n"
        "	float kernel[9] = float[](                                   \n"
        "		-1, -1, -1,                                              \n"
        "		-1,  9, -1,                                              \n"
        "		-1, -1, -1                                               \n"
        "	);                                                           \n"
        "	vec3 sampleTex[9];                                           \n"
        "	for (int i = 0; i < 9; i++){                                 \n"
        "		sampleTex[i] = vec3(texture(tex0, texCoords + offsets[i])); \n"
        "	}                                                            \n"
        "	vec3 sum = vec3(0.0);                                        \n"
        "	for (int i = 0; i < 9; i++)                                  \n"
        "	    sum += sampleTex[i] * kernel[i];                         \n"
        "	color = vec4(sum, 1.0f);                                     \n"
        "}                                                               \0"
    );
}

void GLTextureFilterGammaCorr2D::gen(float gamma) {
    char fragShaderString[1024];
    sprintf(fragShaderString,
        "#version 330 core                                               \n"
        "in vec2 texCoords;                                              \n"
        "uniform sampler2D tex0;                                         \n"
        "out vec4 color;                                                 \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "   const float gammaInv = 1.0f/%.4f;                            \n"
        "	color = texture(tex0, texCoords);                            \n"
        "   color = vec4(pow(color.rgb, vec3(gammaInv)),1.0f);           \n"
        "}                                                               \0"
        , gamma);
    texShader.gen(
        "#version 330 core                                               \n"
        "layout(location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>    \n"
        "uniform mat4 projection;                                        \n"
        "out vec2 texCoords;                                             \n"
        "void main()                                                     \n"
        "{                                                               \n"
        "	gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);        \n"
        "	texCoords = vertex.zw;                                       \n"
        "}                                                               \0"
        ,
        fragShaderString);
}

GLFrameBuffer::GLFrameBuffer() {
    fboID = 0; rboID = 0; w = 0; h = 0;
    for (int i = 0; i < GRAPHICS_MAX_GL_COLOR_ATTACHMENTS; i++) attachedTextureSlots[i] = GRAPHICS_INVALID_TEXTURE_SLOT;
    hasInternalRenderBuffer = false;
}

GLFrameBuffer::~GLFrameBuffer() { unload(); }

void GLFrameBuffer::gen(int w, int h, bool genDepthBuffer) {
    /* If genDepthBuffer == true, then a write-only render buffer object  */
    /* will be created to hold the depth and stencil information, users   */
    /* will not be able to access the content of the render buffer. If    */
    /* genDepthBuffer == false, then users need to attach their own depth */
    /* buffer texture via attachTexture().                                */
    unload();

    glGenFramebuffers(1, &fboID);
    glBindFramebuffer(GL_FRAMEBUFFER, fboID);

    if (genDepthBuffer == true) {
        // create a read-only render buffer object to store depth and 
        // stencil information and attach it to frame buffer. user need 
        // to create their own color texture and attach to frame buffer
        // in order to draw something onto it.
        glGenRenderbuffers(1, &rboID);
        glBindRenderbuffer(GL_RENDERBUFFER, rboID);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, w, h);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rboID);
        this->hasInternalRenderBuffer = true;
    }
    else {
        this->hasInternalRenderBuffer = false;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    this->w = w;
    this->h = h;
}

bool GLFrameBuffer::isComplete() {
    if (fboID == 0) return false;
    // check if this framebuffer is complete
    int lastFboID = _getCurrentFrameBuffer();
    glBindFramebuffer(GL_FRAMEBUFFER, fboID);
    bool complete = (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE) ? true : false;
    glBindFramebuffer(GL_FRAMEBUFFER, lastFboID);
    return complete;
}

void GLFrameBuffer::bind() {
    debugCheck(isWindowMarginSet() == true,
        "Error, window width and height are not manually set.");
    debugCheck(_getCurrentFrameBuffer() == 0,
        "Error, cannot bind new framebuffer when current active "
        "framebuffer ID is not 0.");

    glBindFramebuffer(GL_FRAMEBUFFER, fboID);
    glViewport(0, 0, w, h);
}

void GLFrameBuffer::clearBuffers(GLbitfield mask) {
    if (isComplete() == false) return;
    int lastFboID = _getCurrentFrameBuffer();
    glBindFramebuffer(GL_FRAMEBUFFER, fboID);
    glClear(mask);
    glBindFramebuffer(GL_FRAMEBUFFER, lastFboID);
}

void GLFrameBuffer::unbind() {
    debugCheck((fboID != 0) && (_getCurrentFrameBuffer() == fboID),
        "Error, cannot unbind when framebuffer is not binded.");
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, windowWidth, windowHeight);
}

void GLFrameBuffer::attachTexture(GLTexture2D & texture, unsigned int slot) {
    if (texture.getWidth() != this->w || texture.getHeight() != this->h) {
        printf("error, cannot attach texture on a frame buffer with different size.\n");
        return;
    }
    /* users can attach two types of textures here, one is color texture, */
    /* the other is depth texture. Here we first check the texture type   */
    /* and use different strategies to deal with these different cases.   */
    /* "slot" parameter will be ignored if texture is a depth texture.    */
    GLint texInternalFormat = texture.getInternalFormat();
    if (texInternalFormat == GL_DEPTH_COMPONENT || texInternalFormat == GL_DEPTH_COMPONENT16 ||
        texInternalFormat == GL_DEPTH_COMPONENT24 || texInternalFormat == GL_DEPTH_COMPONENT32F ||
        texInternalFormat == GL_DEPTH_STENCIL) {
        /* this is a depth texture */
        /* we need to check if user has already attached any color texture */
        /* on this framebuffer before attach any depth texture. */
        if (getNumAttachedSlots() > 0) {
            printf("error, cannot attach depth texture if color texture is "
                "already attached. Please attach all color textures to frame "
                "buffer before attach any depth texture.");
            return;
        }
        int lastFboID = _getCurrentFrameBuffer();
        glBindFramebuffer(GL_FRAMEBUFFER, fboID);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, texture.getTextureID(), 0);
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);
        glBindFramebuffer(GL_FRAMEBUFFER, lastFboID);
    }
    else {
        /* is not a depth texture, consider it as a color texture */
        if (slot >= 16 || slot < 0) {
            printf("error, invalid slot value '%d' (0~15 accepted).", slot);
            return;
        }
        int lastFboID = _getCurrentFrameBuffer();
        glBindFramebuffer(GL_FRAMEBUFFER, fboID);
        /* bind texture to slot */
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + slot, GL_TEXTURE_2D, texture.getTextureID(), 0);
        glBindFramebuffer(GL_FRAMEBUFFER, lastFboID);
        for (int i = 0; i < GRAPHICS_MAX_GL_COLOR_ATTACHMENTS; i++) {
            debugCheck(attachedTextureSlots[i] != slot, "this slot has already attached a texture!");
            if (attachedTextureSlots[i] == GRAPHICS_INVALID_TEXTURE_SLOT) {
                attachedTextureSlots[i] = slot;
                break;
            }
        }
        /* activate all color attachments */
        unsigned int colorAttachments[GRAPHICS_MAX_GL_COLOR_ATTACHMENTS];
        for (int i = 0; i < getNumAttachedSlots(); i++) {
            colorAttachments[i] = GL_COLOR_ATTACHMENT0 + attachedTextureSlots[i];
        }
        glDrawBuffers(getNumAttachedSlots(), colorAttachments);
    }
}

void GLFrameBuffer::unload() {
    if (fboID != 0) { glDeleteFramebuffers(1, &fboID); fboID = 0; }
    if (rboID != 0) { glDeleteRenderbuffers(1, &rboID); rboID = 0; }
    for (int i = 0; i < GRAPHICS_MAX_GL_COLOR_ATTACHMENTS; i++) attachedTextureSlots[i] = GRAPHICS_INVALID_TEXTURE_SLOT;
}

unsigned int GLFrameBuffer::getFramebufferID() { return fboID; }

float GLFrameBuffer::getAspectRatio()
{
    return float(this->w) / float(this->h);
}

int GLFrameBuffer::getNumAttachedSlots() {
    int attachedSlots = 0;
    for (int i = 0; i < GRAPHICS_MAX_GL_COLOR_ATTACHMENTS; i++) {
        if (attachedTextureSlots[i] >= 0) attachedSlots++;
    }
    return attachedSlots;
}

GLViewport::GLViewport() { vpx = 0; vpy = 0; vpw = 0; vph = 0; }

GLViewport::~GLViewport() { unload(); }

/* initialize viewport camera in viewport (x,y,w,h) */

void GLViewport::gen(int x, int y, int w, int h) {
    texBlitter.gen();
    setViewportMargin(float(x), float(y), float(w), float(h));
    debugCheck(frameBuffer.isComplete(), "error, framebuffer is not complete!");
    clear();
}

void GLViewport::unload() {
    frameBuffer.unload();
    texViewport.unload();
}

/* resize viewport to a new size */

void GLViewport::setViewportMargin(float x, float y, float w, float h) {
    /* set viewport margin relative to main window */
    vpx = x; vpy = y; vpw = w; vph = h;
    frameBuffer.unload();
    texViewport.unload();
    frameBuffer.gen(int(w), int(h));
    texViewport.gen(int(w), int(h), GL_RGBA, GL_NEAREST);
    frameBuffer.attachTexture(texViewport);
}

void GLViewport::setViewportMargin(glm::vec4 margin) {
    return setViewportMargin(margin.x, margin.y, margin.z, margin.w);
}

glm::vec4 GLViewport::getViewportMargin() {
    return glm::vec4(vpx, vpy, vpw, vph);
}

void GLViewport::setWindowMargin(int w, int h) {
    _GLWindowMargin::setWindowMargin(w, h);
    texBlitter.setWindowMargin(w, h);
    frameBuffer.setWindowMargin(w, h);
}

void GLViewport::resizeWindow(int w, int h) { return setWindowMargin(w, h); }

void GLViewport::resizeViewport(int w, int h) { return setViewportMargin(vpx, vpy, float(w), float(h)); }

float GLViewport::getAspectRatio()
{
    return vpw / vph;
}

/* bind camera to the current active camera */

void GLViewport::bind() {
    frameBuffer.bind();
}

/* unbind framebuffer */

void GLViewport::unbind(bool update) {
    debugCheck(isWindowMarginSet() == true,
        "Error, window width and height are not manually set.");

    frameBuffer.unbind();
    if (update) {
        texBlitter.draw(texViewport, float(vpx), float(vpy));
    }
}

/* update to main window (submit changes to the viewport) */

void GLViewport::submit() {
    int currentFboID = _getCurrentFrameBuffer();
    debugCheck(currentFboID == 0 || currentFboID == frameBuffer.getFramebufferID(),
        "please bind viewport before submit.");
    if (currentFboID == 0) {
        texBlitter.draw(texViewport, float(vpx), float(vpy));
    }
    else {
        unbind(true);
        bind();
    }

}

GLCamera::GLCamera() {
    nearDistance = 0.01f; farDistance = 1e6f;
    physWidth = 1.0f; physHeight = 1.0f;
    fov = 45.0f;
    isOrthographic = false;
    pos = glm::vec3(10, 10, 10);
    look = glm::vec3(0, 0, 0);
    up = glm::vec3(0, 1, 0);
}

GLCamera::~GLCamera() {}

void GLCamera::setPosition(glm::vec3 pos) { this->pos = pos; }

void GLCamera::setLookAt(glm::vec3 look) { this->look = look; }

void GLCamera::setUp(glm::vec3 up) { this->up = up; }

void GLCamera::setCamera(glm::vec3 pos, glm::vec3 look, glm::vec3 up) {
    this->pos = pos; this->look = look; this->up = up;
}

void GLCamera::setNearFar(float nearDistance, float farDistance) {
    this->nearDistance = nearDistance; this->farDistance = farDistance;
}

/* in perspective mode, physical screen size determines the   */
/* aspect ratio of the camera, which is used when calculating */
/* the perspective projection matrix, while in orthographic   */
/* mode, physical screen size sets the bounding box of the    */
/* render area.                                               */

void GLCamera::setPhysicalScreenSize(float width, float height) {
    this->physWidth = width; this->physHeight = height;
}

/* perspective mode */

void GLCamera::setFOV(float degrees) { this->fov = fov; }

/* mode select */

void GLCamera::enableOrthographic() { this->isOrthographic = true; }

void GLCamera::disableOrthographic() { this->isOrthographic = false; }

bool GLCamera::isOrthographicEnabled() { return this->isOrthographic; }

/* calculate the projection*view matrix and this matrix can   */
/* be sent into shader to transform object vertices to device */
/* coordinates.                                               */

glm::mat4 GLCamera::getViewMatrix() {
    return glm::lookAt(pos, look, up);
}

glm::mat4 GLCamera::getProjectionMatrix() {
    glm::mat4 projection;
    if (isOrthographic) {
        projection = glm::ortho(-physWidth / 2.0f, +physWidth / 2.0f,
            -physHeight / 2.0f, +physHeight / 2.0f, nearDistance, farDistance);
    }
    else {
        projection = glm::perspective(glm::radians(fov), physWidth / physHeight,
            nearDistance, farDistance);
    }
    return projection;
}

glm::mat4 GLCamera::getViewProjectionMatrix() {
    glm::mat4 view = getViewMatrix();
    glm::mat4 projection = getProjectionMatrix();
    return projection * view;
}

void GLCamera::moveUp(float amount) {
    glm::vec3 y = glm::vec3(0, 1, 0);
    glm::vec3 diff = y * amount;
    this->pos += diff;
    this->look += diff;
}

void GLCamera::moveDown(float amount) {
    glm::vec3 y = glm::vec3(0, 1, 0);
    glm::vec3 diff = y * amount;
    this->pos -= diff;
    this->look -= diff;
}

void GLCamera::moveLeft(float amount) {
    glm::vec3 front, up, left;
    _buildLocalAxis(front, up, left);
    pos += left * amount;
    look += left * amount;
}

void GLCamera::moveRight(float amount) {
    glm::vec3 front, up, left;
    _buildLocalAxis(front, up, left);
    pos -= left * amount;
    look -= left * amount;
}

void GLCamera::moveForward(float amount) {
    glm::vec3 front, up, left;
    _buildLocalAxis(front, up, left);
    pos += front * amount;
    look += front * amount;
}

void GLCamera::moveBackward(float amount) {
    glm::vec3 front, up, left;
    _buildLocalAxis(front, up, left);
    pos -= front * amount;
    look -= front * amount;
}

void GLCamera::rotateLeft(float degrees) {
    glm::vec3 front, up, left;
    glm::vec3 y = glm::vec3(0.0f, 1.0f, 0.0f);
    _buildLocalAxis(front, up, left);
    front = glm::rotate(front, glm::radians(degrees), y);
    left = glm::rotate(left, glm::radians(degrees), y);
    up = glm::rotate(up, glm::radians(degrees), y);
    look = pos + front;
    this->up = up;
}

void GLCamera::rotateRight(float degrees) {
    glm::vec3 front, up, left;
    glm::vec3 y = glm::vec3(0, 1, 0);
    _buildLocalAxis(front, up, left);
    front = glm::rotate(front, -glm::radians(degrees), y);
    left = glm::rotate(left, -glm::radians(degrees), y);
    up = glm::rotate(up, -glm::radians(degrees), y);
    look = pos + front;
    this->up = up;
}

void GLCamera::rotateUp(float degrees) {
    glm::vec3 front, up, left;
    _buildLocalAxis(front, up, left);
    glm::vec3 y = glm::vec3(0, 1, 0);
    float cosTheta = glm::dot(front, y);
    float yDeg = glm::degrees(glm::acos(cosTheta));
    if (yDeg - degrees < 0.01f) return;
    else {
        front = glm::rotate(front, -glm::radians(degrees), left);
        up = glm::rotate(up, -glm::radians(degrees), left);
        look = pos + front;
        this->up = up;
    }
}

void GLCamera::rotateDown(float degrees) {
    glm::vec3 front, up, left;
    _buildLocalAxis(front, up, left);
    glm::vec3 y = glm::vec3(0, 1, 0);
    float cosTheta = glm::dot(front, y);
    float yDeg = glm::degrees(glm::acos(cosTheta));
    if (yDeg + degrees > 179.9f) return;
    else {
        front = glm::rotate(front, glm::radians(degrees), left);
        up = glm::rotate(up, glm::radians(degrees), left);
        look = pos + front;
        this->up = up;
    }
}

void GLCamera::_buildLocalAxis(glm::vec3 & front, glm::vec3 & up, glm::vec3 & left) {
    front = glm::normalize(look - pos);
    left = glm::normalize(glm::cross(this->up, front));
    up = glm::normalize(glm::cross(front, left));
}

VEC3 GLCamera::getPosition()
{
    return VEC3(pos.x, pos.y, pos.z);
}

VEC3 GLCamera::getLookAt()
{
    return VEC3(look.x, look.y, look.z);
}

VEC3 GLCamera::getUp()
{
    return VEC3(up.x, up.y, up.z);
}
