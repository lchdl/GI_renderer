/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This header defines basic graphic tools using OpenGL, glm and FreeType. */
/*                      by: Chenghao Liu, 2021-11                          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

/* basic type defines */
#include "basedefs.h"

#include <stdio.h>
/* include OpenGL headers */
#include <glew/glew.h>  /* defines gl functions and macros */
#include <glfw/glfw3.h> /* window loop */
/* OpenGL math library */
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp> /* gtc: core, gtx: extension */
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>
/* FreeType2 library */
#include <freetype/ft2build.h>
#include FT_FREETYPE_H
/* stb image: tiny image I/O library */
#include "stb_image/stb_image_load.h"

/* define "glDebugCheck" and "debugCheck" for easier assertion  */
/* DEBUG_BREAK symbol must be correctly defined before defining */
/* the following symbols below.                                 */
#if defined(_DEBUG) || defined(DEBUG) || defined(DEBUG_MODE)
#define glDebugCheck(X)                                         \
{                                                               \
	int err;                                                    \
	glGetError(); /* clear previous error */                    \
	(X);                                                        \
	if ((err=glGetError()) != GL_NO_ERROR) {                    \
		printf("GL errno = %d\n",err);                          \
        printf("%s\n", glewGetErrorString(err));                \
		DEBUG_BREAK;                                            \
	}                                                           \
}
#else
#define glDebugCheck(X) (X);
#endif

class _GLWindowMargin {
public:
    _GLWindowMargin();
    void setWindowMargin(int w, int h);
    glm::vec2 getWindowMargin();
    bool isWindowMarginSet();
protected:
	int windowWidth, windowHeight;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* VertexBuffer : defines basic vertex buffer for drawing.         */
/*                                                                 */
/* NOTE: you must write your own derived class and implement pure  */
/*      virtual function _defineAttribute().                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GLVertexBuffer {
public:
    GLVertexBuffer();
	virtual ~GLVertexBuffer() { _deleteHandles(); }
public:
    void gen();
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* load pure vertex data                                     */
    /*    data: raw data pointer                                 */
    /*    size: buffer size in bytes                             */
    /*    usage: buffer usage, can be GL_DYNAMIC_DRAW (default), */
    /*           GL_STREAM_DRAW or GL_STATIC_DRAW                */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    void loadData(const void* data, int size, GLenum usage = GL_DYNAMIC_DRAW);
    /* load vertex data with indices */
    void loadData(const void* data, int dataSize, const void* elem, int elemSize, GLenum usage = GL_DYNAMIC_DRAW);
    /* update part of the vertex buffer, index buffer (if created) is not changed */
    void updateData(const void* data, int offset, int size);
    /* draw pure vertex data */
    void draw(GLenum mode, GLint first, GLsizei count);
    /* draw vertex data with indices */
    void draw(GLenum mode, GLsizei count, GLenum type, const void* indices);
    void unload();
    bool isInit();

protected:
	unsigned int vboID, vaoID, eboID;

protected:
	virtual void defineAttribute() = 0;

protected:
    void _createHandles();
    void _deleteHandles();
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  some predefined vertex formats and buffer objects  */
/*                                                     */
/* vertex structure that transfers between CPU and GPU */
/* since OpenGL only supports single precision float,  */
/* here we force using "float".                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct GL_VERTEX_XYZ {
    float position[3];
};
struct GL_VERTEX_XYZ_UV {
    float position[3];
    float texcoord[2];
};
struct GL_VERTEX_XYZ_RGB {
    float position[3];
    float color[3];
};
struct GL_VERTEX_XYZ_RGB_UV {
    float position[3];
    float normal[3]; /* or color */
    float texcoord[2];
};
class GLVertexBuffer_XYZ : public GLVertexBuffer {
protected:
    virtual void defineAttribute();
};
class GLVertexBuffer_XYZ_UV : public GLVertexBuffer {
protected:
    virtual void defineAttribute();
};
class GLVertexBuffer_XYZ_RGB : public GLVertexBuffer {
protected:
    virtual void defineAttribute();
};
class GLVertexBuffer_XYZ_RGB_UV : public GLVertexBuffer {
protected:
    virtual void defineAttribute();
};
class GLVertexBuffer_XY_UV : public GLVertexBuffer {
    virtual void defineAttribute() override;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Texture* : OpenGL texture wrapper                               */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GLTexture2D {

public:
    GLTexture2D();
    virtual ~GLTexture2D();

public:
    bool gen(int w, int h, GLint internalFormat, GLint sampling);
    bool gen(const char* file, GLint sampling, bool vertical_flip = true);
    void update(const void* newData); /* update whole texture */
    void updateMipmap();
    void bind(unsigned int slot);
    unsigned int getTextureID();
    int getWidth();
    int getHeight();
    int getChannels();
    void unload();
    bool isInit();
    GLint getInternalFormat();

protected:
	unsigned int texID;
	int width, height, channels;
	GLint internalFormat;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Shader : GLSL shader wrapper                                    */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GLShader {
public:
    GLShader();
    virtual ~GLShader();

public:
    bool gen(const char* vertexShaderSource, const char* fragmentShaderSource);
    void use();

    bool isInit();
    void unload();

public:
	/* set uniform variables */
    void setUniformMatrix2f(const char* varName, const glm::mat2& value, bool needTranspose = false);
    void setUniformMatrix3f(const char* varName, const glm::mat3& value, bool needTranspose = false);
    void setUniformMatrix4f(const char* varName, const glm::mat4& value, bool needTranspose = false);
    void setUniform1f(const char* varName, float value);
    void setUniform2f(const char* varName, glm::vec2 value);
    void setUniform3f(const char* varName, glm::vec3 value);
    void setUniform4f(const char* varName, glm::vec4 value);
    void setUniform1i(const char* varName, int value);
    void setUniform2i(const char* varName, glm::ivec2 value);
    void setUniform3i(const char* varName, glm::ivec3 value);
    void setUniform4i(const char* varName, glm::ivec4 value);
    void setUniformSampler2D(const char* samplerName, GLTexture2D& texture, int slot);


    void setUniform3f(const char* varName, const VEC3& value);
    void setUniformMatrix2f(const char* varName, const float* value, bool needTranspose = false);
    void setUniformMatrix3f(const char* varName, const float* value, bool needTranspose = false);
    void setUniformMatrix4f(const char* varName, const float* value, bool needTranspose = false);
    void setUniformMatrix3f(const char* varName, const MAT3x3& value, bool needTranspose = false);
    void setUniformMatrix4f(const char* varName, const MAT4x4& value, bool needTranspose = false);


protected:
    void _createHandlers();
    void _deleteHandlers();

protected:
	unsigned int shaderProgramID;

};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* TrueTypeFont : utility class to display true type fonts using   */
/*				  FreeType library.                                */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GLTrueTypeFont : public _GLWindowMargin {

	struct SingleChar{           // structure to store each character
		unsigned int texID;      // ID handle of the glyph texture
		glm::ivec2   size;       // Size of glyph
		glm::ivec2   bearing;    // Offset from baseline to left/top of glyph
		unsigned int advance;    // Offset to advance to next glyph
	};

public:

    GLTrueTypeFont();
    virtual ~GLTrueTypeFont();

public:

	/* For more details, please refer to:                 */
	/* https://learnopengl.com/In-Practice/Text-Rendering */

    bool loadFont(const char* path, int pixelSize);
    void drawText(const char* text, float x, float y, float scale, glm::vec3 color);
    void drawText(const char* text, REAL x, REAL y, REAL scale, const VEC3& color);
    bool isLoaded();
    void unload();

protected:
    void _free();

protected:
    SingleChar allChars[128];
	GLShader shader; /* shader used to print font on screen */
	unsigned int vaoID, vboID;
	const int batchSize, numChars;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* TextureBlitter: a utility class to quickly display texture onto */
/* screen. You only need one class instance of TextureBlitter to   */
/* blit all textures with each type of texture format.             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GLTextureBlitter2D_RGBA : public _GLWindowMargin {
public:
    GLTextureBlitter2D_RGBA();
    virtual ~GLTextureBlitter2D_RGBA();
public:
    virtual void gen();

    void draw(GLTexture2D& texture, float x, float y, float xScale, float yScale);
    void draw(GLTexture2D& texture, float x, float y);
    void draw(GLTexture2D& texture, float x, float y, float scale);

    void setVerticalFlip(bool state = true);
    void setHorizontalFlip(bool state = true);

protected:
    virtual bool _checkInternalFormat(GLTexture2D& texture);

protected:
	GLVertexBuffer_XY_UV vbuffer;
	GLShader texShader;
	bool isVerticalFlip, isHorizontalFlip;

};
class GLTextureBlitter2D_Depth : public GLTextureBlitter2D_RGBA {
public:
    virtual void gen();
protected:
    virtual bool _checkInternalFormat(GLTexture2D& texture);
};

/* I defined some texture filters to make some interesing effects. */
class GLTextureFilterBlur2D : public GLTextureBlitter2D_RGBA {

	/* Blur filter: apply convolution to a 2D texture. */
	/*
	kernel is : 1  2  1
	            2  4  2  *  1/16
				1  2  1
	*/
public:
    void gen();
};
class GLTextureFilterSharpen2D : public GLTextureBlitter2D_RGBA {

	/* Blur filter: apply convolution to a 2D texture. */
	/*
	kernel is : -1  -1  -1
				-1   9   1
				-1  -1  -1
	*/
public:
    void gen();
};
class GLTextureFilterGammaCorr2D : public GLTextureBlitter2D_RGBA {
	/* W correction filter */
public:
    void gen(float gamma = 2.2f);
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* FrameBuffer: OpenGL frame buffer object (FBO)                   */
/*              (simple implementation)                            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define GRAPHICS_MAX_GL_COLOR_ATTACHMENTS 16
#define GRAPHICS_INVALID_TEXTURE_SLOT (-1)
class GLFrameBuffer : public _GLWindowMargin{
public:
    GLFrameBuffer();
    virtual ~GLFrameBuffer();

public:
    /* 
    Generate a framebuffer with specific size. You need to attach 
    textures to the framebuffer to actually use it.
    */
    void gen(int w, int h, bool genDepthBuffer = true);
    /*
    Attach texture to a specific slot of a framebuffer. Each slot
    can only accept one unique texture object.
    */
    void attachTexture(GLTexture2D& texture, unsigned int slot = 0);
    /*
    Check if framebuffer is complete and ready for use.
    */
    bool isComplete();
    /*
    Bind current framebuffer to OpenGL, all subsequent draw calls 
    will use currently binded framebuffer as their render target.
    */
    void bind();
    /*
    Clear framebuffer.
    */
    void clearBuffers(GLbitfield mask = GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    /*
    Unbind current framebuffer object.
    */
    void unbind();
    /*
    Release all the resources (useful when window is resized).
    */
    void unload();
    unsigned int getFramebufferID();
    /*
    Get aspect ratio of the framebuffer (useful for calculating projection matrix).
    */
    float getAspectRatio();

protected:
	unsigned int fboID;
	unsigned int rboID;
	int w, h; /* width and height of the frame buffer (in pixels) */
	int attachedTextureSlots[GRAPHICS_MAX_GL_COLOR_ATTACHMENTS];
	bool hasInternalRenderBuffer;

    int getNumAttachedSlots();

};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Viewport : utility class to draw things in a user-defined       */
/*            viewport and display it onto screen.                 */
/*            Can be treated as a wrapper class that wraps         */
/*            glViewport() and Frame Buffer Object (FBO).          */
/*            Also provide basic screen filter effects.            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GLViewport : public _GLWindowMargin{
public:
    GLViewport();
    virtual ~GLViewport();
public:
	/* initialize viewport camera in viewport (x,y,w,h) */
    void gen(int x, int y, int w, int h);
    void unload();

	/* resize viewport to a new size */
    void setViewportMargin(float x, float y, float w, float h);
    void setViewportMargin(glm::vec4 margin);
    glm::vec4 getViewportMargin();
    void setWindowMargin(int w, int h);
    void resizeWindow(int w, int h);
    void resizeViewport(int w, int h);

    float getAspectRatio();

	/* bind camera to the current active camera */
    void bind();

	/* unbind framebuffer */
    void unbind(bool update = false);

	/* update to main window (submit changes to the viewport) */
    void submit();

	/* clear viewport */
	void clear() { frameBuffer.clearBuffers(); }

protected:
	float vpx, vpy, vpw, vph;  /* viewport x, y, w and h */

	/* each viewport camera has its own framebuffer to draw onto */
	/* also, viewport camera needs to know the main window size  */
	GLFrameBuffer frameBuffer;
	GLTexture2D texViewport;
	GLTextureBlitter2D_RGBA texBlitter;

};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Camera : Camera class that wraps the glm::perspective and       */
/*          glm::lookat call                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GLCamera {
public:
    GLCamera();
    virtual ~GLCamera();
public:
    void setPosition(glm::vec3 pos);
    void setLookAt(glm::vec3 look);
    void setUp(glm::vec3 up);
    void setCamera(glm::vec3 pos, glm::vec3 look, glm::vec3 up);
    void setNearFar(float nearDistance, float farDistance);
	/* in perspective mode, physical screen size determines the   */
	/* aspect ratio of the camera, which is used when calculating */
	/* the perspective projection matrix, while in orthographic   */
	/* mode, physical screen size sets the bounding box of the    */
	/* render area.                                               */
    void setPhysicalScreenSize(float width, float height);
	/* perspective mode */
    void setFOV(float degrees);
	/* mode select */
    void enableOrthographic();
    void disableOrthographic();
    bool isOrthographicEnabled();
	/* calculate the projection*view matrix and this matrix can   */
	/* be sent into shader to transform object vertices to device */
	/* coordinates.                                               */
    glm::mat4 getViewMatrix();
    glm::mat4 getProjectionMatrix();
    glm::mat4 getViewProjectionMatrix();

    void moveUp(float amount);
    void moveDown(float amount);
    void moveLeft(float amount);
    void moveRight(float amount);
    void moveForward(float amount);
    void moveBackward(float amount);
    void rotateLeft(float degrees);
    void rotateRight(float degrees);
    void rotateUp(float degrees);
    void rotateDown(float degrees);

protected:
    void _buildLocalAxis(glm::vec3& front, glm::vec3& up, glm::vec3& left);

protected:
	glm::vec3 pos, look, up;
	float nearDistance, farDistance;
	float physWidth, physHeight;
	float fov; // in degrees
	bool isOrthographic;

public:
    VEC3 getPosition();
    VEC3 getLookAt();
    VEC3 getUp();
};

