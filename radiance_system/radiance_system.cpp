#include "mesh.h"
#include "render.h"

int       g_windowWidthPx = 640;
int       g_windowHeightPx = 640;
char      g_windowTitle[1024] = "GI radiance system";
const int numKeys = 512;
bool      g_keyMap[numKeys];

InteractiveRenderer g_interactiveRenderer;

void onWindowResize(GLFWwindow* window, int width, int height);
void onKeyPressRelease(GLFWwindow* window, int key, int scancode, int action, int mods);
void formatWindowTitle(double dt, char* buffer);
void setWindowIcon(GLFWwindow* window, const char* image);

Scene* g_scene = NULL;
Camera* g_camera = NULL;
void makeScene_demo() {
    g_scene = new Scene();
    g_camera = new Camera();

    //BSTMesh mesh;
    //mesh.loadOBJ("res/model/demo/main_wall.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/main_wall.asset");
    //mesh.loadOBJ("res/model/demo/main_light.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/main_light.asset");
    //mesh.loadOBJ("res/model/demo/floor.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/floor.asset");
    //mesh.loadOBJ("res/model/demo/chair.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/chair.asset");
    //mesh.loadOBJ("res/model/demo/paint_main.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/paint_main.asset");
    //mesh.loadOBJ("res/model/demo/paint_frame_black.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/paint_frame_black.asset");
    //mesh.loadOBJ("res/model/demo/paint_frame_white.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/paint_frame_white.asset");
    //mesh.loadOBJ("res/model/demo/lamp_frame.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/lamp_frame.asset");
    //mesh.loadOBJ("res/model/demo/lamp_emissive.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/lamp_emissive.asset");
    //mesh.loadOBJ("res/model/demo/floor_2.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/floor_2.asset");
    //mesh.loadOBJ("res/model/demo/soil.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/soil.asset");
    //mesh.loadOBJ("res/model/demo/stone_01.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/stone_01.asset");
    //mesh.loadOBJ("res/model/demo/stone_02.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/stone_02.asset");
    //mesh.loadOBJ("res/model/demo/stone_03.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/stone_03.asset");
    //mesh.loadOBJ("res/model/demo/stone_04.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/stone_04.asset");
    //mesh.loadOBJ("res/model/demo/plant_01.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/plant_01.asset");
    //mesh.loadOBJ("res/model/demo/plant_02.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/plant_02.asset");
    //mesh.loadOBJ("res/model/demo/plant_03.obj", BSTMesh::defaultBuildSettings());
    //mesh.saveAsset("res/model/demo/plant_03.asset");

    int hMesh, hMaterial, hBitmap, hNormal;

    /* material for main wall */
    Material material;
    material.setDiffuse(rgbToVec3(243,213,174), -1, ONE, -1);
    hMesh = g_scene->addMeshEntity("res/model/demo/main_wall.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* stone */
    hBitmap = g_scene->loadBitmap("res/model/demo/stone_01.png");
    hNormal = g_scene->loadBitmap("res/model/demo/stone_01_normal.png");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, hNormal);
    hMesh = g_scene->addMeshEntity("res/model/demo/stone_01.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    hBitmap = g_scene->loadBitmap("res/model/demo/stone_02.png");
    hNormal = g_scene->loadBitmap("res/model/demo/stone_02_normal.png");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, hNormal);
    hMesh = g_scene->addMeshEntity("res/model/demo/stone_02.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    hBitmap = g_scene->loadBitmap("res/model/demo/stone_03.jpg");
    hNormal = g_scene->loadBitmap("res/model/demo/stone_03_normal.jpg");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, hNormal);
    hMesh = g_scene->addMeshEntity("res/model/demo/stone_03.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    hBitmap = g_scene->loadBitmap("res/model/demo/stone_04.jpeg");
    hNormal = g_scene->loadBitmap("res/model/demo/stone_04_normal.jpeg");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, hNormal);
    hMesh = g_scene->addMeshEntity("res/model/demo/stone_04.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* floor */
    hBitmap = g_scene->loadBitmap("res/model/demo/floor.jpg");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, -1);
    hMesh = g_scene->addMeshEntity("res/model/demo/floor.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);
    
    material.setDiffuse(rgbToVec3(243, 213, 174), -1, ONE, -1);
    hMaterial = g_scene->addMaterial(&material);
    hMesh = g_scene->addMeshEntity("res/model/demo/floor_2.asset");
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* soil */
    hBitmap = g_scene->loadBitmap("res/model/demo/soil.jpg");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, -1);
    hMaterial = g_scene->addMaterial(&material);
    hMesh = g_scene->addMeshEntity("res/model/demo/soil.asset");
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* chair */
    hBitmap = g_scene->loadBitmap("res/model/demo/chair.jpg");
    hNormal = g_scene->loadBitmap("res/model/demo/chair_normal.png");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, hNormal);
    hMesh = g_scene->addMeshEntity("res/model/demo/chair.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* paint */
    hBitmap = g_scene->loadBitmap("res/model/demo/paint.jpg");
    hMesh = g_scene->addMeshEntity("res/model/demo/paint_main.asset");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, -1);
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);
    hMesh = g_scene->addMeshEntity("res/model/demo/paint_frame_black.asset");
    material.setDiffuse(VEC3(REAL(0.02), REAL(0.02), REAL(0.02)), -1, ONE, -1);
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);
    hMesh = g_scene->addMeshEntity("res/model/demo/paint_frame_white.asset");
    material.setDiffuse(VEC3(REAL(0.98), REAL(0.98), REAL(0.98)), -1, ONE, -1);
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* plant */
    hBitmap = g_scene->loadBitmap("res/model/demo/plant_01.png");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, -1, true, rgbToVec3(80, 110, 60));
    hMesh = g_scene->addMeshEntity("res/model/demo/plant_01.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    hBitmap = g_scene->loadBitmap("res/model/demo/plant_02_fix.png");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, -1, true, rgbToVec3(80, 110, 60));
    hMesh = g_scene->addMeshEntity("res/model/demo/plant_02.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    hBitmap = g_scene->loadBitmap("res/model/demo/plant_03.png");
    material.setDiffuse(VEC3::ones(), hBitmap, ONE, -1, true, rgbToVec3(80, 110, 60));
    hMesh = g_scene->addMeshEntity("res/model/demo/plant_03.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* lamp */
    hMesh = g_scene->addMeshEntity("res/model/demo/lamp_frame.asset");
    material.setDiffuse(VEC3(REAL(0.02), REAL(0.02), REAL(0.02)), -1, ONE, -1);
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* main light */
    material.setEmissive(VEC3::all(80), -1, REAL(0.001), 1000);
    hMesh = g_scene->addMeshEntity("res/model/demo/main_light.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    /* lamp emissive part */
    material.setEmissive(VEC3::all(1.1), -1, REAL(0.001), 1000);
    hMesh = g_scene->addMeshEntity("res/model/demo/lamp_emissive.asset");
    hMaterial = g_scene->addMaterial(&material);
    g_scene->setEntityMaterial(hMesh, hMaterial);

    g_scene->buildScene(Scene::defaultBuildSettings());

    g_camera->pos = VEC3(0, 353, 14);
    g_camera->look = VEC3(0, -5, REAL(16.2));
    g_camera->up = VEC3(0, 0, 1);
    g_camera->fov = degToRad(15);
}

int main()
{
    GLFWwindow* window = NULL;
    {
        for (int i = 0; i < numKeys; i++)
            g_keyMap[i] = false;

        if (!glfwInit())
            return -1;

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        
        if ((window = glfwCreateWindow(g_windowWidthPx, g_windowHeightPx, g_windowTitle, NULL, NULL)) == NULL) {
            glfwTerminate();
            return -1;
        }

        glfwSetWindowSizeCallback(window, onWindowResize);
        glfwSetKeyCallback(window, onKeyPressRelease);
        glfwMakeContextCurrent(window);
        glfwWaitEventsTimeout(0.1);

        setWindowIcon(window, "res/icon32x32.png");

        GLenum err = glewInit();
        if (err != GLEW_OK) return -1;

    }
    g_interactiveRenderer.initFrontendDisplay(g_windowWidthPx, g_windowHeightPx);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    makeScene_demo();

    RenderSettings settings = RenderSettings::defaultRenderSettings();
    settings.render.w = 1200; settings.render.h = 500;
    settings.GI.knnSamples_indirectLighting = 100;
    settings.GI.resolution_directLighting_min = -4;
    settings.GI.resolution_directLighting_max = -4;
    settings.GI.resolution_indirectLighting_min = -3;
    settings.GI.resolution_indirectLighting_max = -3;
    settings.GI.knnSamples_directLighting = 40;
    settings.GI.numSamplesPerLight = 128;
    settings.GI.numFinalGatherSamples = 400;

    RenderTask_Albedo task_albedo;
    RenderTask_Normal task_normal;
    RenderTask_Depth task_depth;
    RenderTask_GI task_GI;
    task_albedo.createRenderTask(g_scene, g_camera, &settings, "res/render/albedo.png");
    task_normal.createRenderTask(g_scene, g_camera, &settings, "res/render/normal.png");
    task_depth.createRenderTask(g_scene, g_camera, &settings, "res/render/depth.png");
    task_GI.createRenderTask(g_scene, g_camera, &settings, &task_albedo, &task_normal, &task_depth);
    g_interactiveRenderer.addRenderTask(&task_albedo);
    g_interactiveRenderer.addRenderTask(&task_normal);
    g_interactiveRenderer.addRenderTask(&task_depth);
    g_interactiveRenderer.addRenderTask(&task_GI);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    int z = 0;
    char buf[64] = { 0 };

    Timer fpsTimer;
    double frameTime = 0.0;

    g_interactiveRenderer.checkRenderTasks();

    while (!glfwWindowShouldClose(window))
    {
        g_interactiveRenderer.interactiveRender();

        glfwSwapBuffers(window);
        glfwPollEvents();
        frameTime = fpsTimer.tick();
        formatWindowTitle(frameTime, g_windowTitle);
        glfwSetWindowTitle(window, g_windowTitle);
    }

    glfwTerminate();

    deleteScene();

    return 0;
}

void onWindowResize(GLFWwindow * window, int width, int height)
{
    g_windowWidthPx = width;
    g_windowHeightPx = height;
    printf("INFO: window size changed to (w=%d, h=%d).\n", width, height);
    g_interactiveRenderer.onWindowResize(width, height);
}
void onKeyPressRelease(GLFWwindow * window, int key, int scancode, int action, int mods)
{
    if (key < 0 || key >= numKeys) {
        printf("ignored message from unknown key %d.", key);
    }
    else {
        if (action == GLFW_PRESS)
            g_keyMap[key] = true;
        else if (action == GLFW_RELEASE)
            g_keyMap[key] = false;
    }
}
void formatWindowTitle(double dt, char * buffer)
{
    sprintf(buffer, "GI radiance system | ft=%.2fms, fps=%d", dt * 1000.0, int(1.0 / dt));
}

void setWindowIcon(GLFWwindow * window, const char * image)
{
    GLFWimage icons[1];
    icons[0].pixels = stbi_load(image, &icons[0].width, &icons[0].height, 0, 4);
    glfwSetWindowIcon(window, 1, icons);
    stbi_image_free(icons[0].pixels);
}
