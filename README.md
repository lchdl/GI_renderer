# GI_renderer
A global illumination framework for rendering photorealistic images.
After running the demo you will get the following rendered image:
![final_render](https://user-images.githubusercontent.com/18594210/196041689-847f1da1-84be-4383-947b-ddd3b458b35b.png)

Implemented features:
1.  simple math library (vector math, quaternion, )
2.  robust implementation of basic template classes (`Array<T>`, `Stack<T>`, `Queue<T>`, `MaxHeap<T>`, `PriorityQueue<_storage_t, _priority_t>`)
3.  k-nn search (used in final gathering)
4.  quick sort (for building kd tree)
5.  obj mesh loader/writer
6.  simple Xorwow random number generator
7.  generic map sampler class (for bitmap and procedural textures sampling)
8.  geometry instancing (entity proxy)
9.  rigid transformations (rotation, proportional scaling, translation)
10. normal mapping
11. light cache + final gather techniques to solve the rendering equation (mainly for indirect illumination)
12. interaction between light rays and diffusive surfaces

Work in progress:
1.  refractive/reflective surface rendering
2.  procedural texture support (checkerboard, perlin noise, ...)
3.  caustics effect (may require additional steps, such as photon mapping)
