# Ray-tracing

Multi-level reflection using ray tracing. </br>
Lighting by phong model.

![example](https://github.com/irri094/Ray-tracing/blob/main/examples/output_11.bmp?raw=true)
![example](https://github.com/irri094/Ray-tracing/blob/main/examples/rec5.bmp?raw=true)

### Requirements:
- freeglut3-dev
- g++

### Execute:
compile main.cpp file with -lGL -lGLU -lglut. </br>
navigate with arrows, [1,6] and PageUp, PageDown.  </br>
take capture with 0.  </br>

### Scene:
#### Shapes 
* Triangle
* Sphere
* General Quadric Surfaces
* Floor </br>
Each shape has ambient, diffuse, specular, reflection coefficients and an RGB color

#### Lights 
* Point Light
* Spot Light
