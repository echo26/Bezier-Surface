# Bezier-Surface

Set the path of your VTK on the CMakeList.txt before you run the code.

In this project, I draw splines using de Castejau Algorithm. https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm
To draw splines, I implement Bezier curves of order 3 in this project using C++, VTK.

First, I made a class of Line. It draws a line that consists of two vertices.

Then, I made a class of Curve. It consists of 4 vertices. These 4 vertices will be divided into 2 different curve using de Castejau Algorithm  until they are small enough. When they are small enough, it will render 3 lines using 4 vertices.

Bezier Surface will be divided into two sub surface in respect to rows and they are divided again into two sub surface in respect to columns, again. In this project, I iterated surface dividing process 3 times.

After it finish dividing process, the surface, exactly 8 different curves, will be rendered. To be specific, 8 different curves of the surface will be rendered.
