# Bezier-Surface

Draw splines using de Castejau Algorithm.
https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm

I implement Bezier curves of order 3 in this project.

First, I made a class of Line. It draws a line consisted by two vertices.

Bezier Surface will be divided into two sub surface in respect to rows and they are divided again into two sub surface in respect to columns, again. In this project, I iterated surface dividing process 3 times.

After it finish dividing process, the surface, exactly 8 different curves, will be rendered. In this process, each curves will be divided using de Castejau Algorithm.
