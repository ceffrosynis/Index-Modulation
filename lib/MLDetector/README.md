## ML Detector

The Viterbi-like Algorithm can be visualized by means of a trellis diagram

In a trellis graph, each node corresponds to a distinct state at a given time, and each arrow represents a transition to some new state at the next instant of time. We divide the process into stages. Each stage has a k-combination set of the index pattern which have a size of k and there are n stages in total. 

The image above is an example of a trellis diagram, being composed of four stages, searching for the most likelihood index pattern with length equals to four.  

![alt text](https://github.com/ceffrosynis/Index-Modulation/blob/master/images/trellis%20diagram.png)
