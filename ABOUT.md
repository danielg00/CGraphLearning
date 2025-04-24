This is rudimentary implementation of a structure learning algorithm for directed graphical models in C. 

Our goal is given a dataset, we want to learn an underlying directed graph, where vertices represent random variables and arrows represent conditional dependencies.
The most simplest algorithm is a hill climbing algorithm. We start with a graph with no vertices; for each pair of vertices we use a score function to determine if drawing an arrow from one 
to the other increases our likelihood function - that is the probability that we would observe our data given the current model. If it improves the model, we add the arrow, and continue on to the next until 
convergence or some other suitable stopping criteria. 

An imporn
