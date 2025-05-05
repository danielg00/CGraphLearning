## Structure of Main
	data = load_array(File)
	G = init_graph(data)
      while True
          n1, n2 = choose({Nodes})
	  if check_cyclic(n1 -> n2) == False
	      s1 = score(G(n1).add_parents(n2))
	  if check_cyclic(n2 -> n1) == False
	      s2 = score(G(n).add_parents(n1))
	  else:
	      next
	  if s1 > s2:
	      G(n1).add_parents(n2)
	  else:
	      G(n2).add_parents(n1)

	TODO:
	   Implement linked graph structure  [DONE]
	   Implement check_cycle  [DONE]
	   Implemenet add_children   [DONE]
	   
	   Implement BIC_Score (implement linear regression and calc the variance of residuals)
	   Implement delete child/parent function
	   Fix check cycle function.
	   Find memory leak in graph.c
	
	
## other TODO
* Clean up and comment linalg section
* implement score function for gaussians: Compute linear regression and then get variance of residuals.
* implement directed acyclic graph structure
* If i can implement something like a matrix class that stores data contingously, and make index function.
			
	
	
	
	
