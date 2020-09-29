===================
Cooccurrence Graphs
===================

Cooccurrence graphs are a type of network where edges are created when the 
entities represented by a pair of nodes are considered to have appeared together
in some context.

Create a Cooccurrence Graph
---------------------------

Rhodonite makes creating cooccurrence graphs easy.

The :mod:`cooccurrence_graph` function takes a nested sequence of
integer entities, and creates a :mod:`graph-tool` :mod:`Graph`, where
vertices are the entities and edges represent a cooccurrence. The 
function also returns property maps for the occurrence frequency of 
each entity and the cooccurrence frequency of each edge.

.. testcode::

   from rhodonite.cooccurrence.basic import cooccurrence_graph

   sequences = [
      [0, 1, 2],
      [1, 2, 3],
      [2, 3, 4],
   ]
   g, o, co = cooccurrence_graph(sequences)

   print('Vertex occurrence frequencies:')
   for v in g.vertices():
      print(v, o[v])

   print('\nEdge cooccurrence frequencies:')
   for e in g.edges():
      print(e, co[e])

This will return each vertex and edge, along with their respective occurrence 
and cooccurence frequencies.

.. testoutput::

   Vertex occurrence frequencies:
   0 1
   1 2
   2 3
   3 2
   4 1

   Edge cooccurrence frequencies:
   (0, 1) 1
   (0, 2) 1
   (1, 2) 2
   (1, 3) 1
   (2, 3) 2
   (2, 4) 1
   (3, 4) 1

Create a Cumulative Cooccurrence Graph
--------------------------------------

A cumulative cooccurrence graph is different from a regular cooccurrence graph 
in that it spans multiple steps in a temporal dimension. In Rhodonite, a 
cumulative cooccurrence graph is represented by a single base graph containing
all verticies, and multiple vertex and edge properties for each time step.

To create a cumulative cooccurrence graph, a set of sequences is needed for each
time step. The elements of the sequences must form a contiguous set if they
were to be flattened and sorted.

.. testcode::

   from rhodonite.cooccurrence.cumulative import cumulative_cooccurrence_graph

   sequences = [
      [
         [0, 1, 2],
         [1, 2, 3],
      ],
      [
         [0, 1, 2],
         [1, 2, 3],
         [2, 3, 4],
      ],
      [
         [0, 3, 4],
      ]
   ]

   steps = [1990, 1991, 1992]

   # in this case _c denotes cumulative properties
   g, o, o_c, co, co_c = cumulative_cooccurrence_graph(steps, sequences)

   for step in steps:
   	print('Step:', step)
   	print('{}:{:>15}{:>15}'.format('Vertex', 'Occurrence', 'Cumulative'))
    	for v in g.vertices():
        	print('{}:{:>20}{:>15}'.format(v, o[step][v], o_c[step][v]))

For cumulative cooccurrence graphs, the properties are stored in dictionaries.
The keys of the dictionaries are the input steps to the graph. In the case of
the code above, printing the occurrence values and the cumulative cooccurrence
values should yield

.. testoutput::

Step: 1990
Vertex:     Occurrence     Cumulative
0:                   1              1
1:                   2              2
2:                   2              2
3:                   1              1
4:                   0              0
Step: 1991
Vertex:     Occurrence     Cumulative
0:                   1              2
1:                   2              4
2:                   3              5
3:                   2              3
4:                   1              1
Step: 1992
Vertex:     Occurrence     Cumulative
0:                   1              3
1:                   0              4
2:                   0              5
3:                   1              4
4:                   1              2

And printing the cooccurrence and cumulative cooccurrence values for each edge
should give:

.. testoutput::

Step: 1990
Edge:       Cooccurrence     Cumulative
(0, 1):                1              1
(0, 2):                1              1
(0, 3):                0              0
(0, 4):                0              0
(1, 2):                2              2
(1, 3):                1              1
(2, 3):                1              1
(2, 4):                0              0
(3, 4):                0              0
Step: 1991
Edge:       Cooccurrence     Cumulative
(0, 1):                1              2
(0, 2):                1              2
(0, 3):                0              0
(0, 4):                0              0
(1, 2):                2              4
(1, 3):                1              2
(2, 3):                2              3
(2, 4):                1              1
(3, 4):                1              1
Step: 1992
Edge:       Cooccurrence     Cumulative
(0, 1):                0              2
(0, 2):                0              2
(0, 3):                1              1
(0, 4):                1              1
(1, 2):                0              4
(1, 3):                0              2
(2, 3):                0              3
(2, 4):                0              1
(3, 4):                1              2



