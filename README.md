# Summary
This repository is focused on building efficient CRUD functionality for spatial data. This is done through an interesting paradigm that indexes spatial data with insights from the concepts of polar coordinates.

The value in being able to do more efficient CRUD operations on spatial data is obvious when you consider how many use cases need to conduct CRUD operations on spatial data. Being able to do them faster would not only lead to a better user experience but also would save companies in server costs.

# The Paradigm: Polar Approach to Spatial Data
The following are the basics of how this paradigm lets us store a simple type of spatial data - points in space, as well as how points can be queried in a simple radius query. 


Representing a point in space:

  Points in space can be stored as polar coordinates oriented around some pre defined anchor point in a two part index. The first part o the index for a point in space is an integer value R, representing the radius or distance from the anchor point to the points in space that is being stored. The second part of the index is a double value T representing the angle from north that the point being stored lies at. So,
                Spatial Point = (R, T)
  Together the two parts of the index for a spatial point thus represent a radius (R) and angle (T) value, or polar coordinates for the point, oriented around some pre defined anchor point. 

Querying points in space:

  If spatial points are stored as described above, efficient queries can be done on them by first filtering the integer R values and then the double T values. For example, consider a radius query- one in which we are given some latitude and longitude value and asked to return the stored points within some input distance from the input point. Assuming we have a list of spatial points stored in polar coordinate form around some anchor A, the following steps should perform the desired radius query:

1] Convert the input point into polar coordinate with an R and T value, based on the given anchor point A.

2] Calculate the maximum and minimum R and T values that can be used to filter stored points and return only those within some given input distance of the input point.

3] First, filter the stored points by radius keeping those with a radius >= minR value && radius <= maxR value as calculated in step 3.

4] Filter these already filtered points based on the max and min T values calculated in step 3.
    
5] Return the remaining points after these two filters are done.
    
# Specs
As of right now the only working example of an algorithm based on this paradigm on this repository is a radius query on point spatial data. The information below shows the speeds at which this query can run for given input parameters, as well as some additional information on the queries. All these computations were run on my desktop Mac. 

Query Radius: 25km

Points to query: 50,000

Average execution time: 0.00087 seconds

Points returned: 45

----------------------------------------

Query radius: 25km

Points to query: 250,000

Average execution time: 0.00296 seconds

Points returned: 272

----------------------------------------

Query radius: 25km 

Points to query: 500,000

Average execution time: 0.005924 seconds

Points returned: 498


# ToDo:
The following highlights the end goal in terms of the desired scope of different efficient query types on different spatial data types we hope to support. It also highlights which combinations are supported and have existing implementations already on this repository. Obviously, the end goal is to build a general library on this repo that can be used to conduct all manner of spatial query types on all manner of spatial data types. Furthermore, we hope to support different options for creating, update and deleting spatial data using this paradigm.

_________________________________________________________________________
              |    Points   |   Lines   |   2D Shapes   |   3D Shapes
              
Radius--------|      X      |           |               |

Contains------|             |           |               |

Disjoint------|             |           |               |

Intersects----|             |           |               |

Touches-------|             |           |               |

Crosses-------|             |           |               |

Overlaps------|             |           |               |
