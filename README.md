# Summary
This repository is focused on building efficient CRUD functionality for spatial data. This is done through an interesting paradigm that indexes spatial data with insights from the concepts of polar coordinates.

The value in being able to do more efficient CRUD operations on spatial data is obvious when you consider how many use cases need to conduct CRUD operations on spatial data. Being able to do them faster would not only lead to a better user experience but also would save companies in server costs.

# The Paradigm: Polar Approach to Spatial Data
The following are the basics of how this paragigm lets us store a simple type of spatial data - points in space, as well as how points can be queried in a simple radius query. 


Representing a point in space:

  Points in space can be stored as polar coordiantes oriented around some pre defined anchor point in a two part index. The first part o the index for a point in space is an integer value R, representing the radius or distancefrom the anchor point to the points in space that is being stored. The second part of the index is a double value T representing the angle from north that the point being stored lies at. So,
                Spatial Point = (R, T)
  Together the two parts of the index for a spatial point thus represent a radius (R) and angle (T) value, or polar coordinates for the point, oriented around some pre defined anchor point. 

Querying points in space:

  If spatial points are stored as described above, efficient queries can be done on them by first filtering the integer R values and then the double T values. For example consider a radius query- one in which we are given some latitude and longitude value and asked to return the stored points within some input distance from the input ooint. Assuming we have a list of spatial points stored in polar coordinate form around some anchor A, the following steps should perform the desired radius query:
    
    1) Convert the input point into polar coordinate with an R and T value, based on the given anchor point A.
    2) Calculate the maximum and minimum R and T values that can be used to folter stored points and return only those within some given input distance of the input point.
