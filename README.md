# Vehicle-routing-problem
Assignment for the course [Optimization Methods in Management Science](https://www.dept.aueb.gr/en/dmst/content/optimization-methods-management-science) of the 3rd year winter semester.

## Description
This assignment requires to solve a vehicle routing problem. The business scenario is as follows:
1. The central warehouse of a company receives orders from a total number of 100 customers.
2. Trucks are based in the central warehouse.
3. The purpose is to build a route for each truck in order to achieve best customer service.
4. Each route starts starts from the main warehouse and visits various customers. Then the route ends at the last customer that the truck visited.
5. Every order must be completed by one and only visit of a truck. Therefore, when a truck visits a customer, it transfers to the customer the whole amount of his/her order.
6. Each truck has a specific product capacity, hence the goods transported by a truck must not exceed the truck's maximum capacity.
7. Assume that the vehicles travel at 35 km / hr and that for each customer, the time of unloading the goods is 15 minutes.
8. The following trucks are available in the warehouse:  
a. 15 trucks with a maximum capacity of 1500 kg.  
b. 15 trucks with a maximum capacity of 1200 kg.
9. Each route has a limit of a total duration of 3.5 ℎ𝑟.

<b> The goal is to minimize the sum of distances that all the trucks travel. </b>

## Running the application
Run the following command from the command line:

    java -jar Vehicle-routing-problem-executable.jar

## Solutions
### VND Solution
The produced VND solution can be seen with the sum distance of `1399 km`:

<img src="https://github.com/stef4k/Vehicle-routing-problem/blob/main/images/VND%20end%20solution.png" width="600" height="600" />

### VRP Solution
The produced VRP solution can be seen with the sum distance of `1476 km`

<img src="https://github.com/stef4k/Vehicle-routing-problem/blob/main/images/VRP%20end%20solution.png" width="600" height="600" />


##
[Complete Greek version of the assignment description](https://github.com/stef4k/Vehicle-routing-problem/blob/main/Greek%20Assignment%20Description.pdf)
