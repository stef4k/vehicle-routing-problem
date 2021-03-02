
package vrp;

import java.util.ArrayList;
import java.util.Random;

public class VRP {

	double[][] distanceMatrix;
	ArrayList<Node> allNodes;
	ArrayList<Node> customers;
	Random ran;
	Node depot;
	final int numberOfCustomers = 100;
	int capacity;
	Solution bestSolutionThroughTabuSearch;

	public VRP() {
		capacity = 1500;
		int birthday = 3091997; // Stefanos birthday is on 3 September 1997
		ran = new Random(birthday);
	}

	void GenerateNetworkRandomly() {
		CreateAllNodesAndCustomerLists();
		CalculateDistanceMatrix();
	}

	public void CreateAllNodesAndCustomerLists() {
		// Create the list with the customers
		customers = new ArrayList();

		for (int i = 0; i < numberOfCustomers; i++) {
			Node cust = new Node();
			cust.x = ran.nextInt(100);
			cust.y = ran.nextInt(100);
			cust.demand = 100 * (1 + ran.nextInt(5));
			cust.serviceTime = 0.25;
			customers.add(cust);

		}
		// Build the allNodes array and the corresponding distance matrix
		allNodes = new ArrayList();
		depot = new Node();
		depot.x = 50;
		depot.y = 50;
		depot.demand = 0;
		allNodes.add(depot);
		// builds allNodes arrayList with first element the depot
		for (int i = 0; i < customers.size(); i++) {
			Node cust = customers.get(i);
			allNodes.add(cust);
		}
		// setting up the IDs of all the nodes
		for (int i = 0; i < allNodes.size(); i++) {
			Node nd = allNodes.get(i);
			nd.ID = i;
		}
	}

	public void CalculateDistanceMatrix() {

		distanceMatrix = new double[allNodes.size()][allNodes.size()];
		for (int i = 0; i < allNodes.size(); i++) {
			Node from = allNodes.get(i);

			for (int j = 0; j < allNodes.size(); j++) {
				Node to = allNodes.get(j);

				double Delta_x = (from.x - to.x);
				double Delta_y = (from.y - to.y);
				double distance = Math.sqrt((Delta_x * Delta_x) + (Delta_y * Delta_y));

				//distance = (Math.round(distance*100)) / 100;
				//distance = Math.round(distance); // different results when the distance is rounded
				distanceMatrix[i][j] = distance;
			}
		}
	}

	Solution Solve() {
		Solution s = new Solution();
		ApplyNearestNeighborMethod(s);
		TabuSearch(s);
		return s;
	}


	private void SetRoutedFlagToFalseForAllCustomers() {
		for (int i = 0; i < customers.size(); i++) {
			customers.get(i).isRouted = false;
		}
	}

	private void ApplyNearestNeighborMethod(Solution solution) {
		boolean modelIsFeasible = true;
		ArrayList<Route> routeList = solution.routes;
		SetRoutedFlagToFalseForAllCustomers();

		// As many insertions as the number of customers
		for (int insertions = 0; insertions < customers.size(); /* the insertions will be updated in the for loop */) {
			//  Insertion Identification
			CustomerInsertion bestInsertion = new CustomerInsertion();
			bestInsertion.distance = Double.MAX_VALUE;
			Route lastRoute = GetLastRoute(routeList);
			if (lastRoute != null) {
				IdentifyBestInsertion_NN(bestInsertion, lastRoute);
			}
			//Insertion Application
			// Feasible insertion was identified
			if (bestInsertion.distance < Double.MAX_VALUE) {
				ApplyCustomerInsertion(bestInsertion, solution);
				insertions++;
			} //If no insertion was feasible
			else {
				// There is a customer with demand larger than capacity -> Infeasibility
				if (lastRoute != null && lastRoute.nodes.size() == 1) {
					modelIsFeasible = false;
					break;
				} else {
					CreateAndPushAnEmptyRouteInTheSolution(solution);
				}
			}
		}

		if (modelIsFeasible == false) {
			System.out.println("There is a customer with a demand larger than the capacity");
		}
	}

	private Route GetLastRoute(ArrayList<Route> routeList) {
		if (routeList.isEmpty()) {
			return null;
		} else {
			return routeList.get(routeList.size() - 1);
		}
	}

	private void CreateAndPushAnEmptyRouteInTheSolution(Solution currentSolution) {
		if (currentSolution.routes.size() >= 15) {
			this.capacity = 1200;
		}
		if (currentSolution.routes.size() >= 30) {
			System.out.println("First solution contains more routes than allowed (31 or more)");
		}
		Route rt = new Route(capacity);
		rt.nodes.add(depot);
		currentSolution.routes.add(rt);
	}

	private void ApplyCustomerInsertion(CustomerInsertion insertion, Solution solution) {
		Node insertedCustomer = insertion.customer;
		Route route = insertion.insertionRoute;
		route.nodes.add(insertedCustomer);
		Node beforeInserted = route.nodes.get(route.nodes.size() - 2);
		double distanceAdded = distanceMatrix[beforeInserted.ID][insertedCustomer.ID];
		route.distance = route.distance + (distanceAdded);
		route.load = route.load + insertedCustomer.demand;
		double travel_duration = (distanceAdded / 35);
		route.duration = route.duration + insertedCustomer.serviceTime + travel_duration;
		solution.distance = solution.distance + distanceAdded;
		insertedCustomer.isRouted = true;
	}

	private void IdentifyBestInsertion_NN(CustomerInsertion bestInsertion, Route lastRoute) {
		for (int j = 0; j < customers.size(); j++) {
			// The examined node is called candidate
			Node candidate = customers.get(j);
			// if this candidate has not been pushed in the solution
			if (candidate.isRouted == false) {
				if (lastRoute.load + candidate.demand <= lastRoute.capacity
						&& lastRoute.duration + candidate.serviceTime
								+ (distanceMatrix[lastRoute.nodes.get(lastRoute.nodes.size() - 1).ID][candidate.ID]
										/ 35.0) <= 3.5) {
					// checks that the capacity restriction is not violated and that the duration
					// restriction is not validated
					ArrayList<Node> nodeSequence = lastRoute.nodes;
					Node lastCustomerInTheRoute = nodeSequence.get(nodeSequence.size() - 1);

					double trialDistance = distanceMatrix[lastCustomerInTheRoute.ID][candidate.ID];
					if (trialDistance < bestInsertion.distance) {
						bestInsertion.customer = candidate;
						bestInsertion.insertionRoute = lastRoute;
						bestInsertion.distance = trialDistance;
					}
				}
			}
		}
	}

	private void TabuSearch(Solution sol) {
		bestSolutionThroughTabuSearch = cloneSolution(sol);
		RelocationMove rm = new RelocationMove();
		int loop_counter = 0;
		for (int i = 0; i < 1000; i++) {
			InitializeOperators(rm); 

			FindBestRelocationMove(rm, sol);

			if (LocalOptimumHasBeenReached(rm)) {
				break;
			}

			// Apply move
			ApplyRelocationMove(rm, sol);
			//System.out.print("Iteration "+(i+1));
			//System.out.println(", solution distance:"+ " " + sol.distance +", amount of routes: " +sol.routes.size());
			loop_counter++;
		}
		//draw solution on a png  with the name of "VRP end solution"
		SolutionDrawer.drawRoutes(allNodes, sol, "VRP end solution"); 
		System.out.print("VRP relocation finished after "+loop_counter+" times (loops), contains "+sol.routes.size()+" routes");

	}

	private Solution cloneSolution(Solution sol) {
		Solution cloned = new Solution();
		// No need to clone - basic type
		cloned.distance = sol.distance;

		// Need to clone: Arraylists are objects
		for (int i = 0; i < sol.routes.size(); i++) {
			Route rt = sol.routes.get(i);
			Route clonedRoute = cloneRoute(rt);
			cloned.routes.add(clonedRoute);
		}

		return cloned;
	}

	private Route cloneRoute(Route rt) {
		Route cloned = new Route(rt.capacity);
		cloned.distance = rt.distance;
		cloned.load = rt.load;
		cloned.duration = rt.duration;
		cloned.nodes = new ArrayList();
		for (int i = 0; i < rt.nodes.size(); i++) {
			Node n = rt.nodes.get(i);
			cloned.nodes.add(n);
			// cloned.nodes.add(rt.nodes.get(i));
		}
		// cloned.nodes = rt.nodes.clone();
		return cloned;
	}

	private void InitializeOperators(RelocationMove rm) {
		rm.moveAddedDistance = Double.MAX_VALUE;
	}

	private void FindBestRelocationMove(RelocationMove rm, Solution sol) {
		ArrayList<Route> routes = sol.routes;
		for (int originRouteIndex = 0; originRouteIndex < routes.size(); originRouteIndex++) {
			Route rt1 = routes.get(originRouteIndex);
			for (int targetRouteIndex = 0; targetRouteIndex < routes.size(); targetRouteIndex++) {
				Route rt2 = routes.get(targetRouteIndex);

				for (int originNodeIndex = 1; originNodeIndex < rt1.nodes.size(); originNodeIndex++) {
					for (int targetNodeIndex = 0; targetNodeIndex < rt2.nodes.size(); targetNodeIndex++) {
						// No change for the route involved
						if (originRouteIndex == targetRouteIndex
								&& (targetNodeIndex == originNodeIndex || targetNodeIndex == originNodeIndex - 1)) {
							continue;
						}

						Node a = rt1.nodes.get(originNodeIndex - 1);
						Node b = rt1.nodes.get(originNodeIndex);
						Node c = null;
						//checking if the node we try to relocate is the
						// last node of the route
						if (originNodeIndex != rt1.nodes.size() - 1) { 
							c = rt1.nodes.get(originNodeIndex + 1); 
							// if not, then assign @c with the node after @b
						}

						Node insPoint1 = rt2.nodes.get(targetNodeIndex);
						Node insPoint2 = null;
						// checking if the target node is the last node
						// of its route
						if (targetNodeIndex != rt2.nodes.size() - 1) {
							insPoint2 = rt2.nodes.get(targetNodeIndex + 1); 
							// if not, then assign @insPoin2 with the node after @insPoint1 (in its route)
						}

						if (!check_capacity_duration(a, b, c, insPoint1, insPoint2, rt1, rt2)) {
							continue;
						}

						double moveCost = calculate_moveCost(a, b, c, insPoint1, insPoint2);
					
						StoreBestRelocationMove(originRouteIndex, targetRouteIndex, originNodeIndex, targetNodeIndex,
								moveCost, rm);
					}
				}
			}
		}
	}

	/*
	 * Method that checks the capacity and duration constraints for the specific
	 * relocation move if the constraints are violated return false
	 */
	private boolean check_capacity_duration(Node a, Node b, Node c, Node insPoint1, Node insPoint2, Route rt1,
			Route rt2) {
		// When rt1!= rt2 :No need to check the duration constraint for the
		// rt1(originRoute) because the distance (also time when velocity is stable)
		// from @a to @b added with the distance from @b to @c will always be more or
		// equal to the distance from @a to @c
		// Imagine it like a triangle of those 3 nodes. According to the triangle
		// inequity ab + bc > ac
		// equality happens only when the 3 nodes are in the same line (then no triangle
		// can be formed)
		// Also there there is no need to check the duration constraint in case that rt1
		// and rt2 are the same because
		// if the duration constraint is violated, it means that the sum distance (cost)
		// of the route has increased and as a result
		// the relocation won't happen
		if (insPoint2 == null) {
			// means that we try to relocate a node to another's route as a last destination
			if (rt1 != rt2) {
				if (rt2.load + b.demand > rt2.capacity
						|| rt2.duration + b.serviceTime + (distanceMatrix[insPoint1.ID][b.ID]) / 35 > 3.5) {
					return false;
				}
			}
		} else {// using else in order to check the constraint of duration from @insPoint1 to @b
				// and @b to @insPoint2 added with the servicetime
			if (rt1 != rt2) {
				if (rt2.load + b.demand > rt2.capacity || rt2.duration + b.serviceTime
						+ (distanceMatrix[insPoint1.ID][b.ID] + distanceMatrix[b.ID][insPoint2.ID]) / 35 > 3.5) {
					return false;
				}
			}
		}
		return true;
	}

	/*
	 * method that calculates the change of cost( change of distance of solution)
	 * according to the case
	 */
	private double calculate_moveCost(Node a, Node b, Node c, Node insPoint1, Node insPoint2) {
		double costAdded;
		double costRemoved;
		if (c == null && insPoint2 == null) {
			// means that we try to relocate a node that currently is the last destination
			// of a route
			// to another's route as a last destination
			costAdded = distanceMatrix[insPoint1.ID][b.ID];
			costRemoved = distanceMatrix[a.ID][b.ID];
		} else if (c != null && insPoint2 == null) {
			// means that we try to relocate a node that currently is NOT the last
			// destination of a route
			// to a route as a last destination)
			costAdded = distanceMatrix[a.ID][c.ID] + distanceMatrix[insPoint1.ID][b.ID];
			costRemoved = distanceMatrix[a.ID][b.ID] + distanceMatrix[b.ID][c.ID];
		} else if (c == null && insPoint2 != null) {
			// means that we try to relocate a node that currently is the last destination
			// of a route
			// to another route's but not to a last destination)
			costAdded = distanceMatrix[insPoint1.ID][b.ID] + distanceMatrix[b.ID][insPoint2.ID];
			costRemoved = distanceMatrix[a.ID][b.ID] + distanceMatrix[insPoint1.ID][insPoint2.ID];
		} else {
			// means that we try to relocate a node that currently is NOT the last
			// destination of a route
			// to another route's but NOT to a last destination)
			costAdded = distanceMatrix[a.ID][c.ID] + distanceMatrix[insPoint1.ID][b.ID]
					+ distanceMatrix[b.ID][insPoint2.ID];
			costRemoved = distanceMatrix[a.ID][b.ID] + distanceMatrix[b.ID][c.ID]
					+ distanceMatrix[insPoint1.ID][insPoint2.ID];
		}
		return (costAdded - costRemoved);
	}

	private void StoreBestRelocationMove(int originRouteIndex, int targetRouteIndex, int originNodeIndex,
			int targetNodeIndex, double moveCost, RelocationMove rm) {

		if (moveCost < rm.moveAddedDistance) {
			rm.originNodePosition = originNodeIndex;
			rm.targetNodePosition = targetNodeIndex;
			rm.targetRoutePosition = targetRouteIndex;
			rm.originRoutePosition = originRouteIndex;
			rm.moveAddedDistance = moveCost;
		}
	}

	private boolean LocalOptimumHasBeenReached(RelocationMove rm) {
		if (rm.moveAddedDistance > -0.00001) {
			return true;
		}
		return false;
	}

	private void ApplyRelocationMove(RelocationMove rm, Solution sol) {
        Route originRoute = sol.routes.get(rm.originRoutePosition);
        Route targetRoute = sol.routes.get(rm.targetRoutePosition);

        Node B = originRoute.nodes.get(rm.originNodePosition);

        if (originRoute == targetRoute) {
            originRoute.nodes.remove(rm.originNodePosition);
            if (rm.originNodePosition < rm.targetNodePosition) {
                targetRoute.nodes.add(rm.targetNodePosition, B);
            } else {
                targetRoute.nodes.add(rm.targetNodePosition + 1, B);
            }

            originRoute.distance = originRoute.distance + rm.moveAddedDistance;
            originRoute.duration = originRoute.duration + rm.moveAddedDistance/35; 
            //duration is directly correlated with the distance (v=35 km/h)
        } else {
            Node A = originRoute.nodes.get(rm.originNodePosition - 1);
            Node C = null;
            if (rm.originNodePosition != originRoute.nodes.size() - 1) {
            	C = originRoute.nodes.get(rm.originNodePosition + 1);
            }

            Node F = targetRoute.nodes.get(rm.targetNodePosition);
            Node G = null;
            if (rm.targetNodePosition != targetRoute.nodes.size() - 1) {
            	G = targetRoute.nodes.get(rm.targetNodePosition + 1);
            }
            
            
            double costChangeOrigin = calculate_costChangeOrigin(A, B, C);
            double costChangeTarget = calculate_costChangeTarget(B, F, G);
            originRoute.load = originRoute.load - B.demand;
            targetRoute.load = targetRoute.load + B.demand;

            originRoute.distance = originRoute.distance + costChangeOrigin;
            originRoute.duration = originRoute.duration - B.serviceTime + costChangeOrigin / 35;
            targetRoute.distance = targetRoute.distance + costChangeTarget;
            targetRoute.duration = targetRoute.duration + B.serviceTime + costChangeTarget / 35;

            originRoute.nodes.remove(rm.originNodePosition);
            targetRoute.nodes.add(rm.targetNodePosition + 1, B);

        }
        sol.distance = sol.distance + rm.moveAddedDistance;
    }
	
	/*
	 * method that calculate the change in distance (cost) on the originRoute which is the route that currently is the Node
	 * we are going to relocate. Returns the distance in double
	 */
	private double calculate_costChangeOrigin(Node A, Node B, Node C) {
		if(C == null) {//if B is the last node in the originRoute
			return (- distanceMatrix[A.ID][B.ID]);
		}else {
			return (distanceMatrix[A.ID][C.ID] - distanceMatrix[A.ID][B.ID] - distanceMatrix[B.ID][C.ID]);
		}
	}
	
	/*
	 * method that calculate the change in distance (cost) on the targetRoute which is the route that the Node
	 * is going to be relocated. Returns the distance in double
	 */
	private double calculate_costChangeTarget (Node B, Node F, Node G) {
		if (G == null) {//if F is the last node in the targetRoute 
			return (distanceMatrix[F.ID][B.ID]);
		}else {
			return (distanceMatrix[F.ID][B.ID] + distanceMatrix[B.ID][G.ID] - distanceMatrix[F.ID][G.ID]);
		}
	}

}
