package vrp;

import java.util.ArrayList;
import java.util.Random;

public class VND {
	
	double[][] distanceMatrix;
	ArrayList<Node> allNodes;
	ArrayList<Node> customers;
	Random ran;
	Node depot;
	final int numberOfCustomers = 100;
	int capacity;
	Solution bestSolutionNow;

	public VND() {
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
			// Insertion Identification
			CustomerInsertion bestInsertion = new CustomerInsertion();
			bestInsertion.distance = Double.MAX_VALUE;
			Route lastRoute = GetLastRoute(routeList);
			if (lastRoute != null) {
				IdentifyBestInsertion_NN(bestInsertion, lastRoute);
			}
			// Insertion Application
			// Feasible insertion was identified
			if (bestInsertion.distance < Double.MAX_VALUE) {
				ApplyCustomerInsertion(bestInsertion, solution);
				insertions++;
			} // If no insertion was feasible
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
	
	
	Solution solve() {
        Solution s = new Solution();
        // What would happen if Solution bestSolution = s;
        bestSolutionNow = cloneSolution(s);
        ApplyNearestNeighborMethod(s);

        RelocationMove rm = new RelocationMove();
        SwapMove sm = new SwapMove();
        TwoOptMove top = new TwoOptMove();
        
        int k = 1;
        int kmax = 3;
        int counter = 0;
        while (k <= kmax)
        {
            InitializeMoves(rm, sm, top);
            FindBestNeighbor(k, s, rm, sm, top);
            
            if (MoveIsImproving(k, rm, sm, top))
            {
                ApplyMove(k, s, rm, sm, top);
                k = 1;
            }
            else
            {
                k = k + 1;
            }
            counter++;
            //System.out.println("Solution distance:"+ " " + s.distance +", amount of routes: " +s.routes.size());
			
			
        }
        System.out.print("VND algorithm finished after " +counter +" times (loops), contains "+s.routes.size()+" routes");
        //draw solution on a png with the name of "VND end solution"
		SolutionDrawer.drawRoutes(allNodes, s, "VND end solution"); 
        
        return s;
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
	
	private  void InitializeMoves(RelocationMove rm, SwapMove sm, TwoOptMove top) {
        //Initialize the relocation move rm
        InitializeTheRelocationMove(rm);
        //Initialize the swap move sm
        InitializeTheSwapMove(sm);
        //Initialize the 2 opt move
        InitializeTheTwoOptMove(top);
    }
	
	private void ApplyMove(int k, Solution s, RelocationMove rm, SwapMove sm, TwoOptMove top) {
        
        if (k == 1)
        {
        	ApplyRelocationMove(rm, s);
        }
        else if (k == 2)
        {
            ApplySwapMove(sm, s);
        }
        if (k == 3)
        {
        	ApplyTwoOptMove(top, s);
        }
    }
	
	private  void InitializeTheRelocationMove(RelocationMove rm) {
		rm.originRoutePosition = -1;
        rm.originNodePosition = -1;
        rm.targetRoutePosition = -1;
        rm.targetNodePosition = -1;
        rm.moveAddedDistance = Double.MAX_VALUE;
    }

    private  void InitializeTheSwapMove(SwapMove sm) {
    	sm.firstRoutePosition = -1;
        sm.firstRoutePosition = -1;
        sm.secondRoutePosition = -1;
        sm.secondNodePosition = -1;
        sm.moveAddedDistance = Double.MAX_VALUE;
    }

    private  void InitializeTheTwoOptMove(TwoOptMove top) 
    {
        top.positionOfFirstRoute = -1;
        top.positionOfSecondRoute = -1;
        top.positionOfFirstNode = -1;
        top.positionOfSecondNode = -1;
        top.moveAddedDistance = Double.MAX_VALUE;
    }
    
    private  void FindBestNeighbor(int k, Solution s, RelocationMove rm, SwapMove sm, TwoOptMove top) 
    {
        if (k == 1)
        {
        	FindBestRelocationMove(rm, s);
        }
        else if (k == 2)
        {
        	FindBestSwapMove(sm, s);
        }
        else if (k == 3)
        {
        	FindBestTwoOptMove(top, s);
        }
    }
    
    private void FindBestRelocationMove(RelocationMove rm, Solution sol) {
		ArrayList<Route> routes = sol.routes;
		for (int originRouteIndex = 0; originRouteIndex < routes.size(); originRouteIndex++) {
			Route rt1 = routes.get(originRouteIndex);
			for (int targetRouteIndex = 0; targetRouteIndex < routes.size(); targetRouteIndex++) {
				Route rt2 = routes.get(targetRouteIndex);

				for (int originNodeIndex = 1; originNodeIndex < rt1.nodes.size(); originNodeIndex++) {
					for (int targetNodeIndex = 0; targetNodeIndex < rt2.nodes.size(); targetNodeIndex++) {
						// Why? No change for the route involved
						if (originRouteIndex == targetRouteIndex
								&& (targetNodeIndex == originNodeIndex || targetNodeIndex == originNodeIndex - 1)) {
							continue;
						}
						Node a = rt1.nodes.get(originNodeIndex - 1);
						Node b = rt1.nodes.get(originNodeIndex);
						Node c = null;
						// checking if the node we try to relocate is the
						// last node of the route
						if (originNodeIndex != rt1.nodes.size() - 1) {
							c = rt1.nodes.get(originNodeIndex + 1); // if not, then assign @c with the node after @b
						}

						Node insPoint1 = rt2.nodes.get(targetNodeIndex);
						Node insPoint2 = null;
						// checking if the target node is the last node
						// of its route
						if (targetNodeIndex != rt2.nodes.size() - 1) { 
							insPoint2 = rt2.nodes.get(targetNodeIndex + 1); 
							// if not, then assign @insPoin2 with the node after @insPoint1 (in its route)
						}

						if (!checkRelocation_capacity_duration(a, b, c, insPoint1, insPoint2, rt1, rt2)) {
							continue;
						}

						double moveCost = relocation_calculate_moveCost(a, b, c, insPoint1, insPoint2);

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
	private boolean checkRelocation_capacity_duration(Node a, Node b, Node c, Node insPoint1, Node insPoint2, Route rt1,
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
	private double relocation_calculate_moveCost(Node a, Node b, Node c, Node insPoint1, Node insPoint2) {
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
            
            
            double costChangeOrigin = calculate_swap_costChangeOrigin(A, B, C);
            double costChangeTarget = calculate_swap_costChangeTarget(B, F, G);
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
	private double calculate_swap_costChangeOrigin(Node A, Node B, Node C) {
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
	private double calculate_swap_costChangeTarget (Node B, Node F, Node G) {
		if (G == null) {//if F is the last node in the targetRoute 
			return (distanceMatrix[F.ID][B.ID]);
		}else {
			return (distanceMatrix[F.ID][B.ID] + distanceMatrix[B.ID][G.ID] - distanceMatrix[F.ID][G.ID]);
		}
	}
	
	private boolean MoveIsImproving(int k, RelocationMove rm, SwapMove sm, TwoOptMove top) 
    {
        if (k == 1)
        {
            if (rm.moveAddedDistance < -0.00001)
            {
                return true;
            }
        }
        else if (k == 2)
        {
            if (sm.moveAddedDistance < -0.00001)
            {
                return true;
            }
        }
        else if (k == 3)
        {
            if (top.moveAddedDistance < -0.00001)
            {
                return true;
            }
        }
        
        return false;
    }
	
	
	private void FindBestSwapMove(SwapMove sm, Solution sol) {
        ArrayList<Route> routes = sol.routes;
        for (int firstRouteIndex = 0; firstRouteIndex < routes.size(); firstRouteIndex++) {
            Route rt1 = routes.get(firstRouteIndex);
            for (int secondRouteIndex = firstRouteIndex; secondRouteIndex < routes.size(); secondRouteIndex++) {
                Route rt2 = routes.get(secondRouteIndex);
                for (int firstNodeIndex = 1; firstNodeIndex < rt1.nodes.size(); firstNodeIndex++) {
                    int startOfSecondNodeIndex = 1;
                    if (rt1 == rt2) {
                        startOfSecondNodeIndex = firstNodeIndex + 1;
                    }
                    for (int secondNodeIndex = startOfSecondNodeIndex; secondNodeIndex < rt2.nodes.size() ; secondNodeIndex++) {
                        Node a1 = rt1.nodes.get(firstNodeIndex - 1);
                        Node b1 = rt1.nodes.get(firstNodeIndex);
                        Node c1 = null;
                        // checking if the node we try to swap is the
						// last node of the route
						if (firstNodeIndex != rt1.nodes.size() - 1) { 
							c1 = rt1.nodes.get(firstNodeIndex + 1); 
							// if not, then assign @c1 with the node after @b1
						}

                        Node a2 = rt2.nodes.get(secondNodeIndex - 1);
                        Node b2 = rt2.nodes.get(secondNodeIndex);
                        Node c2 = null;
                        // checking if the second node we try to swap is the
						// last node of its route
                        if (secondNodeIndex != rt2.nodes.size() - 1) { 
                        	c2 = rt2.nodes.get(secondNodeIndex + 1); 
                        	// if not, then assign @c2 with the node after @b2
						}
                        
                        double moveCost = calculate_swapMoveCost(a1, b1, c1, a2, b2, c2, rt1, rt2);
                        StoreBestSwapMove(firstRouteIndex, secondRouteIndex, firstNodeIndex,
                        		secondNodeIndex, moveCost, sm);
                    }
                }
            }
        }
    }
	
	/*
	 * method that calculates the change of cost( change of distance of solution)
	 * according to the case
	 */
	private double calculate_swapMoveCost (Node a1, Node b1, Node c1, Node a2, Node b2, Node c2, Route rt1, Route rt2){
		double moveCost = Double.MAX_VALUE;
        if (rt1 == rt2) // within route 
        {
            if (c1.equals(b2)) { //if the nodes are next to each other
            	double costRemoved;
        		double costAdded;
            	if(c2 == null) { //if b2 is the last node of the route
            		costRemoved = distanceMatrix[a1.ID][b1.ID] + distanceMatrix[b1.ID][b2.ID];
            		costAdded = distanceMatrix[a1.ID][b2.ID] + distanceMatrix[b2.ID][b1.ID];
            	}else {
            		costRemoved = distanceMatrix[a1.ID][b1.ID] + distanceMatrix[b1.ID][b2.ID] 
            				+ distanceMatrix[b2.ID][c2.ID];
            		costAdded = distanceMatrix[a1.ID][b2.ID] + distanceMatrix[b2.ID][b1.ID] 
            				+ distanceMatrix[b1.ID][c2.ID];
            	}
            moveCost = costAdded - costRemoved;
            } else {
            	double costRemoved1;
            	double costAdded1;
            	double costRemoved2;
            	double costAdded2;
            	if (c2 == null) { //if b2 is the last node of the route
            		costRemoved1 = distanceMatrix[a1.ID][b1.ID] + distanceMatrix[b1.ID][c1.ID];
                    costAdded1 = distanceMatrix[a1.ID][b2.ID] + distanceMatrix[b2.ID][c1.ID];
                    costRemoved2 = distanceMatrix[a2.ID][b2.ID] ;
                    costAdded2 = distanceMatrix[a2.ID][b1.ID];
            	}else {
	                costRemoved1 = distanceMatrix[a1.ID][b1.ID] + distanceMatrix[b1.ID][c1.ID];
	                costAdded1 = distanceMatrix[a1.ID][b2.ID] + distanceMatrix[b2.ID][c1.ID];
	                costRemoved2 = distanceMatrix[a2.ID][b2.ID] + distanceMatrix[b2.ID][c2.ID];
	                costAdded2 = distanceMatrix[a2.ID][b1.ID] + distanceMatrix[b1.ID][c2.ID];
            	}
            	moveCost = costAdded1 + costAdded2 - (costRemoved1 + costRemoved2);
            }
        // no need to check for duration constraint because if the distance of the route increases 
        //(danger of violating  the constraint) the solution won't be the stored 
        } else // between routes
        {
            //capacity constraints
            if (rt1.load - b1.demand + b2.demand > rt1.capacity) {
            	return moveCost;
            }
            if (rt2.load - b2.demand + b1.demand > rt2.capacity) {
            	return moveCost;
            }
            
            double costRemoved1;
            double costAdded1;
            double costRemoved2;
            double costAdded2;
            if (c1 == null && c2 == null) {//means that both @b1 and @b2 are the last nodes on their own routes
            	costRemoved1 = distanceMatrix[a1.ID][b1.ID];
                costAdded1 = distanceMatrix[a1.ID][b2.ID];
                costRemoved2 = distanceMatrix[a2.ID][b2.ID];
                costAdded2 = distanceMatrix[a2.ID][b1.ID];
            }else if (c1 != null && c2 == null) { //means that @b2 is the last node on its route
            	costRemoved1 = distanceMatrix[a1.ID][b1.ID] + distanceMatrix[b1.ID][c1.ID];
                costAdded1 = distanceMatrix[a1.ID][b2.ID] + distanceMatrix[b2.ID][c1.ID];
                costRemoved2 = distanceMatrix[a2.ID][b2.ID];
                costAdded2 = distanceMatrix[a2.ID][b1.ID];
            }else if (c1 == null && c2 != null) {//means that @b1 is the last node on its route
            	costRemoved1 = distanceMatrix[a1.ID][b1.ID];
                costAdded1 = distanceMatrix[a1.ID][b2.ID];
                costRemoved2 = distanceMatrix[a2.ID][b2.ID] + distanceMatrix[b2.ID][c2.ID];
                costAdded2 = distanceMatrix[a2.ID][b1.ID] + distanceMatrix[b1.ID][c2.ID];
            }else {
	            costRemoved1 = distanceMatrix[a1.ID][b1.ID] + distanceMatrix[b1.ID][c1.ID];
	            costAdded1 = distanceMatrix[a1.ID][b2.ID] + distanceMatrix[b2.ID][c1.ID];
	            costRemoved2 = distanceMatrix[a2.ID][b2.ID] + distanceMatrix[b2.ID][c2.ID];
	            costAdded2 = distanceMatrix[a2.ID][b1.ID] + distanceMatrix[b1.ID][c2.ID];
            }
            //duration constraint 
            if (rt1.duration - b1.serviceTime + b2.serviceTime + (costAdded1 - costRemoved1)/35 > 3.5) {
            	return moveCost;
            }
            if (rt2.duration - b2.serviceTime + b1.serviceTime + (costAdded2 - costRemoved2)/35 > 3.5) {
            	return moveCost;
            }

            moveCost = costAdded1 + costAdded2 - (costRemoved1 + costRemoved2);
        }
        return moveCost;
	}
	
	private void StoreBestSwapMove(int firstRouteIndex, int secondRouteIndex, int firstNodeIndex, 
			int secondNodeIndex, double moveCost, SwapMove sm) {
        if (moveCost < sm.moveAddedDistance) {
            sm.firstRoutePosition = firstRouteIndex;
            sm.firstNodePosition = firstNodeIndex;
            sm.secondRoutePosition = secondRouteIndex;
            sm.secondNodePosition = secondNodeIndex;
            sm.moveAddedDistance = moveCost;
        }
    }
	
	private void ApplySwapMove(SwapMove sm, Solution sol) {
        Route firstRoute = sol.routes.get(sm.firstRoutePosition);
        Route secondRoute = sol.routes.get(sm.secondRoutePosition);

        if (firstRoute == secondRoute) {
            if (sm.firstNodePosition == sm.secondNodePosition - 1) {
                Node A = firstRoute.nodes.get(sm.firstNodePosition);
                Node B = firstRoute.nodes.get(sm.firstNodePosition + 1);

                firstRoute.nodes.set(sm.firstNodePosition, B);
                firstRoute.nodes.set(sm.firstNodePosition + 1, A);

            } else {
                Node A = firstRoute.nodes.get(sm.firstNodePosition);
                Node B = firstRoute.nodes.get(sm.secondNodePosition);

                firstRoute.nodes.set(sm.firstNodePosition, B);
                firstRoute.nodes.set(sm.secondNodePosition, A);
            }
            firstRoute.distance = firstRoute.distance + sm.moveAddedDistance;
            firstRoute.duration = firstRoute.duration + sm.moveAddedDistance / 35;
        } else {
            Node A = firstRoute.nodes.get(sm.firstNodePosition - 1);
            Node B = firstRoute.nodes.get(sm.firstNodePosition);
            Node C = null;
            if (sm.firstNodePosition != firstRoute.nodes.size() - 1) {
            	C = firstRoute.nodes.get(sm.firstNodePosition + 1);
        	}

            Node E = secondRoute.nodes.get(sm.secondNodePosition - 1);
            Node F = secondRoute.nodes.get(sm.secondNodePosition);
            Node G = null;
            if (sm.secondNodePosition != secondRoute.nodes.size() - 1) {
            	G = secondRoute.nodes.get(sm.secondNodePosition + 1);
            }

            double costChangeFirstRoute = calculate_swapCostChangeFirstRoute (A, B, C, F);
            double costChangeSecondRoute = calculate_swapcostChangeSecondRoute(B, E, F, G);

            firstRoute.distance = firstRoute.distance + costChangeFirstRoute;
            secondRoute.distance = secondRoute.distance + costChangeSecondRoute;

            firstRoute.load = firstRoute.load + F.demand - B.demand;
            secondRoute.load = secondRoute.load + B.demand - F.demand;
            
            firstRoute.duration = firstRoute.duration + F.serviceTime - B.serviceTime + costChangeFirstRoute / 35;
            secondRoute.duration = secondRoute.duration + B.serviceTime - F.serviceTime + costChangeSecondRoute / 35;

            firstRoute.nodes.set(sm.firstNodePosition, F);
            secondRoute.nodes.set(sm.secondNodePosition, B);

        }
        sol.distance = sol.distance + sm.moveAddedDistance;
    }
	
	
	/*
	 * method only for swapMove, calculates the change in distance (cost) on the firstRoute which is the route
	 *  where the first Node is that will be swapped after. Returns the distance in double
	 */
	private double calculate_swapCostChangeFirstRoute( Node A, Node B, Node C, Node F) {
		if (C == null) {//if B is the last node of the route
			return (distanceMatrix[A.ID][F.ID] - distanceMatrix[A.ID][B.ID]);
		}else {
			return (distanceMatrix[A.ID][F.ID] + distanceMatrix[F.ID][C.ID] - distanceMatrix[A.ID][B.ID]
					- distanceMatrix[B.ID][C.ID]);
		}
	}
	
	/*
	 * method only for swapMove, calculates the change in distance (cost) on the secondRoute which is the route 
	 * where the second Node is that will be swapped after. Returns the distance in double
	 */
	private double calculate_swapcostChangeSecondRoute(Node B, Node E, Node F, Node G) {
		if (G == null) { //if F is the last node of the route
			return (distanceMatrix[E.ID][B.ID] - distanceMatrix[E.ID][F.ID]);
		}else {
			return (distanceMatrix[E.ID][B.ID] + distanceMatrix[B.ID][G.ID] - distanceMatrix[E.ID][F.ID] 
					- distanceMatrix[F.ID][G.ID]);
		}
	}
	
	private void FindBestTwoOptMove(TwoOptMove top, Solution sol) {
        for (int rtInd1 = 0; rtInd1 < sol.routes.size(); rtInd1++) {
            Route rt1 = sol.routes.get(rtInd1);

            for (int rtInd2 = rtInd1; rtInd2 < sol.routes.size(); rtInd2++) {
                Route rt2 = sol.routes.get(rtInd2);

                for (int nodeInd1 = 0; nodeInd1 < rt1.nodes.size(); nodeInd1++) {
                    int start2 = 0;
                    if (rt1 == rt2) {
                        start2 = nodeInd1 + 2;
                    }

                    for (int nodeInd2 = start2; nodeInd2 < rt2.nodes.size(); nodeInd2++) 
                    {
                        double moveCost = Double.MAX_VALUE;

                        if (rt1 == rt2) {
                            Node A = rt1.nodes.get(nodeInd1);
                            Node B = rt1.nodes.get(nodeInd1 + 1);
                            Node K = rt2.nodes.get(nodeInd2);
                            Node L = null;
                            //checking if K is the last node of the route
                            if ( nodeInd2 != rt2.nodes.size() - 1) { 
                            	L = rt2.nodes.get(nodeInd2 + 1);
                            }

                            if (nodeInd1 == 0 && nodeInd2 == rt1.nodes.size() - 1) {
                                continue;
                            }
                            
                            double costAdded;
                            double costRemoved;
                            //if K is the last node of the route
                            if (L == null) {
                            	costAdded = distanceMatrix[A.ID][K.ID];
                                costRemoved = distanceMatrix[A.ID][B.ID];
                            }else {
                            costAdded = distanceMatrix[A.ID][K.ID] + distanceMatrix[B.ID][L.ID];
                            costRemoved = distanceMatrix[A.ID][B.ID] + distanceMatrix[K.ID][L.ID];
                            }

                            moveCost = costAdded - costRemoved;

                        } else {
                            Node A = (rt1.nodes.get(nodeInd1));
                            Node B = null;
                            if ( nodeInd1 != rt1.nodes.size() - 1) { 
                            	//checking if B is the last node of the route
                            	B = (rt1.nodes.get(nodeInd1 + 1));
                            }
                            Node K = (rt2.nodes.get(nodeInd2));
                            Node L = null;
                            if ( nodeInd2 != rt2.nodes.size() - 1) {
                            	//checking if K is the last node of the route
                            	L = rt2.nodes.get(nodeInd2 + 1);
                            }

                            if (nodeInd1 == 0 && nodeInd2 == 0) {
                                continue;
                            }
                            if (nodeInd1 == rt1.nodes.size() - 1 && nodeInd2 == rt2.nodes.size() - 1) {
                                continue;
                            }

                            if (TwoOptConstraintsAreViolated(rt1, nodeInd1, rt2, nodeInd2)) {
                                continue;
                            }
                            
                            double costAdded;
                            double costRemoved;
                            if ( B != null && L == null) { //in no case is both B and L null
                            	costAdded = distanceMatrix[B.ID][K.ID];
                                costRemoved = distanceMatrix[A.ID][B.ID];
                            }else if ( B == null && L != null) {
                            	costAdded = distanceMatrix[A.ID][L.ID];
                            	costRemoved = distanceMatrix[K.ID][L.ID];
                            }else {
                            	costAdded = distanceMatrix[A.ID][L.ID] + distanceMatrix[B.ID][K.ID];
                                costRemoved = distanceMatrix[A.ID][B.ID] + distanceMatrix[K.ID][L.ID];
                            }

                            moveCost = costAdded - costRemoved;
                        }

                        if (moveCost < top.moveAddedDistance) 
                        {
                            StoreBestTwoOptMove(rtInd1, rtInd2, nodeInd1, nodeInd2, moveCost, top);
                        }
                    }
                }
            }
        }
    }

    private void StoreBestTwoOptMove(int rtInd1, int rtInd2, int nodeInd1, int nodeInd2, double moveCost, TwoOptMove top) {
        top.positionOfFirstRoute = rtInd1;
        top.positionOfSecondRoute = rtInd2;
        top.positionOfFirstNode = nodeInd1;
        top.positionOfSecondNode = nodeInd2;
        top.moveAddedDistance = moveCost;
    }
    
    private void ApplyTwoOptMove(TwoOptMove top, Solution sol) 
    {
        Route rt1 = sol.routes.get(top.positionOfFirstRoute);
        Route rt2 = sol.routes.get(top.positionOfSecondRoute);

        if (rt1 == rt2) 
        {
            ArrayList modifiedRt = new ArrayList();

            for (int i = 0; i <= top.positionOfFirstNode; i++) 
            {
                modifiedRt.add(rt1.nodes.get(i));
            }
            for (int i = top.positionOfSecondNode; i > top.positionOfFirstNode; i--) 
            {
                modifiedRt.add(rt1.nodes.get(i));
            }
            for (int i = top.positionOfSecondNode + 1; i < rt1.nodes.size(); i++) 
            {
                modifiedRt.add(rt1.nodes.get(i));
            }

            rt1.nodes = modifiedRt;
            
            rt1.distance += top.moveAddedDistance;
            rt1.duration += top.moveAddedDistance /35;
            sol.distance += top.moveAddedDistance;
        }
        else
        {
            ArrayList modifiedRt1 = new ArrayList();
            ArrayList modifiedRt2 = new ArrayList();
            
            Node A = (rt1.nodes.get(top.positionOfFirstNode));
            Node B = null;
            if (top.positionOfFirstNode != rt1.nodes.size() - 1 ) {
            	B = (rt1.nodes.get(top.positionOfFirstNode + 1));
            }
            Node K = (rt2.nodes.get(top.positionOfSecondNode));
            Node L = null;
            if (top.positionOfSecondNode != rt2.nodes.size() - 1) {
            	L = (rt2.nodes.get(top.positionOfSecondNode + 1));
            }
            
           
            for (int i = 0 ; i <= top.positionOfFirstNode; i++)
            {
                modifiedRt1.add(rt1.nodes.get(i));
            }
             for (int i = top.positionOfSecondNode + 1 ; i < rt2.nodes.size(); i++)
            {
                modifiedRt1.add(rt2.nodes.get(i));
            }
             
            for (int i = 0 ; i <= top.positionOfSecondNode; i++)
            {
                modifiedRt2.add(rt2.nodes.get(i));
            }
            for (int i = top.positionOfFirstNode + 1 ; i < rt1.nodes.size(); i++)
            {
                modifiedRt2.add(rt1.nodes.get(i));
            }
            
            double rt1SegmentLoad = 0;
            for (int i = 0 ; i <= top.positionOfFirstNode; i++)
            {
                rt1SegmentLoad += rt1.nodes.get(i).demand;
            }
            
            double rt2SegmentLoad = 0;
            for (int i = 0 ; i <= top.positionOfSecondNode; i++)
            {
                rt2SegmentLoad += rt2.nodes.get(i).demand;
            }
            
            double originalRt1Load = rt1.load;
            rt1.load = rt1SegmentLoad + (rt2.load - rt2SegmentLoad);
            rt2.load = rt2SegmentLoad + (originalRt1Load - rt1SegmentLoad);
            
            rt1.nodes = modifiedRt1;
            rt2.nodes = modifiedRt2;
            
            rt1.distance = UpdateRouteCost(rt1);
            rt2.distance = UpdateRouteCost(rt2);
            
            rt1.duration = UpdateRouteDuration(rt1);
            rt2.duration = UpdateRouteDuration(rt2);
            
            sol.distance += top.moveAddedDistance;
        }
    }
    
    private double UpdateRouteDuration(Route rt) 
    {
    	double totDistance = 0 ;
    	double totServiceTime = 0;
        for (int i = 0 ; i < rt.nodes.size()-1; i++)
        {
            Node A = rt.nodes.get(i);
            Node B = rt.nodes.get(i+1);
            totDistance += distanceMatrix[A.ID][B.ID];
            totServiceTime += B.serviceTime;
        }
        return (totServiceTime + totDistance / 35);
    }
    
    private double UpdateRouteCost(Route rt) 
    {
        double totCost = 0 ;
        for (int i = 0 ; i < rt.nodes.size()-1; i++)
        {
            Node A = rt.nodes.get(i);
            Node B = rt.nodes.get(i+1);
            totCost += distanceMatrix[A.ID][B.ID];
        }
        return totCost;
    }
    
    /*
     * Method that checks if a constraint in the two opt move is violated.
     * The constrains are capacity of a route and the duration of a route
     * If a constraint is 
     * violated the method returns true, otherwise it returns false
     */
    private boolean TwoOptConstraintsAreViolated(Route rt1, int nodeInd1, Route rt2, int nodeInd2) 
    {
	    double rt1FirstSegmentLoad = 0;
	    double rt1FirstSegmentDuration = 0;
	    for (int i = 0 ; i <= nodeInd1; i++)
	    {
	        rt1FirstSegmentLoad += rt1.nodes.get(i).demand;
	    }
	    for (int i = 0 ; i < nodeInd1; i++)
	    {
	    	rt1FirstSegmentDuration = rt1FirstSegmentDuration + rt1.nodes.get(i + 1).serviceTime 
	    			+ distanceMatrix[rt1.nodes.get(i).ID][rt1.nodes.get(i + 1).ID] / 35;
	    }
	    double rt1SecondSegmentLoad = rt1.load - rt1FirstSegmentLoad;
	    double rt1SecondSegmentDuration = 0;
	    // if the first node of the two opt move is not the last node of its route
	    if (nodeInd1 != rt1.nodes.size() - 1) {
	    	rt1SecondSegmentDuration = rt1.duration - rt1FirstSegmentDuration 
	    			- distanceMatrix[rt1.nodes.get(nodeInd1).ID][rt1.nodes.get(nodeInd1 + 1).ID] / 35;
	    }else {
	    	rt1SecondSegmentDuration = rt1.duration - rt1FirstSegmentDuration;
	    }
	    double rt2FirstSegmentLoad = 0;
	    double rt2FirstSegmentDuration = 0;
	    for (int i = 0 ; i <= nodeInd2; i++)
	    {
	        rt2FirstSegmentLoad += rt2.nodes.get(i).demand;
	    }
	    for (int i = 0 ; i < nodeInd2; i++)
	    {
	    	rt2FirstSegmentDuration = rt2FirstSegmentDuration + rt2.nodes.get(i + 1).serviceTime 
	    			+ distanceMatrix[rt2.nodes.get(i).ID][rt2.nodes.get(i + 1).ID] / 35;
	    }
	    double rt2SecondSegmentLoad = rt2.load - rt2FirstSegmentLoad;
	    
	    double rt2SecondSegmentDuration = 0;
	    // if the second node of the two opt move is not the last node of its route
	    if (nodeInd2 != rt2.nodes.size() - 1) {
	    	rt2SecondSegmentDuration = rt2.duration - rt2FirstSegmentDuration 
	    			- distanceMatrix[rt2.nodes.get(nodeInd2).ID][rt2.nodes.get(nodeInd2 + 1).ID] / 35;
	    }else {
	    	rt2SecondSegmentDuration = rt2.duration - rt2FirstSegmentDuration;
	    }
	    
	    if (rt1FirstSegmentLoad +  rt2SecondSegmentLoad > rt1.capacity)
	    {
	        return true;
	    }
	    
	    if (rt2FirstSegmentLoad +  rt1SecondSegmentLoad > rt2.capacity)
	    {
	        return true;
	    }
	    
	    //duration constraints for first route 
	    if (nodeInd2 !=rt2.nodes.size() - 1 ) {
	    	if (rt1FirstSegmentDuration + rt2SecondSegmentDuration
	    			+ distanceMatrix[rt1.nodes.get(nodeInd1).ID][rt2.nodes.get(nodeInd2 + 1).ID] / 35 > 3.5) {
	    	return true;
	    	}
	    }else {
	    	if (rt1FirstSegmentDuration + rt2SecondSegmentDuration > 3.5) {
	    		return true;
	    	}
	    }
	  //duration constraints for second route 
	    if (nodeInd1 !=rt1.nodes.size() - 1 ) {
	    	if (rt2FirstSegmentDuration + rt1SecondSegmentDuration
	    			+ distanceMatrix[rt2.nodes.get(nodeInd2).ID][rt1.nodes.get(nodeInd1 + 1).ID] / 35 > 3.5) {
	    	return true;
	    	}
	    }else {
	    	if (rt2FirstSegmentDuration + rt1SecondSegmentDuration > 3.5) {
	    		return true;
	    	}
	    }
	    return false;
    }
	
}
