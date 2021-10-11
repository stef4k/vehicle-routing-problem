package vrp;

public class mainClass {

    
    /**
     * @param args the command line arguments
     * @authors stef4k & AdamAvrantinis
     */
    public static void main(String[] args) {
        
        VRP vrp = new VRP();
        vrp.GenerateNetworkRandomly();
        Solution s1 =vrp.Solve();
        System.out.println(" and the cost of the objective function is: " +calculate_distance_solution(s1) +" km"); 
        
        VND vnd = new VND();
        vnd.GenerateNetworkRandomly();
        Solution s2 =vnd.solve();
        System.out.println(" and the cost of the objective function is: " +calculate_distance_solution(s2) +" km");
    }
    
    /**
     * Method that calculates the value of the objective function for any solution. Adds up all the 
     * distances of the routes and returns a double
     */
    public static double calculate_distance_solution(Solution s) {
		double distance = 0;
		for (Route route : s.routes) {
			distance += route.distance;
		}
		return distance;
	}
    
}
