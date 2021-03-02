package vrp;

import java.util.ArrayList;

public class Route {
        
    ArrayList <Node> nodes = new ArrayList();
        double distance;
        double load;
        double capacity;
        double duration;
        
        public Route(double cap) 
        {
            distance = 0;
            nodes = new ArrayList();
            load = 0;
            capacity = cap;
            duration=0;
        }
}
