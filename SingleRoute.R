# Load required library
library("lpSolve")

# Function to calculate the Euclidean distance matrix between points
Distmat = function(X){
  # Get the number of points
  d = length(X[,1])
  # Calculate the distance matrix using the formula
  D = sqrt((X[,1]%*%t(rep(1,d)) - rep(1,d)%*%t(X[,1]))^2 + (X[,2]%*%t(rep(1,d)) - rep(1,d)%*%t(X[,2]))^2)
  # Return the distance matrix
  return(D)
}

# Function to find the shortest route that visits all points in the given order
SingleRoute = function(dist,X=NULL,plot=TRUE){
  # Get the number of points
  d = length(dist[,1])
  # Extract the distance values (excluding the diagonal and the last row and column)
  c = dist[-c((d+1)*(1:d)-d)]
  # Set the constraint matrix (each point must be visited exactly once)
  b = rep(1,2*d)
  A = matrix(0,2*d,length(c))
  I = (1:(d*d))[-c((d+1)*(1:d)-d)]
  for(i in 1:length(I)){
    A[((I[i]-1)%%d+1),i] = 1
    A[((I[i]-1)%/%d+1)+d,i] = 1
  }
  I2 = I
  w = TRUE;  lap = 0
  while (w){
    w = FALSE
    # Solve the linear programming problem
    res = lp(objective.in = c,const.mat = A,const.dir = c(rep("=",2*d),rep("<=",lap)),const.rhs = b,int.vec = 1:length(c))
    # Get the indices of the chosen edges
    e = res$solution;  I = I2[(round(e)==1)];  totaldist = sum(dist[I])
    # Construct the route based on the chosen edges
    route = c(((I[1]-1)%%d+1),((I[1]-1)%/%d+1),numeric(d-2))
    #return(I)
    
    # Print the order in which the points are visited
    #print(rbind(((I-1)%%d+1),((I-1)%/%d+1)))
    for(i in 3:d){
      r = route[i-1]
      for(j in 1:d){
        if(r == ((I[j]-1)%%d+1)){
          route[i] = ((I[j]-1)%/%d+1)
          break
        }
      }
      if(route[1] == route[i]){
        J = route[1:(i-1)];  a = rep(0,length(c))
        b = c(b,length(J)-1)
        lap = lap + 1
        for(j in 1:length(J)){
          a = a + ((((I2-1)%%d+1) == J[j]) + (((I2-1)%/%d+1) == J[j]))
        }
        a = 1*(a==2)
        A = rbind(A,a)
        w = TRUE
        break
      }
    }
  }
  if(is.null(X)){ plot = FALSE }
  # Plot the points and the shortest route (if points are given)
  if(plot){
    X = data.frame(X,row.names = 1:d)
    colnames(X) = c("x","y")
    g <- ggplot(X[c(route,route[1]), ], aes_string(x = names(X)[1], y = names(X)[2])) + geom_path(lineend = "round", linetype = 2, show.legend = TRUE) + labs(title = "Solution with Linear Programming")
    g <- g + annotate("text", x = X[, 1], y = X[, 2], label = 1:d)
    plot(g)
  }
  # Return the shortest route
  return(route)
}

ClusterSingleRoute = function(Locations,Cluster,title="Cluster LP",plot=TRUE){
  X = Locations[,-1];  C = Cluster;  k = max(Cluster)
  Routes = list();  d = length(X[,1])
  for(i in 1:k){
    I = c(TRUE,C==i)
    D = Distmat(X[I,])
    route = SingleRoute(D,X[I,],FALSE)
    route = ((1:d)[I])[route];  j = match(1,route)
    if(j > 1){ route = c(route[j:length(route)],route[1:(j-1)]) }
    Routes = append(Routes,list(c(route,route[1])))
  }
  if(plot){
    g <- ggplot(X[unlist(Routes), ], aes_string(x = names(X)[1], y = names(X)[2])) + geom_path(lineend = "round", linetype = 2, show.legend = TRUE) + labs(title = title) 
    g <- g + annotate("text", x = X[, 1], y = X[, 2], label = 1:d, col=c(1,C+1))
    plot(g)
  }
  #print(Routes)
  return(Routes)
}

Cluster = kmeans(locations[-1,2:3],6)[[1]]
Cluster
locations
Total_Cost(ClusterSingleRoute(locations,Cluster,title="Linear Programming with Clustering"),DMat)
SingleRoute(locations)

# Example usage
# Create a distance matrix for a set of random points
X = matrix(runif(50),25,2)
DM2 = Distmat(X)
# Find the shortest route that visits all points in the given order
DFDVV = SingleRoute(DM2,X)

cbind((DFDVV-1)%%25+1,(DFDVV-1)%/%25+1)


#Function used to plot sub-routes, generated by linear programming

DDDDD = function(){
  X = data.frame(X,row.names = 1:25)
  colnames(X) = c("x","y")
  g <- ggplot() + geom_path(aes(x=X[c(1,4,12,1),1],y=X[c(1,4,12,1),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(2,11,2),1],y=X[c(2,11,2),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(3,13,3),1],y=X[c(3,13,3),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(5,24,8,5),1],y=X[c(5,24,8,5),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(6,25,6),1],y=X[c(6,25,6),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(7,9,7),1],y=X[c(7,9,7),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(10,18,10),1],y=X[c(10,18,10),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(14,17,14),1],y=X[c(14,17,14),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(15,23,15),1],y=X[c(15,23,15),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(16,20,19,16),1],y=X[c(16,20,19,16),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + geom_path(aes(x=X[c(21,22,21),1],y=X[c(21,22,21),2]),lineend = "round", linetype = 2, show.legend = TRUE)
  g <- g + annotate("text", x = X[, 1], y = X[, 2], label = 1:25) + xlab("x") + ylab("y") + labs(title = "Subroutes Created by Linear Programming")
  plot(g)
}


X2
X2 = cbind(1:25,X)
rownames(X2) = 1:25;  colnames(X2) = c("ID","X","Y")
X2 = as.data.frame(X2)
ClusterSingleRoute(X2,kmeans(X2[-1,2:3],2)[[1]],title="Linear Programming with Clustering")
